#!/usr/bin/env python3
"""Evaluate Monte Carlo ensemble convergence for Net_Emissions.

The largest completed, current ensemble is the reference distribution only; it
cannot pass against itself. Candidate N values are evaluated by repeated random
subsets without replacement. A recommended minimum N is reported only when that
candidate AND every larger tested candidate pass the thresholds, preventing a
single non-monotonic lucky PASS at small N from being selected.

The script also validates that ``ensemble_summary.json`` and
``net_emissions_samples.parquet`` still correspond to the current MC config,
parameters, inputs, sample table, and result-affecting code. Run
``tools.summarize_mc`` again if this check fails.

Python 3.8 compatible.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Dict, List

import numpy as np
import pandas as pd

CODE_DIR = Path(__file__).resolve().parent.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from src.mc_runner import (  # noqa: E402
    current_mc_identity,
    mc_output_root,
    mc_sample_table_path,
)
from src.mc_sampling import canonical_json  # noqa: E402
from src.run_config import load_run_config, sha256_file  # noqa: E402


STAT_NAMES = ("mean", "sd", "p05", "p50", "p95")


def _same(a: Any, b: Any) -> bool:
    return canonical_json(a) == canonical_json(b)


def _stats(values: np.ndarray) -> Dict[str, np.ndarray]:
    return {
        "mean": np.mean(values, axis=0),
        "sd": (
            np.std(values, axis=0, ddof=1)
            if values.shape[0] > 1
            else np.zeros(values.shape[1:] or ())
        ),
        "p05": np.quantile(values, 0.05, axis=0),
        "p50": np.quantile(values, 0.50, axis=0),
        "p95": np.quantile(values, 0.95, axis=0),
    }


def _safe_relative_error(candidate, reference, fallback_scale) -> float:
    scale = abs(float(reference))
    if scale < 1.0e-12:
        scale = max(abs(float(fallback_scale)), 1.0e-12)
    return 100.0 * abs(float(candidate) - float(reference)) / scale


def _parse_sizes(text: str, n_total: int) -> List[int]:
    sizes = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        value = int(token)
        if value < 2:
            raise ValueError("All convergence sizes must be >= 2")
        if value < n_total:
            sizes.append(value)
    return sorted(set(sizes))


def _validate_summary_identity(config, meta: Dict[str, Any], frame: pd.DataFrame) -> None:
    sample_table = mc_sample_table_path(config)
    if not sample_table.is_file():
        raise FileNotFoundError("MC sample table not found: {}".format(sample_table))

    table_sha = sha256_file(sample_table)
    if str(meta.get("sample_table_sha256", "")) != table_sha:
        raise RuntimeError(
            "Ensemble summary is stale: sample_table_sha256 changed. "
            "Run tools.summarize_mc again."
        )

    sample_ids = [int(v) for v in meta.get("sample_ids", [])]
    if not sample_ids:
        raise RuntimeError("ensemble_summary.json has no sample_ids")

    sample_info = pd.read_parquet(sample_table, engine="pyarrow")
    first = sample_info.loc[
        sample_info["sample_id"].astype(int) == sample_ids[0]
    ]
    if len(first) != 1:
        raise RuntimeError(
            "Cannot resolve first summarized sample {} in current sample table"
            .format(sample_ids[0])
        )

    override_sha = str(first.iloc[0].get("override_sha256", ""))
    expected = current_mc_identity(config, table_sha, override_sha)
    expected.pop("override_sha256", None)
    for key, value in expected.items():
        if not _same(meta.get(key), value):
            raise RuntimeError(
                "Ensemble summary is stale: {!r} changed. "
                "Run tools.merge_all_mc and tools.summarize_mc again."
                .format(key)
            )

    sample_columns = [
        column for column in frame.columns if column.startswith("sample_")
    ]
    expected_columns = ["sample_{:06d}".format(v) for v in sample_ids]
    if sample_columns != expected_columns:
        raise RuntimeError(
            "net_emissions_samples.parquet columns do not match ensemble_summary.json"
        )
    if int(meta.get("n_ensemble_samples", -1)) != len(sample_ids):
        raise RuntimeError("ensemble_summary.json sample count is inconsistent")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/run_mc_1deg.yml")
    parser.add_argument(
        "--sizes",
        default="25,50,100,200,300,400",
        help=(
            "Candidate ensemble sizes. Values >= the completed ensemble size "
            "are ignored because the full ensemble is the reference."
        ),
    )
    parser.add_argument("--start-year", type=int, default=1850)
    parser.add_argument("--end-year", type=int, default=2020)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--mean-threshold-pct", type=float, default=2.0)
    parser.add_argument("--spread-threshold-pct", type=float, default=5.0)
    parser.add_argument("--annual-nrmse-threshold-pct", type=float, default=5.0)
    args = parser.parse_args()

    if args.replicates < 1:
        raise ValueError("--replicates must be >= 1")

    config = load_run_config(args.config)
    out_root = mc_output_root(config)
    matrix_path = out_root / "net_emissions_samples.parquet"
    meta_path = out_root / "ensemble_summary.json"

    if not matrix_path.is_file() or not meta_path.is_file():
        raise FileNotFoundError(
            "Run tools.summarize_mc first; convergence inputs are missing."
        )

    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    frame = pd.read_parquet(matrix_path, engine="pyarrow")
    _validate_summary_identity(config, meta, frame)

    years = frame["year"].to_numpy(dtype=int)
    sample_cols = [
        column for column in frame.columns if column.startswith("sample_")
    ]
    if len(sample_cols) < 3:
        raise RuntimeError("Need at least three MC samples for convergence analysis.")

    mask = (years >= int(args.start_year)) & (years <= int(args.end_year))
    if not np.any(mask):
        raise ValueError(
            "No years within {}-{}".format(args.start_year, args.end_year)
        )

    selected_years = years[mask]
    values = frame.loc[mask, sample_cols].to_numpy(dtype=np.float64).T
    n_total = values.shape[0]
    candidate_sizes = _parse_sizes(args.sizes, n_total)
    if not candidate_sizes:
        raise ValueError(
            "No candidate size is smaller than completed ensemble N={}. "
            "Provide smaller --sizes values or add more MC samples."
            .format(n_total)
        )

    rng = np.random.default_rng(args.seed)
    full_annual = _stats(values)
    cumulative = values.sum(axis=1)
    full_cum_raw = _stats(cumulative[:, None])
    full_cum = {
        key: float(np.asarray(value).reshape(-1)[0])
        for key, value in full_cum_raw.items()
    }

    annual_scale = float(
        np.mean(np.abs(full_annual["p95"] - full_annual["p05"]))
    )
    if annual_scale < 1.0e-12:
        annual_scale = float(np.mean(np.abs(full_annual["mean"])))
    annual_scale = max(annual_scale, 1.0e-12)

    cumulative_scale = max(
        abs(full_cum["mean"]),
        abs(full_cum["p95"] - full_cum["p05"]),
        1.0e-12,
    )

    rows = []
    for n in candidate_sizes:
        rep_metrics = []
        for _ in range(int(args.replicates)):
            idx = rng.choice(n_total, size=n, replace=False)
            subset = values[idx, :]
            annual = _stats(subset)
            cum_raw = _stats(subset.sum(axis=1)[:, None])
            cum = {
                key: float(np.asarray(value).reshape(-1)[0])
                for key, value in cum_raw.items()
            }

            metric: Dict[str, float] = {}
            for stat in ("mean", "p05", "p50", "p95"):
                rmse = float(
                    np.sqrt(np.mean((annual[stat] - full_annual[stat]) ** 2))
                )
                metric["annual_{}_nrmse_pct".format(stat)] = (
                    100.0 * rmse / annual_scale
                )
            for stat in STAT_NAMES:
                metric["cumulative_{}_error_pct".format(stat)] = (
                    _safe_relative_error(
                        cum[stat], full_cum[stat], cumulative_scale
                    )
                )
            rep_metrics.append(metric)

        summary: Dict[str, Any] = {
            "n": int(n),
            "replicates": int(args.replicates),
            "is_reference": False,
        }
        for key in sorted(rep_metrics[0]):
            array = np.asarray([m[key] for m in rep_metrics], dtype=float)
            summary[key + "_median"] = float(np.median(array))
            summary[key + "_p95"] = float(np.quantile(array, 0.95))

        pass_mean = (
            summary["cumulative_mean_error_pct_p95"]
            <= args.mean_threshold_pct
        )
        pass_spread = all(
            summary["cumulative_{}_error_pct_p95".format(stat)]
            <= args.spread_threshold_pct
            for stat in ("sd", "p05", "p95")
        )
        pass_annual = all(
            summary["annual_{}_nrmse_pct_p95".format(stat)]
            <= args.annual_nrmse_threshold_pct
            for stat in ("mean", "p05", "p50", "p95")
        )
        summary["pass"] = bool(pass_mean and pass_spread and pass_annual)
        rows.append(summary)

    # A candidate is a stable pass only if it and ALL larger tested candidates pass.
    all_larger_pass = True
    for row in reversed(rows):
        all_larger_pass = bool(row["pass"]) and all_larger_pass
        row["stable_pass"] = bool(all_larger_pass)

    stable = [int(row["n"]) for row in rows if row["stable_pass"]]
    recommended_n = min(stable) if stable else None
    convergence_demonstrated = recommended_n is not None

    reference_row: Dict[str, Any] = {
        "n": int(n_total),
        "replicates": 0,
        "is_reference": True,
        "pass": False,
        "stable_pass": False,
    }
    metric_names = [
        key for key in rows[0]
        if key.endswith("_median") or key.endswith("_p95")
    ]
    for key in metric_names:
        reference_row[key] = 0.0

    result = pd.DataFrame(rows + [reference_row]).sort_values("n").reset_index(drop=True)
    parquet_path = out_root / "convergence_summary.parquet"
    csv_path = out_root / "convergence_summary.csv"
    json_path = out_root / "convergence_summary.json"
    result.to_parquet(parquet_path, index=False, engine="pyarrow")
    result.to_csv(csv_path, index=False)

    metadata = {
        "variable": "Net_Emissions",
        "period": [int(selected_years[0]), int(selected_years[-1])],
        "n_total_reference": int(n_total),
        "candidate_sizes": candidate_sizes,
        "replicates": int(args.replicates),
        "random_seed": int(args.seed),
        "reference": (
            "full completed current ensemble; reference is never allowed "
            "to pass against itself"
        ),
        "selection_rule": (
            "recommended_minimum_n is the smallest candidate for which that N "
            "and every larger tested candidate pass"
        ),
        "thresholds_pct": {
            "cumulative_mean": float(args.mean_threshold_pct),
            "cumulative_sd_p05_p95": float(args.spread_threshold_pct),
            "annual_nrmse": float(args.annual_nrmse_threshold_pct),
        },
        "convergence_demonstrated": bool(convergence_demonstrated),
        "recommended_minimum_n": recommended_n,
        "ensemble_summary_sample_table_sha256": meta.get("sample_table_sha256"),
    }
    json_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )

    print(result.to_string(index=False))
    if convergence_demonstrated:
        print(
            "[RESULT] convergence demonstrated; recommended_minimum_n={}"
            .format(recommended_n)
        )
    else:
        print(
            "[RESULT] convergence NOT demonstrated with a stable-pass suffix. "
            "Increase ensemble size and repeat."
        )
    print("[OK] {}".format(parquet_path))
    print("[OK] {}".format(csv_path))
    print("[OK] {}".format(json_path))


if __name__ == "__main__":
    main()
