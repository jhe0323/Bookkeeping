#!/usr/bin/env python3
"""Evaluate Monte Carlo ensemble convergence for Net_Emissions.

Repeated random subsets are compared with the full completed ensemble for both
annual trajectories and cumulative emissions over a selected period.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Dict, List

import numpy as np
import pandas as pd

CODE_DIR = Path(__file__).resolve().parent.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from src.mc_runner import mc_output_root  # noqa: E402
from src.run_config import load_run_config  # noqa: E402


STAT_NAMES = ("mean", "sd", "p05", "p50", "p95")


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


def _safe_relative_error(candidate, reference, fallback_scale):
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
        if value <= n_total:
            sizes.append(value)
    sizes.append(n_total)
    return sorted(set(sizes))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/run_mc_1deg.yml")
    parser.add_argument("--sizes", default="25,50,100,200,300,500")
    parser.add_argument("--start-year", type=int, default=1850)
    parser.add_argument("--end-year", type=int, default=2020)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--mean-threshold-pct", type=float, default=2.0)
    parser.add_argument("--spread-threshold-pct", type=float, default=5.0)
    parser.add_argument(
        "--annual-nrmse-threshold-pct",
        type=float,
        default=5.0,
    )
    args = parser.parse_args()

    config = load_run_config(args.config)
    out_root = mc_output_root(config)
    matrix_path = out_root / "net_emissions_samples.parquet"
    meta_path = out_root / "ensemble_summary.json"

    if not matrix_path.is_file() or not meta_path.is_file():
        raise FileNotFoundError(
            "Run tools.summarize_mc first; missing convergence input."
        )

    frame = pd.read_parquet(matrix_path, engine="pyarrow")
    years = frame["year"].to_numpy(dtype=int)
    sample_cols = [
        c for c in frame.columns if c.startswith("sample_")
    ]
    if len(sample_cols) < 2:
        raise RuntimeError(
            "Need at least two MC samples for convergence analysis."
        )

    mask = (
        (years >= int(args.start_year))
        & (years <= int(args.end_year))
    )
    if not np.any(mask):
        raise ValueError(
            "No years within {}-{}".format(
                args.start_year, args.end_year
            )
        )

    selected_years = years[mask]
    values = frame.loc[
        mask, sample_cols
    ].to_numpy(dtype=np.float64).T

    n_total = values.shape[0]
    sizes = _parse_sizes(args.sizes, n_total)
    rng = np.random.default_rng(args.seed)

    full_annual = _stats(values)
    cumulative = values.sum(axis=1)
    full_cum_raw = _stats(cumulative[:, None])
    full_cum = {
        k: float(np.asarray(v).reshape(-1)[0])
        for k, v in full_cum_raw.items()
    }

    annual_scale = float(
        np.mean(
            np.abs(
                full_annual["p95"] - full_annual["p05"]
            )
        )
    )
    if annual_scale < 1.0e-12:
        annual_scale = float(
            np.mean(np.abs(full_annual["mean"]))
        )
    annual_scale = max(annual_scale, 1.0e-12)

    cumulative_scale = max(
        abs(full_cum["mean"]),
        abs(full_cum["p95"] - full_cum["p05"]),
        1.0e-12,
    )

    rows = []
    for n in sizes:
        n_reps = 1 if n == n_total else int(args.replicates)
        rep_metrics = []

        for _ in range(n_reps):
            if n == n_total:
                idx = np.arange(n_total)
            else:
                idx = rng.choice(
                    n_total, size=n, replace=False
                )

            subset = values[idx, :]
            annual = _stats(subset)
            cum_raw = _stats(
                subset.sum(axis=1)[:, None]
            )
            cum = {
                k: float(np.asarray(v).reshape(-1)[0])
                for k, v in cum_raw.items()
            }

            metric = {}
            for stat in ("mean", "p05", "p50", "p95"):
                rmse = float(
                    np.sqrt(
                        np.mean(
                            (
                                annual[stat]
                                - full_annual[stat]
                            ) ** 2
                        )
                    )
                )
                metric[
                    "annual_{}_nrmse_pct".format(stat)
                ] = 100.0 * rmse / annual_scale

            for stat in STAT_NAMES:
                metric[
                    "cumulative_{}_error_pct".format(stat)
                ] = _safe_relative_error(
                    cum[stat],
                    full_cum[stat],
                    cumulative_scale,
                )

            rep_metrics.append(metric)

        summary = {
            "n": int(n),
            "replicates": int(n_reps),
        }
        for key in sorted(rep_metrics[0]):
            arr = np.asarray(
                [m[key] for m in rep_metrics],
                dtype=float,
            )
            summary[key + "_median"] = float(
                np.median(arr)
            )
            summary[key + "_p95"] = float(
                np.quantile(arr, 0.95)
            )

        pass_mean = (
            summary[
                "cumulative_mean_error_pct_p95"
            ]
            <= args.mean_threshold_pct
        )
        pass_spread = all(
            summary[
                "cumulative_{}_error_pct_p95".format(stat)
            ]
            <= args.spread_threshold_pct
            for stat in ("sd", "p05", "p95")
        )
        pass_annual = all(
            summary[
                "annual_{}_nrmse_pct_p95".format(stat)
            ]
            <= args.annual_nrmse_threshold_pct
            for stat in ("mean", "p05", "p50", "p95")
        )
        summary["pass"] = bool(
            pass_mean and pass_spread and pass_annual
        )
        rows.append(summary)

    result = pd.DataFrame(rows)
    parquet_path = out_root / "convergence_summary.parquet"
    csv_path = out_root / "convergence_summary.csv"
    json_path = out_root / "convergence_summary.json"

    result.to_parquet(
        parquet_path, index=False, engine="pyarrow"
    )
    result.to_csv(csv_path, index=False)

    passing = result.loc[
        result["pass"] == True, "n"  # noqa: E712
    ].tolist()
    recommended_n = (
        int(min(passing)) if passing else None
    )

    metadata = {
        "variable": "Net_Emissions",
        "period": [
            int(selected_years[0]),
            int(selected_years[-1]),
        ],
        "n_total": int(n_total),
        "sizes": sizes,
        "replicates": int(args.replicates),
        "random_seed": int(args.seed),
        "reference": "full completed ensemble",
        "thresholds_pct": {
            "cumulative_mean": float(
                args.mean_threshold_pct
            ),
            "cumulative_sd_p05_p95": float(
                args.spread_threshold_pct
            ),
            "annual_nrmse": float(
                args.annual_nrmse_threshold_pct
            ),
        },
        "recommended_minimum_n": recommended_n,
    }
    json_path.write_text(
        json.dumps(
            metadata,
            indent=2,
            ensure_ascii=False,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    print(result.to_string(index=False))
    print(
        "[RESULT] recommended_minimum_n={}".format(
            recommended_n
        )
    )


if __name__ == "__main__":
    main()
