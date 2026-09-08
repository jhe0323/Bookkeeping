#!/usr/bin/env python3
"""Evaluate Monte Carlo ensemble convergence for Net_Emissions.

The largest completed ensemble is used only as the reference distribution.
It is NOT allowed to "pass" against itself when determining the recommended
minimum ensemble size.

For each candidate N < N_total, repeated random subsets are compared with the
full completed ensemble for:
- annual mean / P05 / P50 / P95 trajectories;
- cumulative mean / SD / P05 / P50 / P95 over a selected period.

If no candidate N smaller than the full ensemble passes, convergence is reported
as NOT DEMONSTRATED. The correct response is to add more Monte Carlo samples and
repeat the test, rather than declaring N_total converged because it matches
itself exactly.

Python 3.8 compatible.
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


STAT_NAMES = (
    "mean",
    "sd",
    "p05",
    "p50",
    "p95",
)


def _stats(
    values: np.ndarray,
) -> Dict[str, np.ndarray]:
    return {
        "mean": np.mean(values, axis=0),
        "sd": (
            np.std(values, axis=0, ddof=1)
            if values.shape[0] > 1
            else np.zeros(
                values.shape[1:] or ()
            )
        ),
        "p05": np.quantile(
            values, 0.05, axis=0
        ),
        "p50": np.quantile(
            values, 0.50, axis=0
        ),
        "p95": np.quantile(
            values, 0.95, axis=0
        ),
    }


def _safe_relative_error(
    candidate,
    reference,
    fallback_scale,
) -> float:
    scale = abs(float(reference))

    if scale < 1.0e-12:
        scale = max(
            abs(float(fallback_scale)),
            1.0e-12,
        )

    return (
        100.0
        * abs(
            float(candidate)
            - float(reference)
        )
        / scale
    )


def _parse_sizes(
    text: str,
    n_total: int,
) -> List[int]:
    sizes = []

    for token in text.split(","):
        token = token.strip()
        if not token:
            continue

        value = int(token)
        if value < 2:
            raise ValueError(
                "All convergence sizes must be >= 2"
            )

        if value < n_total:
            sizes.append(value)

    return sorted(set(sizes))


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        default="config/run_mc_1deg.yml",
    )
    parser.add_argument(
        "--sizes",
        default="25,50,100,200,300,400",
        help=(
            "Candidate ensemble sizes. Values >= the "
            "completed ensemble size are ignored because "
            "the full ensemble is the reference."
        ),
    )
    parser.add_argument(
        "--start-year",
        type=int,
        default=1850,
    )
    parser.add_argument(
        "--end-year",
        type=int,
        default=2020,
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260908,
    )
    parser.add_argument(
        "--mean-threshold-pct",
        type=float,
        default=2.0,
    )
    parser.add_argument(
        "--spread-threshold-pct",
        type=float,
        default=5.0,
    )
    parser.add_argument(
        "--annual-nrmse-threshold-pct",
        type=float,
        default=5.0,
    )

    args = parser.parse_args()

    if args.replicates < 1:
        raise ValueError(
            "--replicates must be >= 1"
        )

    config = load_run_config(
        args.config
    )
    out_root = mc_output_root(
        config
    )

    matrix_path = (
        out_root
        / "net_emissions_samples.parquet"
    )
    meta_path = (
        out_root
        / "ensemble_summary.json"
    )

    if (
        not matrix_path.is_file()
        or not meta_path.is_file()
    ):
        raise FileNotFoundError(
            "Run tools.summarize_mc first; "
            "missing convergence input."
        )

    frame = pd.read_parquet(
        matrix_path,
        engine="pyarrow",
    )

    years = frame[
        "year"
    ].to_numpy(dtype=int)

    sample_cols = [
        column
        for column in frame.columns
        if column.startswith("sample_")
    ]

    if len(sample_cols) < 3:
        raise RuntimeError(
            "Need at least three MC samples "
            "for convergence analysis."
        )

    mask = (
        (years >= int(args.start_year))
        & (years <= int(args.end_year))
    )

    if not np.any(mask):
        raise ValueError(
            "No years within {}-{}".format(
                args.start_year,
                args.end_year,
            )
        )

    selected_years = years[mask]

    values = frame.loc[
        mask,
        sample_cols,
    ].to_numpy(
        dtype=np.float64
    ).T

    n_total = values.shape[0]

    candidate_sizes = _parse_sizes(
        args.sizes,
        n_total,
    )

    if not candidate_sizes:
        raise ValueError(
            "No candidate size is smaller than "
            "the completed ensemble size N={}. "
            "Provide smaller --sizes values.".format(
                n_total
            )
        )

    rng = np.random.default_rng(
        args.seed
    )

    # Full completed ensemble = reference only.
    full_annual = _stats(
        values
    )

    cumulative = values.sum(
        axis=1
    )

    full_cum_raw = _stats(
        cumulative[:, None]
    )

    full_cum = {
        key: float(
            np.asarray(value).reshape(-1)[0]
        )
        for key, value
        in full_cum_raw.items()
    }

    # Annual normalization scale:
    # average full-ensemble P05-P95 width.
    annual_scale = float(
        np.mean(
            np.abs(
                full_annual["p95"]
                - full_annual["p05"]
            )
        )
    )

    if annual_scale < 1.0e-12:
        annual_scale = float(
            np.mean(
                np.abs(
                    full_annual["mean"]
                )
            )
        )

    annual_scale = max(
        annual_scale,
        1.0e-12,
    )

    # Fallback scale for cumulative statistics
    # whose full-reference value is near zero.
    cumulative_scale = max(
        abs(full_cum["mean"]),
        abs(
            full_cum["p95"]
            - full_cum["p05"]
        ),
        1.0e-12,
    )

    rows = []

    for n in candidate_sizes:
        rep_metrics = []

        for _ in range(
            int(args.replicates)
        ):
            idx = rng.choice(
                n_total,
                size=n,
                replace=False,
            )

            subset = values[
                idx,
                :,
            ]

            annual = _stats(
                subset
            )

            cum_raw = _stats(
                subset.sum(
                    axis=1
                )[:, None]
            )

            cum = {
                key: float(
                    np.asarray(value)
                    .reshape(-1)[0]
                )
                for key, value
                in cum_raw.items()
            }

            metric = {}

            for stat in (
                "mean",
                "p05",
                "p50",
                "p95",
            ):
                rmse = float(
                    np.sqrt(
                        np.mean(
                            (
                                annual[stat]
                                - full_annual[stat]
                            )
                            ** 2
                        )
                    )
                )

                metric[
                    "annual_{}_nrmse_pct".format(
                        stat
                    )
                ] = (
                    100.0
                    * rmse
                    / annual_scale
                )

            for stat in STAT_NAMES:
                metric[
                    "cumulative_{}_error_pct".format(
                        stat
                    )
                ] = _safe_relative_error(
                    cum[stat],
                    full_cum[stat],
                    cumulative_scale,
                )

            rep_metrics.append(
                metric
            )

        summary = {
            "n": int(n),
            "replicates": int(
                args.replicates
            ),
            "is_reference": False,
        }

        for key in sorted(
            rep_metrics[0]
        ):
            array = np.asarray(
                [
                    metric[key]
                    for metric
                    in rep_metrics
                ],
                dtype=float,
            )

            summary[
                key + "_median"
            ] = float(
                np.median(array)
            )

            summary[
                key + "_p95"
            ] = float(
                np.quantile(
                    array,
                    0.95,
                )
            )

        pass_mean = (
            summary[
                "cumulative_mean_error_pct_p95"
            ]
            <= args.mean_threshold_pct
        )

        pass_spread = all(
            summary[
                "cumulative_{}_error_pct_p95".format(
                    stat
                )
            ]
            <= args.spread_threshold_pct
            for stat in (
                "sd",
                "p05",
                "p95",
            )
        )

        pass_annual = all(
            summary[
                "annual_{}_nrmse_pct_p95".format(
                    stat
                )
            ]
            <= args.annual_nrmse_threshold_pct
            for stat in (
                "mean",
                "p05",
                "p50",
                "p95",
            )
        )

        summary["pass"] = bool(
            pass_mean
            and pass_spread
            and pass_annual
        )

        rows.append(
            summary
        )

    # Add the full ensemble as an explicitly labelled reference row.
    # It is NOT evaluated for convergence and cannot become recommended_minimum_n.
    reference_row = {
        "n": int(n_total),
        "replicates": 0,
        "is_reference": True,
        "pass": False,
    }

    metric_names = [
        key
        for key in rows[0].keys()
        if key.endswith("_median")
        or key.endswith("_p95")
    ]

    for key in metric_names:
        reference_row[key] = 0.0

    rows.append(
        reference_row
    )

    result = pd.DataFrame(
        rows
    ).sort_values(
        "n"
    ).reset_index(
        drop=True
    )

    passing = result.loc[
        (result["is_reference"] == False)  # noqa: E712
        & (result["pass"] == True),  # noqa: E712
        "n",
    ].tolist()

    recommended_n = (
        int(min(passing))
        if passing
        else None
    )

    convergence_demonstrated = (
        recommended_n is not None
    )

    parquet_path = (
        out_root
        / "convergence_summary.parquet"
    )
    csv_path = (
        out_root
        / "convergence_summary.csv"
    )
    json_path = (
        out_root
        / "convergence_summary.json"
    )

    result.to_parquet(
        parquet_path,
        index=False,
        engine="pyarrow",
    )
    result.to_csv(
        csv_path,
        index=False,
    )

    metadata = {
        "variable": "Net_Emissions",
        "period": [
            int(selected_years[0]),
            int(selected_years[-1]),
        ],
        "n_total_reference": int(
            n_total
        ),
        "candidate_sizes": candidate_sizes,
        "replicates": int(
            args.replicates
        ),
        "random_seed": int(
            args.seed
        ),
        "reference": (
            "full completed ensemble; "
            "reference is never allowed "
            "to pass against itself"
        ),
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
        "convergence_demonstrated": bool(
            convergence_demonstrated
        ),
        "recommended_minimum_n": (
            recommended_n
        ),
        "interpretation": (
            "If convergence_demonstrated is false, "
            "increase the Monte Carlo ensemble and "
            "repeat the analysis. Do not interpret "
            "N_total as converged merely because it "
            "matches the reference ensemble."
        ),
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

    print(
        result.to_string(
            index=False
        )
    )

    if convergence_demonstrated:
        print(
            "[RESULT] convergence demonstrated; "
            "recommended_minimum_n={}".format(
                recommended_n
            )
        )
    else:
        print(
            "[RESULT] convergence NOT demonstrated "
            "within tested candidate sizes. "
            "Increase ensemble size and repeat."
        )

    print(
        "[OK] {}".format(
            parquet_path
        )
    )
    print(
        "[OK] {}".format(
            csv_path
        )
    )
    print(
        "[OK] {}".format(
            json_path
        )
    )


if __name__ == "__main__":
    main()
