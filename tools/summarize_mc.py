#!/usr/bin/env python3
"""Calculate ensemble uncertainty statistics from merged MC global summaries."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Dict, List

import numpy as np
import pandas as pd

CODE_DIR = Path(__file__).resolve().parent.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from src.mc_runner import global_output_path, mc_output_root, mc_sample_table_path  # noqa: E402
from src.run_config import load_run_config  # noqa: E402


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".{}.{}.tmp.parquet".format(path.stem, os.getpid()))
    try:
        frame.to_parquet(temporary, index=False, engine="pyarrow")
        os.replace(str(temporary), str(path))
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/run_mc_1deg.yml")
    parser.add_argument(
        "--include-baseline-in-ensemble",
        action="store_true",
        help="Normally the deterministic baseline is excluded from MC quantiles.",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Summarize completed samples even if some expected MC runs are missing.",
    )
    args = parser.parse_args()

    config = load_run_config(args.config, expected_resolution="1deg")
    sample_table = mc_sample_table_path(config)
    sample_info = pd.read_parquet(sample_table, engine="pyarrow")

    selected_rows = []
    missing: List[int] = []
    for _, row in sample_info.sort_values("sample_id").iterrows():
        sample_id = int(row["sample_id"])
        kind = str(row.get("sample_kind", "mc"))
        if kind == "baseline" and not args.include_baseline_in_ensemble:
            continue
        path = global_output_path(config, sample_id)
        if not path.is_file():
            missing.append(sample_id)
            continue
        selected_rows.append((sample_id, kind, path))

    if missing and not args.allow_missing:
        raise RuntimeError(
            "Missing merged global outputs for {} samples. First IDs: {}".format(
                len(missing), missing[:20]
            )
        )
    if not selected_rows:
        raise RuntimeError("No merged MC global outputs were found.")

    reference = pd.read_parquet(selected_rows[0][2], engine="pyarrow")
    years = reference["year"].to_numpy(dtype=int)
    variables = list(reference.columns[1:])
    arrays: Dict[str, List[np.ndarray]] = {name: [] for name in variables}
    sample_ids: List[int] = []

    for sample_id, _, path in selected_rows:
        frame = pd.read_parquet(path, engine="pyarrow")
        if list(frame.columns) != list(reference.columns):
            raise RuntimeError("Column mismatch in sample {}".format(sample_id))
        if not np.array_equal(frame["year"].to_numpy(dtype=int), years):
            raise RuntimeError("Year mismatch in sample {}".format(sample_id))
        sample_ids.append(sample_id)
        for variable in variables:
            arrays[variable].append(frame[variable].to_numpy(dtype=np.float64))

    summary_rows = []
    for variable in variables:
        values = np.stack(arrays[variable], axis=0)  # sample x year
        mean = np.mean(values, axis=0)
        sd = np.std(values, axis=0, ddof=1) if values.shape[0] > 1 else np.zeros(len(years))
        p05 = np.quantile(values, 0.05, axis=0)
        p50 = np.quantile(values, 0.50, axis=0)
        p95 = np.quantile(values, 0.95, axis=0)
        for idx, year in enumerate(years):
            summary_rows.append(
                {
                    "year": int(year),
                    "variable": variable,
                    "n": int(values.shape[0]),
                    "mean": float(mean[idx]),
                    "sd": float(sd[idx]),
                    "p05": float(p05[idx]),
                    "p50": float(p50[idx]),
                    "p95": float(p95[idx]),
                }
            )

    out_root = mc_output_root(config)
    summary_path = out_root / "ensemble_summary.parquet"
    _atomic_parquet(pd.DataFrame(summary_rows), summary_path)

    # A compact sample x year table for the primary ELUC diagnostic is useful for
    # convergence checks and plotting without reopening every sample directory.
    if "Net_Emissions" in variables:
        matrix = np.stack(arrays["Net_Emissions"], axis=0).T
        net_frame = pd.DataFrame(
            matrix,
            columns=["sample_{:06d}".format(value) for value in sample_ids],
        )
        net_frame.insert(0, "year", years)
        _atomic_parquet(net_frame, out_root / "net_emissions_samples.parquet")

    metadata = {
        "n_ensemble_samples": len(sample_ids),
        "sample_ids": sample_ids,
        "missing_sample_ids": missing,
        "baseline_included": bool(args.include_baseline_in_ensemble),
        "statistics": ["mean", "sd", "p05", "p50", "p95"],
        "units": "same as model output (carbon variables are t C; area variables are ha)",
    }
    (out_root / "ensemble_summary.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("[OK] ensemble summary: {}".format(summary_path))
    print("[INFO] ensemble samples={}, missing={}".format(len(sample_ids), len(missing)))


if __name__ == "__main__":
    main()
