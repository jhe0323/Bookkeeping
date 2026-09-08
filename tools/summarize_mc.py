#!/usr/bin/env python3
"""Calculate ensemble uncertainty statistics from current merged MC summaries."""
from __future__ import annotations

import argparse
import json
import os
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
    global_output_path,
    mc_output_root,
    mc_sample_table_path,
    sample_manifest_path,
)
from src.mc_sampling import canonical_json  # noqa: E402
from src.run_config import load_run_config, sha256_file  # noqa: E402


def _same(a: Any, b: Any) -> bool:
    return canonical_json(a) == canonical_json(b)


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(".{}.{}.tmp.parquet".format(path.stem, os.getpid()))
    try:
        frame.to_parquet(tmp, index=False, engine="pyarrow")
        os.replace(str(tmp), str(path))
    finally:
        if tmp.exists():
            tmp.unlink()


def _atomic_json(data: Dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(".{}.{}.tmp.json".format(path.stem, os.getpid()))
    try:
        tmp.write_text(
            json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True),
            encoding="utf-8",
        )
        os.replace(str(tmp), str(path))
    finally:
        if tmp.exists():
            tmp.unlink()


def _manifest_is_current(config, table_sha: str, row, path: Path):
    if not path.is_file():
        return False, "missing manifest"

    sample_id = int(row["sample_id"])
    override_sha = str(row.get("override_sha256", ""))
    expected = current_mc_identity(config, table_sha, override_sha)

    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return False, "manifest read error: {}".format(exc)

    if int(manifest.get("sample_id", -1)) != sample_id:
        return False, "sample_id mismatch"

    for key, value in expected.items():
        if not _same(manifest.get(key), value):
            return False, "{} mismatch".format(key)

    return True, ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/run_mc_1deg.yml")
    parser.add_argument(
        "--include-baseline-in-ensemble",
        action="store_true",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Skip missing/stale samples and report their IDs.",
    )
    args = parser.parse_args()

    config = load_run_config(args.config)
    sample_table = mc_sample_table_path(config)
    table_sha = sha256_file(sample_table)
    sample_info = pd.read_parquet(sample_table, engine="pyarrow")

    selected = []
    missing: List[int] = []
    stale: Dict[int, str] = {}

    for _, row in sample_info.sort_values("sample_id").iterrows():
        sample_id = int(row["sample_id"])
        kind = str(row.get("sample_kind", "mc"))

        if kind == "baseline" and not args.include_baseline_in_ensemble:
            continue

        global_path = global_output_path(config, sample_id)
        manifest_path = sample_manifest_path(config, sample_id)

        if not global_path.is_file():
            missing.append(sample_id)
            continue

        ok, reason = _manifest_is_current(
            config, table_sha, row, manifest_path
        )
        if not ok:
            stale[sample_id] = reason
            continue

        selected.append((sample_id, kind, global_path))

    if (missing or stale) and not args.allow_missing:
        raise RuntimeError(
            "Cannot summarize: missing={} stale={}; first missing={}, first stale={}".format(
                len(missing),
                len(stale),
                missing[:20],
                list(stale.items())[:10],
            )
        )
    if not selected:
        raise RuntimeError("No current merged MC outputs were found.")

    reference = pd.read_parquet(selected[0][2], engine="pyarrow")
    years = reference["year"].to_numpy(dtype=int)
    variables = list(reference.columns[1:])
    arrays = {name: [] for name in variables}
    sample_ids = []

    for sample_id, _, path in selected:
        frame = pd.read_parquet(path, engine="pyarrow")
        if list(frame.columns) != list(reference.columns):
            raise RuntimeError(
                "Column mismatch in sample {}".format(sample_id)
            )
        if not np.array_equal(
            frame["year"].to_numpy(dtype=int), years
        ):
            raise RuntimeError(
                "Year mismatch in sample {}".format(sample_id)
            )

        sample_ids.append(sample_id)
        for variable in variables:
            arrays[variable].append(
                frame[variable].to_numpy(dtype=np.float64)
            )

    summary_rows = []
    for variable in variables:
        values = np.stack(arrays[variable], axis=0)
        mean = np.mean(values, axis=0)
        sd = (
            np.std(values, axis=0, ddof=1)
            if values.shape[0] > 1
            else np.zeros(len(years))
        )
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

    if "Net_Emissions" in variables:
        matrix = np.stack(arrays["Net_Emissions"], axis=0).T
        net = pd.DataFrame(
            matrix,
            columns=[
                "sample_{:06d}".format(v) for v in sample_ids
            ],
        )
        net.insert(0, "year", years)
        _atomic_parquet(
            net, out_root / "net_emissions_samples.parquet"
        )

    first_row = sample_info.loc[
        sample_info["sample_id"].astype(int) == sample_ids[0]
    ].iloc[0]
    identity = current_mc_identity(
        config,
        table_sha,
        str(first_row.get("override_sha256", "")),
    )
    identity.pop("override_sha256", None)

    metadata = {
        "schema_version": 2,
        "n_ensemble_samples": len(sample_ids),
        "sample_ids": sample_ids,
        "missing_sample_ids": missing,
        "stale_sample_ids": stale,
        "baseline_included": bool(args.include_baseline_in_ensemble),
        "statistics": ["mean", "sd", "p05", "p50", "p95"],
        "sample_table_sha256": table_sha,
        "units": (
            "same as model output; carbon variables are t C, "
            "area variables are ha"
        ),
    }
    metadata.update(identity)
    _atomic_json(metadata, out_root / "ensemble_summary.json")

    print("[OK] ensemble summary: {}".format(summary_path))
    print(
        "[INFO] samples={}, missing={}, stale={}".format(
            len(sample_ids), len(missing), len(stale)
        )
    )


if __name__ == "__main__":
    main()
