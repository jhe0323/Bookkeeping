#!/usr/bin/env python3
"""Merge all longitude-band summaries for one Monte Carlo realization."""
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
    band_output_paths,
    current_mc_identity,
    global_output_path,
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


def _load_sample_row(sample_table: Path, sample_id: int) -> Dict[str, Any]:
    frame = pd.read_parquet(sample_table, engine="pyarrow")
    selected = frame.loc[frame["sample_id"].astype(int) == int(sample_id)]
    if len(selected) != 1:
        raise RuntimeError(
            "Could not uniquely resolve sample {} in {}".format(
                sample_id, sample_table
            )
        )
    return selected.iloc[0].to_dict()


def _global_is_current(
    out_path: Path,
    manifest_path: Path,
    *,
    sample_id: int,
    identity: Dict[str, Any],
) -> bool:
    if not out_path.is_file() or not manifest_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if int(manifest.get("sample_id", -1)) != int(sample_id):
            return False
        for key, expected in identity.items():
            if not _same(manifest.get(key), expected):
                return False
        frame = pd.read_parquet(out_path, engine="pyarrow")
        return "year" in frame.columns
    except Exception:
        return False


def merge_one(config_path: str, sample_id: int, force: bool = False) -> Path:
    config = load_run_config(config_path)

    bands_total_float = 360.0 / float(config.band_size_deg)
    bands_total = int(round(bands_total_float))
    if abs(bands_total_float - bands_total) > 1.0e-8:
        raise ValueError("server.band_size_deg must divide 360")

    sample_table = mc_sample_table_path(config)
    sample_table_sha = sha256_file(sample_table)
    sample_row = _load_sample_row(sample_table, sample_id)
    override_sha = str(sample_row.get("override_sha256", ""))
    identity = current_mc_identity(config, sample_table_sha, override_sha)

    out_path = global_output_path(config, sample_id)
    manifest_path = sample_manifest_path(config, sample_id)

    if (
        not force
        and _global_is_current(
            out_path,
            manifest_path,
            sample_id=sample_id,
            identity=identity,
        )
    ):
        print("[SKIP] current global output exists: {}".format(out_path))
        return out_path

    if out_path.exists():
        out_path.unlink()
    if manifest_path.exists():
        manifest_path.unlink()

    missing: List[str] = []
    metadata_list: List[Dict[str, Any]] = []
    frames: List[pd.DataFrame] = []

    for band_id in range(bands_total):
        parquet_path, metadata_path = band_output_paths(
            config, sample_id, band_id
        )
        if not parquet_path.is_file() or not metadata_path.is_file():
            missing.append("band_{:03d}".format(band_id))
            continue

        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if int(metadata.get("sample_id", -1)) != int(sample_id):
            raise RuntimeError("Sample ID mismatch in {}".format(metadata_path))
        if int(metadata.get("band_id", -1)) != band_id:
            raise RuntimeError("Band ID mismatch in {}".format(metadata_path))

        for key, expected in identity.items():
            if not _same(metadata.get(key), expected):
                raise RuntimeError(
                    "Stale MC band {}: mismatch for {!r}".format(
                        band_id, key
                    )
                )

        metadata_list.append(metadata)
        frames.append(pd.read_parquet(parquet_path, engine="pyarrow"))

    if missing:
        raise RuntimeError(
            "Cannot merge sample {}: missing {} of {} bands: {}".format(
                sample_id,
                len(missing),
                bands_total,
                ", ".join(missing[:20]),
            )
        )

    reference = metadata_list[0]
    identity_keys = (
        "override_sha256",
        "run_config_sha256",
        "model_code_sha256",
        "start_year",
        "end_year",
        "output_variables",
        "parameter_sha256",
        "experiment_sha256",
        "state_fingerprint",
        "transition_fingerprint",
        "pft_fingerprint",
    )
    for band_id, metadata in enumerate(metadata_list[1:], start=1):
        for key in identity_keys:
            if not _same(metadata.get(key), reference.get(key)):
                raise RuntimeError(
                    "Band metadata mismatch for {!r}: band 0 vs {}".format(
                        key, band_id
                    )
                )

    columns = list(frames[0].columns)
    if not columns or columns[0] != "year":
        raise ValueError("Band output must begin with a year column")

    years = frames[0]["year"].to_numpy(dtype=int)
    variables = columns[1:]
    total = np.zeros((len(years), len(variables)), dtype=np.float64)

    for band_id, frame in enumerate(frames):
        if list(frame.columns) != columns:
            raise RuntimeError("Column mismatch in band {}".format(band_id))
        if not np.array_equal(
            frame["year"].to_numpy(dtype=int), years
        ):
            raise RuntimeError("Year mismatch in band {}".format(band_id))
        total += frame[variables].to_numpy(dtype=np.float64)

    merged = pd.DataFrame(total, columns=variables)
    merged.insert(0, "year", years)
    _atomic_parquet(merged, out_path)

    manifest = {
        "schema_version": 2,
        "sample_id": int(sample_id),
        "sample_kind": str(sample_row.get("sample_kind", "mc")),
        "sample_seed": int(sample_row.get("sample_seed", -1)),
        "bands_total": int(bands_total),
        "output_variables": variables,
        "start_year": int(years[0]),
        "end_year": int(years[-1]),
        "global_output": str(out_path),
        "overrides_json": str(sample_row.get("overrides_json", "{}")),
        "draws_json": str(sample_row.get("draws_json", "{}")),
        "sample_table_sha256": sample_table_sha,
        "git_commit": reference.get("git_commit", "unknown"),
    }
    manifest.update(identity)
    _atomic_json(manifest, manifest_path)

    print("[OK] merged sample {} -> {}".format(sample_id, out_path))
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/run_mc_1deg.yml")
    parser.add_argument("--sample-id", type=int, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    merge_one(args.config, args.sample_id, force=args.force)


if __name__ == "__main__":
    main()
