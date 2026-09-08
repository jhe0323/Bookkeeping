#!/usr/bin/env python3
"""Generate reproducible Monte Carlo parameter realizations.

The output Parquet file contains one row per realization with a nested override
mapping serialized as JSON. These overrides are consumed by ``src.mc_runner``.
This script does not open the large LULCC NetCDF inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Dict

import numpy as np
import pandas as pd
import yaml

CODE_DIR = Path(__file__).resolve().parent.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from src.mc_sampling import (  # noqa: E402
    canonical_json,
    deep_merge,
    sample_overrides,
    validate_model_parameters,
)


def _load_yaml(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("Top-level YAML must be a mapping: {}".format(path))
    return data


def _resolve(repo_root: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = repo_root / path
    return path.resolve()


def _effective_baseline(repo_root: Path, run_config_path: Path) -> Dict[str, Any]:
    run_cfg = _load_yaml(run_config_path)
    parameter_path = _resolve(
        repo_root,
        str((run_cfg.get("parameters", {}) or {}).get("path", "config/config.yml")),
    )
    experiment_path = _resolve(
        repo_root,
        str((run_cfg.get("experiment", {}) or {}).get(
            "path", "config/experiments/harvest_area.yml"
        )),
    )

    baseline = _load_yaml(parameter_path)
    if experiment_path.is_file():
        baseline = deep_merge(baseline, _load_yaml(experiment_path))
    baseline = deep_merge(baseline, run_cfg.get("model_overrides", {}) or {})
    validate_model_parameters(baseline)
    return baseline


def _has_active_sampling(spec: Dict[str, Any]) -> bool:
    for rule in spec.get("parameters", []) or []:
        if bool(rule.get("enabled", True)):
            return True
    for group in spec.get("simplex_groups", []) or []:
        if bool(group.get("enabled", True)):
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--spec",
        required=True,
        help="Monte Carlo sampling specification YAML.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing sample table and manifest.",
    )
    parser.add_argument(
        "--allow-unconfirmed-ranges",
        action="store_true",
        help="Allow a spec whose scientific_ranges_confirmed flag is false.",
    )
    args = parser.parse_args()

    repo_root = CODE_DIR.resolve()
    spec_path = _resolve(repo_root, args.spec)
    spec_text = spec_path.read_text(encoding="utf-8")
    spec = _load_yaml(spec_path)

    mc = spec.get("monte_carlo", {}) or {}
    n_samples = int(mc.get("n_samples", 500))
    if n_samples <= 0:
        raise ValueError("monte_carlo.n_samples must be > 0")
    master_seed = int(mc.get("seed", 20260831))
    include_baseline = bool(mc.get("include_baseline", True))

    confirmed = bool(mc.get("scientific_ranges_confirmed", False))
    if not confirmed and not args.allow_unconfirmed_ranges:
        raise RuntimeError(
            "The sampling spec has scientific_ranges_confirmed=false. "
            "Confirm/replace the uncertainty ranges first, or use "
            "--allow-unconfirmed-ranges only for a pipeline smoke test."
        )
    if not _has_active_sampling(spec):
        raise RuntimeError(
            "No active Monte Carlo parameter or simplex rules were configured."
        )

    run_config_path = _resolve(
        repo_root,
        str(mc.get("run_config", "config/run_mc_1deg.yml")),
    )
    baseline = _effective_baseline(repo_root, run_config_path)

    dynamic_enabled = bool(
        (baseline.get("dynamic_carbon_density", {}) or {}).get("enabled", False)
    )
    if dynamic_enabled:
        target_text = canonical_json({
            "parameters": spec.get("parameters", []),
            "simplex_groups": spec.get("simplex_groups", []),
        })
        if "carbon_density" in target_text:
            raise RuntimeError(
                "Static carbon_density is being sampled while dynamic carbon "
                "density is enabled. Disable dynamic density for this MC design "
                "or implement a transient-density multiplier experiment separately."
            )

    output_path = _resolve(
        repo_root,
        str(mc.get("output", "config/mc_samples.parquet")),
    )
    manifest_path = output_path.with_suffix(output_path.suffix + ".manifest.json")
    if (output_path.exists() or manifest_path.exists()) and not args.force:
        raise FileExistsError(
            "Sample output already exists. Use --force to replace it: {}".format(
                output_path
            )
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seed_sequence = np.random.SeedSequence(master_seed)
    child_sequences = seed_sequence.spawn(n_samples)

    rows = []
    if include_baseline:
        rows.append({
            "sample_id": 0,
            "sample_kind": "baseline",
            "sample_seed": -1,
            "override_sha256": hashlib.sha256(b"{}").hexdigest(),
            "overrides_json": "{}",
            "draws_json": "{}",
        })
        first_mc_id = 1
    else:
        first_mc_id = 0

    for offset, child in enumerate(child_sequences):
        sample_id = first_mc_id + offset
        sample_seed = int(child.generate_state(1, dtype=np.uint32)[0])
        rng = np.random.default_rng(sample_seed)
        overrides, draws = sample_overrides(
            baseline_config=baseline,
            rng=rng,
            spec=spec,
        )
        overrides_json = canonical_json(overrides)
        rows.append({
            "sample_id": sample_id,
            "sample_kind": "mc",
            "sample_seed": sample_seed,
            "override_sha256": hashlib.sha256(
                overrides_json.encode("utf-8")
            ).hexdigest(),
            "overrides_json": overrides_json,
            "draws_json": canonical_json(draws),
        })

    frame = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    frame.to_parquet(output_path, index=False, engine="pyarrow")

    manifest = {
        "schema_version": 1,
        "sampling_spec": str(spec_path),
        "sampling_spec_sha256": hashlib.sha256(
            spec_text.encode("utf-8")
        ).hexdigest(),
        "run_config": str(run_config_path),
        "n_random_samples": n_samples,
        "include_baseline": include_baseline,
        "n_rows": int(len(frame)),
        "master_seed": master_seed,
        "scientific_ranges_confirmed": confirmed,
        "sample_table": str(output_path),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )

    print("[OK] samples: {}".format(output_path))
    print("[OK] manifest: {}".format(manifest_path))
    print("[INFO] rows={}, random={}, baseline={}".format(
        len(frame), n_samples, include_baseline
    ))
    print("[INFO] sample_id range: {}..{}".format(
        int(frame.sample_id.min()), int(frame.sample_id.max())
    ))


if __name__ == "__main__":
    main()
