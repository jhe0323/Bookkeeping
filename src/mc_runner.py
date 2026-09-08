"""Summary-only Monte Carlo runner for the Bookkeeping model.

Existing bands are reused only when sample, config, parameters, experiment,
large-input fingerprints, and result-affecting Python source all match the
current run. The code hash catches uncommitted local edits as well as commits.

Python 3.8 compatible.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from src.LULCCSimulator import FORMAL_OUTPUT_KEYS, STATE_OUTPUT_KEYS, LULCCSimulator
from src.mc_sampling import canonical_json, deep_merge
from src.run_config import build_run_metadata, load_run_config, sha256_file
from src.run_manager import _band_slice


DEFAULT_MC_OUTPUT_VARIABLES = [
    "Gross_Sources", "Gross_Sinks", "Net_Emissions",
    "Flux_Clearing", "Flux_Abandonment", "Flux_Harvest_Net", "Flux_Other",
    "Flux_Products_Total", "Flux_FD", "Flux_NFC", "Flux_FR", "Flux_NFR",
    "Flux_CAL", "Flux_WHp",
    "Closure_Error", "Atmosphere_Closure_Error",
    "Carbon_Density_Adjustment", "PFT_Remap_Adjustment",
    "biomass_total", "soil_total", "P1", "P10", "P100",
    "atmosphere", "system_carbon_total",
    "area_clearing", "area_abandonment", "area_other", "area_harvest",
    "area_deforestation", "transition_requested_area",
    "transition_applied_area", "transition_clipped_area",
    "harvest_luh2_area", "harvest_effective_area", "harvest_extra_area",
]

_ALL_OUTPUT_KEYS = set(STATE_OUTPUT_KEYS + FORMAL_OUTPUT_KEYS)
_NON_ADDITIVE_DEFAULT_EXCLUSIONS = {
    "harvest_luh2_area_frac",
    "harvest_effective_area_frac",
}

_MODEL_CODE_FILES = (
    "src/LULCCSimulator.py",
    "src/parameter_loader.py",
    "src/carbon_pools_init.py",
    "src/events.py",
    "src/harvest.py",
    "src/transition.py",
    "src/summary_yearly.py",
    "src/file_loader.py",
    "src/run_config.py",
    "src/run_manager.py",
    "src/mc_sampling.py",
    "src/mc_runner.py",
)


def _resolve(repo_root: Path, value: Union[str, Path]) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = repo_root / path
    return path.resolve()


def mc_section(config) -> Dict[str, Any]:
    section = config.raw.get("monte_carlo", {}) or {}
    if not isinstance(section, dict):
        raise ValueError("monte_carlo section must be a mapping")
    return section


def mc_sample_table_path(config) -> Path:
    section = mc_section(config)
    return _resolve(
        config.repo_root,
        section.get("sample_table", "config/mc_samples.parquet"),
    )


def mc_output_root(config) -> Path:
    section = mc_section(config)
    configured = section.get("output_dir")
    if configured:
        out_root = _resolve(config.repo_root, configured)
    else:
        out_root = config.output_root / "mc"
    name = str(section.get("name", config.run_name))
    return out_root / config.resolution / name


def sample_directory(config, sample_id: int) -> Path:
    return mc_output_root(config) / "sample_{:06d}".format(int(sample_id))


def band_output_paths(config, sample_id: int, band_id: int) -> Tuple[Path, Path]:
    band_dir = sample_directory(config, sample_id) / "bands"
    stem = "band_{:03d}".format(int(band_id))
    return band_dir / (stem + ".parquet"), band_dir / (stem + ".json")


def global_output_path(config, sample_id: int) -> Path:
    return sample_directory(config, sample_id) / "global.parquet"


def sample_manifest_path(config, sample_id: int) -> Path:
    return sample_directory(config, sample_id) / "sample_manifest.json"


def model_code_sha256(repo_root: Path) -> str:
    digest = hashlib.sha256()
    found = 0
    for rel in _MODEL_CODE_FILES:
        path = repo_root / rel
        if not path.is_file():
            continue
        found += 1
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("ascii"))
        digest.update(b"\0")
    if found == 0:
        raise RuntimeError("Could not hash model code under {}".format(repo_root))
    return digest.hexdigest()


def _load_sample_row(sample_table: Path, sample_id: int) -> Dict[str, Any]:
    if not sample_table.is_file():
        raise FileNotFoundError(
            "Monte Carlo sample table not found: {}".format(sample_table)
        )
    frame = pd.read_parquet(sample_table, engine="pyarrow")
    if "sample_id" not in frame.columns or "overrides_json" not in frame.columns:
        raise ValueError(
            "Sample table must contain sample_id and overrides_json columns: {}".format(
                sample_table
            )
        )
    selected = frame.loc[frame["sample_id"].astype(int) == int(sample_id)]
    if len(selected) != 1:
        raise ValueError(
            "Expected exactly one row for sample_id={}, found {}".format(
                sample_id, len(selected)
            )
        )
    return selected.iloc[0].to_dict()


def _output_variables(config) -> List[str]:
    section = mc_section(config)
    variables = [
        str(value)
        for value in section.get("output_variables", DEFAULT_MC_OUTPUT_VARIABLES)
    ]
    if not variables:
        raise ValueError("monte_carlo.output_variables cannot be empty")
    if len(set(variables)) != len(variables):
        raise ValueError("Duplicate names in monte_carlo.output_variables")

    unknown = [name for name in variables if name not in _ALL_OUTPUT_KEYS]
    if unknown:
        raise ValueError("Unknown MC output variables: {}".format(unknown))

    unsafe = [name for name in variables if name in _NON_ADDITIVE_DEFAULT_EXCLUSIONS]
    if unsafe and not bool(section.get("allow_nonadditive_variables", False)):
        raise ValueError(
            "Non-additive fraction variables cannot be spatially summed: {}".format(
                unsafe
            )
        )
    return variables


def current_mc_identity(
    config,
    sample_table_sha256: str,
    override_sha256: str,
) -> Dict[str, Any]:
    # sample_table_sha256 and git_commit are stored as provenance but are not
    # result-validity keys. This lets an ensemble be extended with new samples,
    # or an identical working tree be committed, without invalidating unchanged
    # sample results. The sample-specific override hash and direct code hash are
    # the result-affecting identities.
    base = build_run_metadata(config)
    return {
        "override_sha256": str(override_sha256),
        "run_config_sha256": hashlib.sha256(
            config.raw_text.encode("utf-8")
        ).hexdigest(),
        "model_code_sha256": model_code_sha256(config.repo_root),
        "parameter_sha256": base.get("parameter_sha256"),
        "experiment_sha256": base.get("experiment_sha256"),
        "state_fingerprint": base.get("state_fingerprint"),
        "transition_fingerprint": base.get("transition_fingerprint"),
        "pft_fingerprint": base.get("pft_fingerprint"),
    }


def _same(a: Any, b: Any) -> bool:
    return canonical_json(a) == canonical_json(b)


def _completion_ok(
    parquet_path: Path,
    metadata_path: Path,
    *,
    sample_id: int,
    band_id: int,
    identity: Dict[str, Any],
    expected_years: Sequence[int],
    variables: Sequence[str],
) -> bool:
    if not parquet_path.is_file() or not metadata_path.is_file():
        return False
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if int(metadata.get("sample_id", -1)) != int(sample_id):
            return False
        if int(metadata.get("band_id", -1)) != int(band_id):
            return False

        for key, expected in identity.items():
            if not _same(metadata.get(key), expected):
                return False

        if list(metadata.get("output_variables", [])) != list(variables):
            return False

        frame = pd.read_parquet(parquet_path, engine="pyarrow")
        if list(frame.columns) != ["year"] + list(variables):
            return False
        return np.array_equal(
            frame["year"].to_numpy(dtype=int),
            np.asarray(expected_years, dtype=int),
        )
    except Exception:
        return False


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


def run_mc_band(
    config_path: Union[str, Path],
    *,
    sample_id: Optional[int] = None,
    band_id: Optional[int] = None,
    expected_resolution: Optional[str] = None,
) -> Path:
    config = load_run_config(
        config_path, expected_resolution=expected_resolution
    )
    if sample_id is None:
        sample_id = int(os.environ.get("MC_SAMPLE_ID", "0"))
    if band_id is None:
        band_id = int(
            os.environ.get(
                "BAND_ID",
                os.environ.get("SLURM_ARRAY_TASK_ID", "0"),
            )
        )

    sample_table = mc_sample_table_path(config)
    sample_table_sha = sha256_file(sample_table)
    row = _load_sample_row(sample_table, sample_id)

    sample_overrides = json.loads(str(row["overrides_json"]))
    if not isinstance(sample_overrides, dict):
        raise ValueError("Sample overrides_json must decode to a mapping")

    override_sha = str(
        row.get(
            "override_sha256",
            hashlib.sha256(
                canonical_json(sample_overrides).encode("utf-8")
            ).hexdigest(),
        )
    )
    sample_seed = int(row.get("sample_seed", -1))
    sample_kind = str(row.get("sample_kind", "mc"))

    parameter_overrides = deep_merge(
        config.raw.get("model_overrides", {}) or {},
        sample_overrides,
    )

    lat_slice, lon_slice, bands_total = _band_slice(config, int(band_id))
    variables = _output_variables(config)
    years = np.arange(
        config.start_year, config.end_year + 1, dtype=np.int32
    )

    identity = current_mc_identity(config, sample_table_sha, override_sha)
    base_metadata = build_run_metadata(config)
    parquet_path, metadata_path = band_output_paths(
        config, sample_id, int(band_id)
    )

    if _completion_ok(
        parquet_path,
        metadata_path,
        sample_id=int(sample_id),
        band_id=int(band_id),
        identity=identity,
        expected_years=years,
        variables=variables,
    ):
        print("[SKIP] complete and current MC band exists: {}".format(parquet_path))
        return parquet_path

    if parquet_path.exists():
        parquet_path.unlink()
    if metadata_path.exists():
        metadata_path.unlink()

    run_metadata = dict(base_metadata)
    run_metadata.update(
        {
            "mc_mode": "summary_only",
            "mc_sample_id": int(sample_id),
            "mc_sample_kind": sample_kind,
            "mc_sample_seed": sample_seed,
            "mc_override_sha256": override_sha,
            "mc_sample_table": str(sample_table),
            "mc_sample_table_sha256": sample_table_sha,
            "mc_model_code_sha256": identity["model_code_sha256"],
            "band_id": int(band_id),
            "bands_total": int(bands_total),
            "band_size_deg": float(config.band_size_deg),
            "server_mode": 1,
        }
    )

    print(
        "[MC RUN] sample={} ({}) band={}/{} years={}-{} variables={}".format(
            sample_id,
            sample_kind,
            band_id,
            bands_total - 1,
            config.start_year,
            config.end_year,
            len(variables),
        ),
        flush=True,
    )

    simulator = LULCCSimulator(
        config_path=str(config.base_parameter_path),
        experiment_path=str(config.experiment_path),
        LULC_path=str(config.state_path),
        trans_path=str(config.transition_path),
        pft_path=str(config.pft_path),
        input_format=config.input_format,
        input_base_year=config.input_base_year,
        start_year=config.start_year,
        end_year=config.end_year,
        time_encoding=config.time_encoding,
        pft_var=config.pft_variable,
        pft_base_year=config.pft_base_year,
        pft_time_encoding=config.pft_time_encoding,
        pft_update_mode=config.pft_update_mode,
        pft_min_year_policy=config.pft_min_year_policy,
        pft_max_year_policy=config.pft_max_year_policy,
        pft_missing_cell_policy=config.pft_missing_cell_policy,
        pft_nearest_search_radius=config.pft_nearest_search_radius,
        pft_default_index=config.pft_default_index,
        lat_slice=lat_slice,
        lon_slice=lon_slice,
        area_unit=config.area_unit,
        run_metadata=run_metadata,
        parameter_overrides=parameter_overrides,
        compression_level=config.compression_level,
    )

    nlat, nlon = simulator._area_grid.shape
    totals = np.zeros((len(years), len(variables)), dtype=np.float64)
    total_cells = nlat * nlon
    completed = 0
    land_cells = 0
    started = time.time()

    for i in range(nlat):
        for j in range(nlon):
            yearly = simulator.run_simulation(i, j)
            nonzero_cell = False
            for t, record in enumerate(yearly):
                values = [float(record.get(name, 0.0)) for name in variables]
                totals[t, :] += values
                if not nonzero_cell and any(value != 0.0 for value in values):
                    nonzero_cell = True
            if nonzero_cell:
                land_cells += 1

            completed += 1
            if completed % 250 == 0 or completed == total_cells:
                elapsed = max(time.time() - started, 1.0e-9)
                rate = completed / elapsed
                eta = (
                    (total_cells - completed) / rate
                    if rate > 0.0
                    else 0.0
                )
                print(
                    "\r[MC band] {}/{} ({:.1f}%) ETA={:.2f}h".format(
                        completed,
                        total_cells,
                        100.0 * completed / total_cells,
                        eta / 3600.0,
                    ),
                    end="",
                    flush=True,
                )
    print("", flush=True)

    frame = pd.DataFrame(totals, columns=variables)
    frame.insert(0, "year", years)
    _atomic_parquet(frame, parquet_path)

    metadata = {
        "schema_version": 2,
        "mode": "summary_only",
        "sample_id": int(sample_id),
        "sample_kind": sample_kind,
        "sample_seed": sample_seed,
        "band_id": int(band_id),
        "bands_total": int(bands_total),
        "band_size_deg": float(config.band_size_deg),
        "output_variables": variables,
        "start_year": int(config.start_year),
        "end_year": int(config.end_year),
        "n_cells": int(total_cells),
        "n_nonzero_output_cells": int(land_cells),
        "elapsed_seconds": float(time.time() - started),
        "parameter_overrides": parameter_overrides,
        "sample_overrides": sample_overrides,
        "sample_table_sha256": sample_table_sha,
        "git_commit": base_metadata.get("git_commit", "unknown"),
    }
    metadata.update(identity)
    _atomic_json(metadata, metadata_path)

    print("[DONE] {}".format(parquet_path), flush=True)
    return parquet_path
