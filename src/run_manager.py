"""Common local/server run orchestration."""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Union
import random
import time

import netCDF4 as nc
import numpy as np

from src.LULCCSimulator import LULCCSimulator
from src.file_loader import FileLoader
from src.input_validator import validate_inputs
from src.run_config import (
    build_run_metadata,
    load_run_config,
    manifest_identity,
    write_run_manifest,
)


def _canonical_metadata_value(value: Any) -> str:
    """Normalize Python and NetCDF attribute values for stable comparison."""

    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.startswith("{") or stripped.startswith("["):
            try:
                value = json.loads(stripped)
            except json.JSONDecodeError:
                pass
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )


def output_complete(
    path: Path,
    expected_time: int,
    expected_lat: int,
    expected_lon: int,
    expected_metadata: Optional[Dict[str, Any]] = None,
) -> bool:
    if not path.exists():
        return False
    try:
        with nc.Dataset(path, "r") as ds:
            actual = tuple(len(ds.dimensions[name]) for name in ("time", "lat", "lon"))
            if actual != (expected_time, expected_lat, expected_lon):
                return False

            if expected_metadata is not None:
                for key, expected in manifest_identity(expected_metadata).items():
                    if key not in ds.ncattrs():
                        return False
                    actual_value = ds.getncattr(key)
                    if _canonical_metadata_value(actual_value) != _canonical_metadata_value(expected):
                        return False

            if "done" not in ds.variables:
                return False
            done = np.asarray(ds.variables["done"][:])
            return done.shape == (expected_lat, expected_lon) and bool(np.all(done == 1))
    except Exception:
        return False


_RUNTIME_METADATA_KEYS = (
    "pft_variable",
    "pft_source_mode",
    "pft_dynamic",
    "transition_channel_count",
    "state_channel_count",
)


def _read_runtime_metadata(path: Path) -> Dict[str, Any]:
    """Read runtime-resolved metadata from an existing NetCDF output."""

    runtime: Dict[str, Any] = {}
    try:
        with nc.Dataset(path, "r") as ds:
            for key in _RUNTIME_METADATA_KEYS:
                if key in ds.ncattrs():
                    value = ds.getncattr(key)
                    if isinstance(value, np.generic):
                        value = value.item()
                    if isinstance(value, bytes):
                        value = value.decode("utf-8", errors="replace")
                    runtime[key] = value
    except Exception:
        return {}
    return runtime


def _persist_runtime_manifests(
    *,
    config,
    output_dir: Path,
    base_metadata: Dict[str, Any],
    runtime_metadata: Dict[str, Any],
    server_mode: bool,
    band_id: Optional[int],
) -> None:
    """Persist runtime-resolved metadata without shared-file write races."""

    if server_mode:
        if band_id is None:
            raise ValueError("band_id is required in server mode")
        write_run_manifest(
            config,
            output_dir,
            runtime_metadata,
            manifest_name="run_manifest.rank{:03d}.json".format(band_id),
            overwrite_manifest=True,
        )
        if band_id == 0:
            shared_runtime_metadata = dict(base_metadata)
            for key in _RUNTIME_METADATA_KEYS:
                shared_runtime_metadata[key] = runtime_metadata.get(key)
            write_run_manifest(
                config,
                output_dir,
                shared_runtime_metadata,
                manifest_name="run_manifest.json",
                overwrite_manifest=True,
            )
    else:
        write_run_manifest(
            config,
            output_dir,
            runtime_metadata,
            manifest_name="run_manifest.json",
            overwrite_manifest=True,
        )


def _band_slice(config, band_id: int):
    info = FileLoader().inspect_dataset(config.state_path)
    lat_res = abs(float(np.median(np.diff(info.lat))))
    lon_res = abs(float(np.median(np.diff(info.lon))))
    nlat, nlon = len(info.lat), len(info.lon)
    expected_nlon = int(round(360.0 / lon_res))
    if nlon != expected_nlon:
        raise ValueError(
            f"Longitude grid is not global/regular: nlon={nlon}, expected={expected_nlon}"
        )
    if config.band_size_deg <= 0 or abs(360.0 / config.band_size_deg - round(360.0 / config.band_size_deg)) > 1e-8:
        raise ValueError(f"server.band_size_deg must divide 360: {config.band_size_deg}")
    bands_total = int(round(360.0 / config.band_size_deg))
    if not (0 <= band_id < bands_total):
        raise ValueError(f"band_id={band_id} outside 0..{bands_total-1}")

    lat_internal = np.linspace(-90 + 0.5 * lat_res, 90 - 0.5 * lat_res, nlat)
    lon_internal = np.linspace(-180 + 0.5 * lon_res, 180 - 0.5 * lon_res, nlon)
    lon_min = -180.0 + band_id * config.band_size_deg
    lon_max = lon_min + config.band_size_deg - 1e-9
    i0 = int(np.searchsorted(lat_internal, -90.0, side="left"))
    i1 = int(np.searchsorted(lat_internal, 90.0, side="right"))
    j0 = int(np.searchsorted(lon_internal, lon_min, side="left"))
    j1 = int(np.searchsorted(lon_internal, lon_max, side="right"))
    return slice(i0, i1), slice(j0, j1), bands_total


def run_from_config(
    config_path: Union[str, Path],
    *,
    expected_resolution: Optional[str] = None,
    server_mode: bool = False,
) -> Path:
    config = load_run_config(config_path, expected_resolution=expected_resolution)
    base_metadata = build_run_metadata(config)
    metadata = dict(base_metadata)
    output_dir = config.output_root / config.resolution / config.run_name
    output_dir.mkdir(parents=True, exist_ok=True)

    if config.validate_before_run and server_mode:
        raise RuntimeError(
            "Automatic full input validation is disabled for server array jobs "
            "because every band would otherwise rescan the same global inputs. "
            "Run `python3 -m tools.validate_inputs` once before sbatch, confirm "
            "PASS, and set validation.run_before_simulation=false in the server "
            "run YAML."
        )

    if server_mode:
        band_id = int(os.environ.get("BAND_ID", os.environ.get("SLURM_ARRAY_TASK_ID", "0")))
        if config.stagger_max > 0:
            random.seed(os.getpid() + band_id * 10007)
            delay = random.uniform(0.0, config.stagger_max)
            print(f"[STAGGER] sleeping {delay:.1f}s", flush=True)
            time.sleep(delay)
        lat_slice, lon_slice, bands_total = _band_slice(config, band_id)
        output_path = output_dir / f"{config.output_prefix}.rank{band_id:03d}.nc"
        metadata.update({
            "server_mode": 1,
            "band_id": band_id,
            "bands_total": bands_total,
            "band_size_deg": config.band_size_deg,
        })
    else:
        band_id = None
        lat_slice = None
        lon_slice = None
        output_path = output_dir / f"{config.output_prefix}.global.nc"
        metadata["server_mode"] = 0

    metadata["output_file"] = str(output_path)

    if config.validate_before_run:
        report = validate_inputs(
            config,
            report_path=output_dir / "input_validation_report.json",
        )
        if report["status"] != "PASS":
            raise RuntimeError(
                "Input validation failed. See input_validation_report.json in the run directory."
            )

    if server_mode:
        # Shared manifest: common information for the complete global run.
        write_run_manifest(
            config,
            output_dir,
            base_metadata,
            manifest_name="run_manifest.json",
            overwrite_manifest=False,
        )

        # Band-specific manifest: no shared-file overwrite between array tasks.
        write_run_manifest(
            config,
            output_dir,
            metadata,
            manifest_name="run_manifest.rank{:03d}.json".format(band_id),
            overwrite_manifest=True,
        )
    else:
        # Verify an existing local run before deciding whether output can be reused.
        write_run_manifest(
            config,
            output_dir,
            metadata,
            manifest_name="run_manifest.json",
            overwrite_manifest=False,
        )

    state_info = FileLoader().inspect_dataset(config.state_path)
    expected_lat = len(state_info.lat) if lat_slice is None else lat_slice.stop - lat_slice.start
    expected_lon = len(state_info.lon) if lon_slice is None else lon_slice.stop - lon_slice.start
    expected_time = config.years + 1

    if output_complete(
        output_path,
        expected_time,
        expected_lat,
        expected_lon,
        expected_metadata=base_metadata,
    ):
        # Backfill runtime metadata from the completed NetCDF when an existing
        # run is skipped, so manifests remain complete without reloading inputs.
        existing_runtime_metadata = dict(metadata)
        existing_runtime_metadata.update(_read_runtime_metadata(output_path))
        _persist_runtime_manifests(
            config=config,
            output_dir=output_dir,
            base_metadata=base_metadata,
            runtime_metadata=existing_runtime_metadata,
            server_mode=server_mode,
            band_id=band_id,
        )
        print(f"[SKIP] complete output exists: {output_path}")
        return output_path
    if output_path.exists():
        print(f"[RESTART] removing incomplete output: {output_path}")
        output_path.unlink()

    print(
        f"[RUN] name={config.run_name}, resolution={config.resolution}, "
        f"format={config.input_format}, years={config.start_year}-{config.end_year}, "
        f"band={band_id if band_id is not None else 'global'}",
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
        lat_slice=lat_slice,
        lon_slice=lon_slice,
        area_unit=config.area_unit,
        run_metadata=metadata,
        parameter_overrides=config.raw.get("model_overrides", {}) or {},
        compression_level=config.compression_level,
    )

    # The simulator adds runtime-resolved information such as the actual PFT
    # variable, PFT layout, and compiled channel counts. Persist those values
    # after initialization instead of leaving them only in NetCDF attributes.
    runtime_metadata = dict(simulator.run_metadata)
    _persist_runtime_manifests(
        config=config,
        output_dir=output_dir,
        base_metadata=base_metadata,
        runtime_metadata=runtime_metadata,
        server_mode=server_mode,
        band_id=band_id,
    )

    simulator.run_simulation_grid(
        out_nc=str(output_path),
        lat_slice=None,
        lon_slice=None,
        sync_every=config.sync_every,
    )
    if not output_complete(
        output_path,
        expected_time,
        expected_lat,
        expected_lon,
        expected_metadata=base_metadata,
    ):
        raise RuntimeError(f"Output failed completion check: {output_path}")
    return output_path
