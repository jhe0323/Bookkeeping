"""One-time validation for model inputs and configuration."""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Union

import netCDF4 as nc
import numpy as np

from src.file_loader import FileLoader, calendarize_time_values
from src.parameter_loader import ParameterLoader
from src.run_config import ResolvedRunConfig, load_run_config


_TRANSITION_RE = re.compile(r"^(.+)_to_(.+?)(?:_frac)?$")
_HARVEST_FAMILIES = ("primf", "primn", "secmf", "secyf", "secnf")


def _coord_equal(a, b, *, longitude=False, atol=1e-7) -> bool:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if longitude:
        a = ((a + 180) % 360) - 180
        b = ((b + 180) % 360) - 180
    return a.size == b.size and np.allclose(np.sort(a), np.sort(b), atol=atol, rtol=0)


def _variable_stats(variable, *, full_scan: bool, chunk_time: int = 16, chunk_lat: int = 128):
    shape = variable.shape
    minimum = np.inf
    maximum = -np.inf
    nonfinite = 0
    count = 0

    def consume(data):
        nonlocal minimum, maximum, nonfinite, count
        if hasattr(data, "filled"):
            data = data.filled(np.nan)
        data = np.asarray(data, dtype=np.float64)
        finite = np.isfinite(data)
        nonfinite += int(data.size - finite.sum())
        count += int(data.size)
        if finite.any():
            minimum = min(minimum, float(data[finite].min()))
            maximum = max(maximum, float(data[finite].max()))

    if not full_scan:
        index = tuple(slice(0, min(size, 2)) for size in shape)
        consume(variable[index])
    else:
        dims = variable.dimensions
        time_axis = next((i for i, d in enumerate(dims) if d.lower() in {"time", "year", "time_counter"}), None)
        lat_axis = next((i for i, d in enumerate(dims) if d.lower() in {"lat", "latitude", "y"}), None)
        t_len = shape[time_axis] if time_axis is not None else 1
        y_len = shape[lat_axis] if lat_axis is not None else 1
        for t0 in range(0, t_len, chunk_time):
            for y0 in range(0, y_len, chunk_lat):
                index = [slice(None)] * len(shape)
                if time_axis is not None:
                    index[time_axis] = slice(t0, min(t_len, t0 + chunk_time))
                if lat_axis is not None:
                    index[lat_axis] = slice(y0, min(y_len, y0 + chunk_lat))
                consume(variable[tuple(index)])

    return {
        "min": None if not np.isfinite(minimum) else minimum,
        "max": None if not np.isfinite(maximum) else maximum,
        "nonfinite": nonfinite,
        "count": count,
    }


def _pft_structure_stats(variable, *, pft_axis: Optional[int], full_scan: bool):
    shape = variable.shape
    dims = variable.dimensions
    time_axis = next((i for i, d in enumerate(dims) if d.lower() in {"time", "year", "time_counter"}), None)
    lat_axis = next((i for i, d in enumerate(dims) if d.lower() in {"lat", "latitude", "y"}), None)
    t_len = shape[time_axis] if time_axis is not None else 1
    y_len = shape[lat_axis] if lat_axis is not None else 1
    t_step = 16 if full_scan else min(2, t_len)
    y_step = 128 if full_scan else min(2, y_len)
    max_sum_error = 0.0
    nonzero_cells = 0
    max_integer_error = 0.0

    for t0 in range(0, t_len, t_step):
        for y0 in range(0, y_len, y_step):
            index = [slice(None)] * len(shape)
            if time_axis is not None:
                index[time_axis] = slice(t0, min(t_len, t0 + t_step))
            if lat_axis is not None:
                index[lat_axis] = slice(y0, min(y_len, y0 + y_step))
            data = variable[tuple(index)]
            if hasattr(data, "filled"):
                data = data.filled(np.nan)
            data = np.asarray(data, dtype=np.float64)
            if pft_axis is not None:
                sums = np.nansum(data, axis=pft_axis)
                active = sums > 0.0
                if active.any():
                    max_sum_error = max(max_sum_error, float(np.max(np.abs(sums[active] - 1.0))))
                    nonzero_cells += int(active.sum())
            else:
                finite = np.isfinite(data)
                if finite.any():
                    max_integer_error = max(
                        max_integer_error,
                        float(np.max(np.abs(data[finite] - np.rint(data[finite])))),
                    )
            if not full_scan:
                return {
                    "max_sum_error": max_sum_error,
                    "nonzero_cells": nonzero_cells,
                    "max_integer_error": max_integer_error,
                }
    return {
        "max_sum_error": max_sum_error,
        "nonzero_cells": nonzero_cells,
        "max_integer_error": max_integer_error,
    }


def _required_state_variables(cfg: ResolvedRunConfig, params: ParameterLoader) -> List[str]:
    if cfg.input_format == "vscp":
        return ["v", "s", "p", "c"]
    return list((params.config.get("LUH2toLULC", {}) or {}).keys())


def validate_inputs(
    config: ResolvedRunConfig,
    *,
    report_path: Optional[Union[str, Path]] = None,
    full_scan: Optional[bool] = None,
) -> dict:
    full_scan = config.validation_full_scan if full_scan is None else bool(full_scan)
    loader = FileLoader()
    params = ParameterLoader(
        str(config.base_parameter_path),
        experiment_path=str(config.experiment_path),
        override_config=config.raw.get("model_overrides", {}) or {},
    )
    state_info = loader.inspect_dataset(config.state_path)
    trans_info = loader.inspect_dataset(config.transition_path)
    pft_info = loader.inspect_dataset(config.pft_path)

    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_config": str(config.config_path),
        "full_scan": full_scan,
        "status": "PASS",
        "errors": [],
        "warnings": [],
        "summary": {},
        "variables": {},
    }

    def error(message: str):
        report["errors"].append(message)
        report["status"] = "FAIL"

    def warning(message: str):
        report["warnings"].append(message)

    if not _coord_equal(state_info.lat, trans_info.lat):
        error("States and transitions latitude coordinates differ.")
    if not _coord_equal(state_info.lon, trans_info.lon, longitude=True):
        error("States and transitions longitude coordinates differ.")

    state_time_len = state_info.dimensions.get(state_info.time_name, 1) if state_info.time_name else 1
    trans_time_len = trans_info.dimensions.get(trans_info.time_name, 1) if trans_info.time_name else 1
    try:
        state_years = calendarize_time_values(
            state_info.time_values,
            length=state_time_len,
            input_base_year=config.input_base_year,
            encoding=config.time_encoding,
        )
        trans_years = calendarize_time_values(
            trans_info.time_values,
            length=trans_time_len,
            input_base_year=config.input_base_year,
            encoding=config.time_encoding,
        )
        if config.start_year not in set(state_years.tolist()):
            error(f"Start year {config.start_year} is absent from states time axis.")
        required_transition_years = set(range(config.start_year, config.end_year))
        missing_years = sorted(required_transition_years - set(trans_years.tolist()))
        if missing_years:
            error(
                f"Transitions time axis lacks {len(missing_years)} required years; "
                f"first missing={missing_years[:5]}"
            )
        report["summary"]["state_year_range"] = [int(state_years.min()), int(state_years.max())]
        report["summary"]["transition_year_range"] = [int(trans_years.min()), int(trans_years.max())]
    except Exception as exc:
        error(f"Time-axis validation failed: {exc}")

    state_required = _required_state_variables(config, params)
    state_missing = sorted(set(state_required) - set(state_info.variables))
    if state_missing:
        error(f"States file is missing variables: {state_missing}")

    transition_channels = []
    mapping = params.config.get("LUH2toLULC", {}) or {}
    for name in trans_info.variables:
        match = _TRANSITION_RE.match(name)
        if not match:
            continue
        src, dst = match.groups()
        if config.input_format == "vscp":
            if src in {"v", "s", "p", "c", "U"} and dst in {"v", "s", "p", "c", "U"} and src != dst:
                transition_channels.append(name)
        elif src in mapping and dst in mapping and mapping[src] != mapping[dst]:
            transition_channels.append(name)
    if not transition_channels:
        error("No usable transition channels were found.")
    report["summary"]["transition_channel_count"] = len(transition_channels)

    area_vars = [f"{family}_harv" for family in _HARVEST_FAMILIES]
    bio_vars = [f"{family}_bioh" for family in _HARVEST_FAMILIES]
    if not any(name in trans_info.variables for name in area_vars):
        error("No LUH2 harvest-area variable was found.")
    if params.harvest_use_bioh and not any(name in trans_info.variables for name in bio_vars):
        error("The selected biomass-demand experiment has no *_bioh variables.")

    with nc.Dataset(config.state_path, "r") as ds:
        for name in state_required:
            if name not in ds.variables:
                continue
            stats = _variable_stats(ds.variables[name], full_scan=full_scan)
            report["variables"][f"states:{name}"] = stats
            if stats["min"] is not None and stats["min"] < -1e-8:
                error(f"State variable {name} contains negative values: min={stats['min']}")
            if stats["max"] is not None and stats["max"] > 1 + 1e-6:
                error(f"State variable {name} exceeds 1: max={stats['max']}")
            if stats["nonfinite"]:
                warning(f"State variable {name} has {stats['nonfinite']} masked/nonfinite values.")

        # State-sum check, performed in chunks to avoid loading the whole file.
        present = [name for name in state_required if name in ds.variables]
        if present:
            first = ds.variables[present[0]]
            dims = first.dimensions
            time_axis = next((i for i, d in enumerate(dims) if d.lower() in {"time", "year", "time_counter"}), None)
            lat_axis = next((i for i, d in enumerate(dims) if d.lower() in {"lat", "latitude", "y"}), None)
            t_len = first.shape[time_axis] if time_axis is not None else 1
            y_len = first.shape[lat_axis] if lat_axis is not None else 1
            min_sum = np.inf
            max_sum = -np.inf
            for t0 in range(0, t_len, 16):
                for y0 in range(0, y_len, 128):
                    index = [slice(None)] * len(first.shape)
                    if time_axis is not None:
                        index[time_axis] = slice(t0, min(t_len, t0 + 16))
                    if lat_axis is not None:
                        index[lat_axis] = slice(y0, min(y_len, y0 + 128))
                    total = None
                    for name in present:
                        data = ds.variables[name][tuple(index)]
                        if hasattr(data, "filled"):
                            data = data.filled(0.0)
                        data = np.asarray(data, dtype=np.float64)
                        total = data if total is None else total + data
                    finite = np.isfinite(total)
                    if finite.any():
                        min_sum = min(min_sum, float(total[finite].min()))
                        max_sum = max(max_sum, float(total[finite].max()))
            report["summary"]["state_sum_min"] = None if not np.isfinite(min_sum) else min_sum
            report["summary"]["state_sum_max"] = None if not np.isfinite(max_sum) else max_sum
            if np.isfinite(max_sum) and max_sum > 1 + 1e-5:
                error(f"Land-use state sum exceeds 1: max={max_sum}")

    with nc.Dataset(config.transition_path, "r") as ds:
        scan_names = transition_channels + [name for name in area_vars if name in ds.variables]
        if params.harvest_use_bioh:
            scan_names += [name for name in bio_vars if name in ds.variables]
        for name in scan_names:
            stats = _variable_stats(ds.variables[name], full_scan=full_scan)
            report["variables"][f"transitions:{name}"] = stats
            if stats["min"] is not None and stats["min"] < -1e-8:
                error(f"Transition variable {name} contains negative values: min={stats['min']}")
            if not name.endswith("_bioh") and stats["max"] is not None and stats["max"] > 1 + 1e-6:
                error(f"Fraction variable {name} exceeds 1: max={stats['max']}")

    # PFT structure and range.
    with nc.Dataset(config.pft_path, "r") as ds:
        pft_name = config.pft_variable
        if pft_name is None:
            preferred = (
                "pft_fraction", "pft_frac", "pft_fractions", "vegetation_fraction",
                "vegetfrac", "maxvegetfrac", "pft", "PFT", "dominant_pft",
                "pft_map", "pft_dominant",
            )
            pft_name = next((name for name in preferred if name in ds.variables), None)
        if pft_name is None or pft_name not in ds.variables:
            error("Cannot identify PFT data variable; set inputs.pft_variable.")
        else:
            variable = ds.variables[pft_name]
            stats = _variable_stats(variable, full_scan=full_scan)
            report["variables"][f"pft:{pft_name}"] = stats
            non_spatial = [
                (dim, size)
                for dim, size in zip(variable.dimensions, variable.shape)
                if dim.lower() not in {"lat", "latitude", "y", "lon", "longitude", "x", "time", "year", "time_counter"}
            ]
            pft_axis = next(
                (
                    axis for axis, (dim, size) in enumerate(zip(variable.dimensions, variable.shape))
                    if dim.lower() in {"pft", "pfts", "npft", "veget", "veg", "vegtype", "vegetation", "class", "classes", "type", "types"}
                    or (
                        dim.lower() not in {"lat", "latitude", "y", "lon", "longitude", "x", "time", "year", "time_counter"}
                        and int(size) == params.n_pft
                    )
                ),
                None,
            )
            if non_spatial and pft_axis is None:
                error(
                    f"PFT variable non-spatial dimensions {non_spatial} do not contain "
                    f"configured n_pft={params.n_pft}."
                )
            if stats["min"] is not None and stats["min"] < -1e-8:
                error(f"PFT variable contains negative values: min={stats['min']}")

            structure = _pft_structure_stats(
                variable,
                pft_axis=pft_axis,
                full_scan=full_scan,
            )
            report["summary"]["pft_structure"] = structure
            if pft_axis is not None:
                if stats["max"] is not None and stats["max"] > 1 + 1e-6:
                    error(f"Fractional PFT variable exceeds 1: max={stats['max']}")
                if structure["max_sum_error"] > 1e-4:
                    error(
                        "Fractional PFT sums differ from 1 in active cells: "
                        f"max_error={structure['max_sum_error']}"
                    )
            else:
                if stats["max"] is not None and stats["max"] > params.n_pft + 1e-6:
                    error(
                        f"Dominant PFT code exceeds configured class count {params.n_pft}: "
                        f"max={stats['max']}"
                    )
                if structure["max_integer_error"] > 1e-5:
                    error(
                        "Dominant PFT map contains non-integer codes: "
                        f"max_error={structure['max_integer_error']}"
                    )

            time_axis = next(
                (axis for axis, dim in enumerate(variable.dimensions) if dim.lower() in {"time", "year", "time_counter"}),
                None,
            )
            if time_axis is not None and variable.shape[time_axis] > 1:
                time_name = variable.dimensions[time_axis]
                raw_time = np.asarray(ds.variables[time_name][:]) if time_name in ds.variables else None
                pft_years = calendarize_time_values(
                    raw_time,
                    length=variable.shape[time_axis],
                    input_base_year=config.pft_base_year,
                    encoding=config.pft_time_encoding,
                )
                report["summary"]["pft_year_range"] = [int(pft_years.min()), int(pft_years.max())]
                if config.start_year < pft_years.min() or config.end_year > pft_years.max():
                    warning(
                        "Dynamic PFT years do not fully cover the simulation and will use "
                        f"configured clip policies: source={pft_years.min()}..{pft_years.max()}"
                    )

    report["summary"].update({
        "input_format": config.input_format,
        "configured_pft_count": params.n_pft,
        "forest_pft_count": len(params.get_forest_pfts()),
        "state_grid": [len(state_info.lat), len(state_info.lon)],
        "transition_grid": [len(trans_info.lat), len(trans_info.lon)],
        "pft_grid": [len(pft_info.lat), len(pft_info.lon)],
    })

    if report_path is None:
        report_path = config.repo_root / "validation_report.json"
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


def validate_run_config_file(path: Union[str, Path]) -> dict:
    config = load_run_config(path)
    output_dir = config.output_root / config.resolution / config.run_name
    return validate_inputs(
        config,
        report_path=output_dir / "input_validation_report.json",
    )
