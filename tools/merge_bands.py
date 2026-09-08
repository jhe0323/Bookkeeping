#!/usr/bin/env python3
"""Safely merge deterministic longitude-band NetCDF outputs.

Before writing the global file this tool verifies:
- all expected/selected band files exist in the supplied glob;
- every band has ``done == 1`` for every grid cell;
- time/latitude coordinates and variable schemas agree;
- longitude bands are non-overlapping and form a regular global grid;
- critical run/model metadata agree across bands.

The merge is streamed in time blocks and never loads the full global 0.25-degree
output into memory.

Python 3.8 compatible.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import glob
import json
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import netCDF4 as nc
import numpy as np


IDENTITY_ATTRS = (
    "run_name",
    "resolution",
    "input_format",
    "input_base_year",
    "start_year",
    "end_year",
    "run_config_sha256",
    "parameter_sha256",
    "experiment_sha256",
    "model_code_sha256",
    "state_fingerprint",
    "transition_fingerprint",
    "pft_fingerprint",
    "dynamic_density_file",
    "dynamic_density_fingerprint",
    "pft_update_mode",
    "pft_missing_cell_policy",
)


def _canon(value: Any) -> str:
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


def _wrap_lon(values: np.ndarray) -> np.ndarray:
    return ((np.asarray(values, dtype=float) + 180.0) % 360.0) - 180.0


def _copy_attrs(src, dst, *, exclude: Sequence[str] = ()) -> None:
    excluded = set(exclude)
    for name in src.ncattrs():
        if name in excluded:
            continue
        dst.setncattr(name, src.getncattr(name))


def _critical_metadata(ds: nc.Dataset) -> Dict[str, Any]:
    return {
        key: ds.getncattr(key)
        for key in IDENTITY_ATTRS
        if key in ds.ncattrs()
    }


def _check_done(ds: nc.Dataset, path: Path) -> None:
    if "done" not in ds.variables:
        raise RuntimeError("Missing done variable: {}".format(path))
    done = np.asarray(ds.variables["done"][:])
    if done.ndim != 2 or not np.all(done == 1):
        completed = int(np.count_nonzero(done == 1))
        total = int(done.size)
        raise RuntimeError(
            "Incomplete band {}: done={}/{}".format(path.name, completed, total)
        )


def _variable_schema(ds: nc.Dataset) -> Dict[str, Tuple[str, Tuple[str, ...]]]:
    return {
        name: (str(var.dtype), tuple(var.dimensions))
        for name, var in ds.variables.items()
    }


def _inspect_bands(paths: Sequence[Path]):
    if not paths:
        raise FileNotFoundError("No band files matched the pattern")

    records = []
    ref_time = None
    ref_lat = None
    ref_schema = None
    ref_meta = None

    for path in paths:
        with nc.Dataset(path, "r") as ds:
            for coord in ("time", "lat", "lon"):
                if coord not in ds.variables:
                    raise RuntimeError("{} has no {} coordinate".format(path, coord))
            _check_done(ds, path)

            time = np.asarray(ds.variables["time"][:])
            lat = np.asarray(ds.variables["lat"][:], dtype=float)
            lon_raw = np.asarray(ds.variables["lon"][:], dtype=float)
            lon = _wrap_lon(lon_raw)
            schema = _variable_schema(ds)
            meta = _critical_metadata(ds)

            if ref_time is None:
                ref_time = time
                ref_lat = lat
                ref_schema = schema
                ref_meta = meta
            else:
                if not np.array_equal(time, ref_time):
                    raise RuntimeError("Time coordinate mismatch: {}".format(path))
                if not np.allclose(lat, ref_lat, atol=1e-7, rtol=0):
                    raise RuntimeError("Latitude coordinate mismatch: {}".format(path))
                if schema != ref_schema:
                    raise RuntimeError("Variable schema mismatch: {}".format(path))
                # Compare keys present in either reference or current metadata.
                for key in sorted(set(ref_meta) | set(meta)):
                    if _canon(ref_meta.get(key)) != _canon(meta.get(key)):
                        raise RuntimeError(
                            "Band metadata mismatch for {!r}: {}".format(key, path)
                        )

            records.append((path, lon))

    records.sort(key=lambda item: float(np.min(item[1])))
    all_lon = np.concatenate([lon for _, lon in records])

    if all_lon.size < 2:
        raise RuntimeError("Merged longitude axis is too short")
    diffs = np.diff(all_lon)
    positive = diffs[diffs > 1.0e-9]
    if positive.size == 0:
        raise RuntimeError("Cannot infer longitude resolution")
    lon_res = float(np.median(positive))
    if not np.allclose(diffs, lon_res, atol=max(1.0e-7, lon_res * 1.0e-6), rtol=0):
        raise RuntimeError(
            "Longitude bands have gaps, overlaps, duplicates, or irregular spacing"
        )
    expected_nlon = int(round(360.0 / lon_res))
    if all_lon.size != expected_nlon:
        raise RuntimeError(
            "Longitude coverage is not global: nlon={}, expected={} at resolution {}"
            .format(all_lon.size, expected_nlon, lon_res)
        )

    expected_first = -180.0 + 0.5 * lon_res
    expected_last = 180.0 - 0.5 * lon_res
    if abs(float(all_lon[0]) - expected_first) > 1.0e-5 or abs(float(all_lon[-1]) - expected_last) > 1.0e-5:
        raise RuntimeError(
            "Longitude coverage does not span the expected global cell centers: "
            "first={}, last={}, expected=({}, {})".format(
                all_lon[0], all_lon[-1], expected_first, expected_last
            )
        )

    return records, ref_time, ref_lat, all_lon, ref_schema


def _create_variable_like(
    out: nc.Dataset,
    name: str,
    src_var,
    *,
    compression: int,
    time_block: int,
) -> Any:
    fill_value = getattr(src_var, "_FillValue", None)
    kwargs: Dict[str, Any] = {}

    if name not in {"time", "lat", "lon"}:
        kwargs["zlib"] = True
        kwargs["complevel"] = int(compression)
        if tuple(src_var.dimensions) == ("time", "lat", "lon"):
            kwargs["chunksizes"] = (
                min(len(out.dimensions["time"]), max(1, int(time_block))),
                min(len(out.dimensions["lat"]), 16),
                min(len(out.dimensions["lon"]), 16),
            )
        elif tuple(src_var.dimensions) == ("lat", "lon"):
            kwargs["chunksizes"] = (
                min(len(out.dimensions["lat"]), 64),
                min(len(out.dimensions["lon"]), 64),
            )

    if fill_value is not None:
        kwargs["fill_value"] = fill_value

    dst = out.createVariable(name, src_var.dtype, src_var.dimensions, **kwargs)
    _copy_attrs(src_var, dst, exclude=("_FillValue",))
    return dst


def merge_bands(
    *,
    pattern: str,
    output: Path,
    time_min: int,
    time_block: int,
    compression: int,
    overwrite: bool,
) -> Path:
    paths = [Path(value).resolve() for value in sorted(glob.glob(pattern))]
    records, full_time, lat, lon, _schema = _inspect_bands(paths)

    time_values = np.asarray(full_time)
    selected_idx = np.flatnonzero(time_values >= int(time_min))
    if selected_idx.size == 0:
        raise ValueError("No time values >= {}".format(time_min))
    selected_time = time_values[selected_idx]

    output = Path(output).resolve()
    if output.exists():
        if not overwrite:
            raise FileExistsError("Output exists; use --overwrite: {}".format(output))
        output.unlink()
    output.parent.mkdir(parents=True, exist_ok=True)

    first_path = records[0][0]
    with nc.Dataset(first_path, "r") as first, nc.Dataset(output, "w", format="NETCDF4") as out:
        out.createDimension("time", len(selected_time))
        out.createDimension("lat", len(lat))
        out.createDimension("lon", len(lon))

        for name, src_var in first.variables.items():
            dims = tuple(src_var.dimensions)
            if dims not in {
                ("time",), ("lat",), ("lon",),
                ("time", "lat", "lon"), ("lat", "lon"),
            }:
                raise RuntimeError(
                    "Unsupported deterministic output variable dimensions {}: {}"
                    .format(name, dims)
                )
            _create_variable_like(
                out,
                name,
                src_var,
                compression=compression,
                time_block=time_block,
            )

        _copy_attrs(
            first,
            out,
            exclude=("band_id", "server_mode", "output_file", "history"),
        )
        out.setncattr("server_mode", 0)
        out.setncattr("merged_from_bands", int(len(records)))
        out.setncattr("merge_pattern", str(pattern))
        out.setncattr(
            "history",
            "Merged {} longitude bands at {}".format(
                len(records), datetime.now(timezone.utc).isoformat()
            ),
        )

        out.variables["time"][:] = selected_time
        out.variables["lat"][:] = lat
        out.variables["lon"][:] = lon.astype(out.variables["lon"].dtype)

        lon_offset = 0
        for band_number, (path, band_lon) in enumerate(records):
            band_width = len(band_lon)
            target_lon = slice(lon_offset, lon_offset + band_width)
            print(
                "[MERGE] band {}/{}: {} -> lon[{}:{}]".format(
                    band_number + 1,
                    len(records),
                    path.name,
                    lon_offset,
                    lon_offset + band_width,
                ),
                flush=True,
            )

            with nc.Dataset(path, "r") as src:
                for name, src_var in src.variables.items():
                    dims = tuple(src_var.dimensions)
                    if dims in {("time",), ("lat",), ("lon",)}:
                        continue
                    if dims == ("lat", "lon"):
                        out.variables[name][:, target_lon] = src_var[:, :]
                        continue
                    if dims == ("time", "lat", "lon"):
                        for t0 in range(0, len(selected_idx), max(1, int(time_block))):
                            t1 = min(len(selected_idx), t0 + max(1, int(time_block)))
                            source_indices = selected_idx[t0:t1]
                            out.variables[name][t0:t1, :, target_lon] = (
                                src_var[source_indices, :, :]
                            )
                        continue
                lon_offset += band_width
            out.sync()

    # Final integrity check.
    with nc.Dataset(output, "r") as ds:
        if "done" not in ds.variables or not np.all(np.asarray(ds.variables["done"][:]) == 1):
            raise RuntimeError("Merged output failed done==1 integrity check")
        if len(ds.dimensions["lon"]) != len(lon):
            raise RuntimeError("Merged longitude dimension is incomplete")
        if len(ds.dimensions["time"]) != len(selected_time):
            raise RuntimeError("Merged time dimension is incomplete")

    print("[OK] merged {} bands -> {}".format(len(records), output))
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--time-min", type=int, default=1850)
    parser.add_argument("--time-block", type=int, default=32)
    parser.add_argument("--compression", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if not (0 <= args.compression <= 9):
        raise ValueError("--compression must be between 0 and 9")
    merge_bands(
        pattern=args.pattern,
        output=Path(args.output),
        time_min=args.time_min,
        time_block=args.time_block,
        compression=args.compression,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
