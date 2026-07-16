#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Run a 0.25° longitude-band bookkeeping simulation on the server."""

from __future__ import annotations

import argparse
import os
import random
import sys
import time
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml


# ----------------------------------------------------------------------
# Directory layout
#
# /mnt/beegfs/product/lulc0120/
# ├── Code/                  <- Git repository
# │   ├── src/
# │   ├── config/
# │   └── server/            <- this file
# ├── In_ncfile/
# ├── Out_ncfile/
# ├── logs/
# └── site-packages/
# ----------------------------------------------------------------------
SERVER_DIR = Path(__file__).resolve().parent
CODE_DIR = SERVER_DIR.parent
PROJECT_DIR = CODE_DIR.parent

# Make direct execution and module execution both able to import src.
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from src.LULCCSimulator import LULCCSimulator  # noqa: E402


DATA_DIR = PROJECT_DIR / "In_ncfile"
OUT_ROOT = PROJECT_DIR / "Out_ncfile"
CONFIG_DIR = CODE_DIR / "config"
EXPERIMENT_DIR = CONFIG_DIR / "experiments"
CONFIG_PATH = CONFIG_DIR / "config.yml"

EXPERIMENT_ALIASES = {
    "area": EXPERIMENT_DIR / "harvest_area.yml",
    "bio-strict": EXPERIMENT_DIR / "harvest_bio_strict.yml",
    "bio-forced": EXPERIMENT_DIR / "harvest_bio_forced.yml",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one 0.25° longitude band on the server."
    )
    parser.add_argument(
        "--experiment",
        default=os.environ.get("EXPERIMENT", "area"),
        help=(
            "Experiment alias (area, bio-strict, bio-forced) or YAML path. "
            "Default: environment variable EXPERIMENT or area."
        ),
    )
    parser.add_argument("--state-file", default="states.nc")
    parser.add_argument("--transition-file", default="transitions.nc")
    parser.add_argument("--pft-file", default="IBIS_PFT_dominant_0.25deg.nc")
    parser.add_argument("--start-year", type=int, default=850)
    parser.add_argument("--end-year", type=int, default=2020)
    parser.add_argument("--sync-every", type=int, default=200)
    return parser.parse_args()


def resolve_experiment(value: str) -> Path:
    alias = EXPERIMENT_ALIASES.get(value)
    path = alias if alias is not None else Path(value)
    if not path.is_absolute():
        path = CODE_DIR / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Experiment configuration not found: {path}")
    return path


def experiment_name(path: Path) -> str:
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    name = cfg.get("experiment", {}).get("name", path.stem)
    return str(name).strip().replace(" ", "_")


def resolve_input(value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = DATA_DIR / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    return path


def coord_to_1d(arr: np.ndarray, name: str) -> np.ndarray:
    """Convert a regular 1-D or 2-D coordinate variable to a 1-D axis."""
    arr = np.asarray(arr, dtype=float)

    if arr.ndim == 1:
        return arr

    if arr.ndim == 2:
        if name.lower().startswith("lat"):
            col0 = arr[:, 0]
            if not np.allclose(arr, col0[:, None], equal_nan=True):
                raise ValueError(f"Irregular 2-D latitude coordinate: {arr.shape}")
            return col0

        if name.lower().startswith("lon"):
            row0 = arr[0, :]
            if not np.allclose(arr, row0[None, :], equal_nan=True):
                raise ValueError(f"Irregular 2-D longitude coordinate: {arr.shape}")
            return row0

    raise ValueError(f"Unsupported {name} coordinate shape: {arr.shape}")


def read_coordinates(state_path: Path) -> tuple[np.ndarray, np.ndarray]:
    with nc.Dataset(str(state_path), "r") as ds:
        lat_name = "lat" if "lat" in ds.variables else "latitude"
        lon_name = "lon" if "lon" in ds.variables else "longitude"

        if lat_name not in ds.variables:
            raise ValueError(f"Cannot find lat/latitude in {state_path}")
        if lon_name not in ds.variables:
            raise ValueError(f"Cannot find lon/longitude in {state_path}")

        lat = coord_to_1d(ds.variables[lat_name][:], "lat")
        lon = coord_to_1d(ds.variables[lon_name][:], "lon")

    if lat.size < 2 or lon.size < 2:
        raise ValueError(f"Invalid coordinate lengths: lat={lat.size}, lon={lon.size}")

    return lat, lon


def output_complete(
    path: Path,
    expected_time: int,
    expected_lat: int,
    expected_lon: int,
) -> bool:
    """Return True only when dimensions are correct and all done flags equal 1."""
    if not path.exists():
        return False

    try:
        with nc.Dataset(path, "r") as ds:
            for name in ("time", "lat", "lon"):
                if name not in ds.dimensions:
                    return False

            actual = (
                len(ds.dimensions["time"]),
                len(ds.dimensions["lat"]),
                len(ds.dimensions["lon"]),
            )
            expected = (expected_time, expected_lat, expected_lon)
            if actual != expected:
                return False

            if "done" not in ds.variables:
                return False

            done = np.asarray(ds.variables["done"][:])
            return done.shape == (expected_lat, expected_lon) and bool(
                np.all(done == 1)
            )
    except Exception as exc:
        print(f"[CHECK] Cannot validate {path}: {exc}", flush=True)
        return False


def main() -> None:
    args = parse_args()

    state_path = resolve_input(args.state_file)
    transition_path = resolve_input(args.transition_file)
    pft_path = resolve_input(args.pft_file)
    experiment_path = resolve_experiment(args.experiment)
    exp_name = experiment_name(experiment_path)

    if not CONFIG_PATH.is_file():
        raise FileNotFoundError(f"Base configuration not found: {CONFIG_PATH}")
    if args.end_year <= args.start_year:
        raise ValueError("end-year must be greater than start-year")

    total_years = args.end_year - args.start_year
    expected_time = total_years + 1

    band_size_deg = float(os.environ.get("BAND_SIZE", "3"))
    band_id = int(
        os.environ.get(
            "BAND_ID",
            os.environ.get("SLURM_ARRAY_TASK_ID", "0"),
        )
    )
    stagger_max = float(os.environ.get("STAGGER_MAX", "0"))

    if band_size_deg <= 0 or 360.0 % band_size_deg != 0:
        raise ValueError(f"BAND_SIZE must divide 360 exactly: {band_size_deg}")

    bands_total = int(round(360.0 / band_size_deg))
    if band_id < 0 or band_id >= bands_total:
        raise SystemExit(
            f"BAND_ID out of range: {band_id}; expected 0-{bands_total - 1}"
        )

    if stagger_max > 0:
        random.seed(os.getpid() + band_id * 10007)
        sleep_time = random.uniform(0.0, stagger_max)
        print(
            f"[STAGGER] band_id={band_id}; sleeping {sleep_time:.1f} seconds",
            flush=True,
        )
        time.sleep(sleep_time)

    lat, lon = read_coordinates(state_path)
    lat_res = abs(float(np.median(np.diff(lat))))
    lon_res = abs(float(np.median(np.diff(lon))))
    nlat = lat.size
    nlon = lon.size

    expected_global_nlon = int(round(360.0 / lon_res))
    if nlon != expected_global_nlon:
        raise ValueError(
            f"Longitude axis is not global/regular: nlon={nlon}, "
            f"expected={expected_global_nlon}, resolution={lon_res}"
        )

    lat_internal = np.linspace(
        -90.0 + 0.5 * lat_res,
        90.0 - 0.5 * lat_res,
        nlat,
    )
    lon_internal = np.linspace(
        -180.0 + 0.5 * lon_res,
        180.0 - 0.5 * lon_res,
        nlon,
    )

    lon_min = -180.0 + band_id * band_size_deg
    lon_max = lon_min + band_size_deg
    lon_max_adjusted = lon_max - 1e-9

    def wrap180(x: float) -> float:
        return ((x + 180.0) % 360.0) - 180.0

    i0 = int(np.searchsorted(lat_internal, -90.0, side="left"))
    i1 = int(np.searchsorted(lat_internal, 90.0, side="right"))
    j0 = int(np.searchsorted(lon_internal, wrap180(lon_min), side="left"))
    j1 = int(
        np.searchsorted(lon_internal, wrap180(lon_max_adjusted), side="right")
    )

    i0 = max(0, min(i0, nlat))
    i1 = max(0, min(i1, nlat))
    j0 = max(0, min(j0, nlon))
    j1 = max(0, min(j1, nlon))

    if i1 <= i0 or j1 <= j0:
        raise ValueError(
            f"Empty band slice: lat=({i0},{i1}), lon=({j0},{j1}), "
            f"band_id={band_id}"
        )

    lat_slice = slice(i0, i1)
    lon_slice = slice(j0, j1)
    expected_lat = i1 - i0
    expected_lon = j1 - j0

    out_dir = OUT_ROOT / "025deg" / exp_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_nc = out_dir / f"summary_025deg.rank{band_id:03d}.nc"

    print(
        "[MAIN] "
        f"experiment={exp_name}, band_id={band_id}/{bands_total - 1}, "
        f"lon=({lon_min},{lon_max}), "
        f"lat_slice=({i0},{i1}), lon_slice=({j0},{j1}), "
        f"expected_shape=({expected_time},{expected_lat},{expected_lon})",
        flush=True,
    )
    print(f"[PATH] states={state_path}", flush=True)
    print(f"[PATH] transitions={transition_path}", flush=True)
    print(f"[PATH] pft={pft_path}", flush=True)
    print(f"[PATH] output={out_nc}", flush=True)

    if output_complete(out_nc, expected_time, expected_lat, expected_lon):
        print(f"[SKIP] Complete output found: {out_nc}", flush=True)
        return

    if out_nc.exists():
        print(f"[RESTART] Removing incomplete output: {out_nc}", flush=True)
        out_nc.unlink()

    print("[INIT] Creating LULCCSimulator...", flush=True)
    simulator = LULCCSimulator(
        config_path=str(CONFIG_PATH),
        experiment_path=str(experiment_path),
        LULC_path=str(state_path),
        trans_path=str(transition_path),
        pft_path=str(pft_path),
        lat_slice=lat_slice,
        lon_slice=lon_slice,
        area_unit="ha",
    )

    simulator.run_simulation_grid(
        years=total_years,
        out_nc=str(out_nc),
        lat_slice=None,
        lon_slice=None,
        start_year_idx=0,
        sync_every=args.sync_every,
    )

    if not output_complete(out_nc, expected_time, expected_lat, expected_lon):
        raise RuntimeError(f"Output was written but did not pass completion check: {out_nc}")

    print(f"[DONE] Finished band_id={band_id}: {out_nc}", flush=True)


if __name__ == "__main__":
    main()
