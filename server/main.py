#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import annotations

from pathlib import Path
from LULCCSimulator import LULCCSimulator


# ==================================================
# Paths
# ==================================================
SERVER_ROOT = Path("/mnt/beegfs/product/lulc0120")
DATA_DIR = SERVER_ROOT
OUT_DIR = SERVER_ROOT / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STATE_PATH = DATA_DIR / "states.nc"
TRANS_PATH = DATA_DIR / "transitions.nc"
PFT_PATH = DATA_DIR / "IBIS_PFT_dominant_0.25deg.nc"

BASE_DIR = Path(__file__).resolve().parent
CONFIG_PATH = BASE_DIR / "config.yml"


def output_complete(path: Path, expected_time: int, expected_lat: int, expected_lon: int) -> bool:
    """
    Check whether an output NetCDF file is complete.

    A file is considered complete only if:
    - it can be opened normally;
    - time/lat/lon dimensions match expectations;
    - it contains the 'done' variable;
    - all grid cells in done[:, :] are marked as 1.
    """
    if not path.exists():
        return False

    try:
        import netCDF4 as nc
        import numpy as np

        with nc.Dataset(path, "r") as ds:
            if "time" not in ds.dimensions:
                return False
            if "lat" not in ds.dimensions:
                return False
            if "lon" not in ds.dimensions:
                return False

            nt = len(ds.dimensions["time"])
            nlat = len(ds.dimensions["lat"])
            nlon = len(ds.dimensions["lon"])

            if nt != expected_time or nlat != expected_lat or nlon != expected_lon:
                return False

            if "done" not in ds.variables:
                return False

            done = np.asarray(ds.variables["done"][:])
            return int(np.nansum(done)) == nlat * nlon

    except Exception:
        return False


if __name__ == "__main__":
    import os
    import time
    import random
    import netCDF4 as nc
    import numpy as np

    # ==================================================
    # Simulation period
    # years means annual transitions.
    # T = years + 1.
    # 850-2020 inclusive gives 1171 time steps.
    # ==================================================
    start_year = 850
    end_year = 2020
    total_years = end_year - start_year
    expected_time = total_years + 1

    # ==================================================
    # Band settings
    # ==================================================
    lat_min, lat_max = -90.0, 90.0

    # BAND_SIZE=3 means 120 longitude bands globally at 3 degrees each.
    band_size_deg = float(os.environ.get("BAND_SIZE", "3"))

    # Prefer BAND_ID, fallback to SLURM_ARRAY_TASK_ID.
    band_id = int(
        os.environ.get(
            "BAND_ID",
            os.environ.get("SLURM_ARRAY_TASK_ID", "0")
        )
    )

    # --------------------------------------------------
    # Staggered start
    # Avoid all array tasks reading NetCDF files at the same time.
    # Set STAGGER_MAX=0 to disable.
    # --------------------------------------------------
    stagger_max = float(os.environ.get("STAGGER_MAX", "0"))

    if stagger_max > 0:
        random.seed(os.getpid() + band_id * 10007)
        sleep_time = random.uniform(0.0, stagger_max)
        print(
            f"[STAGGER] band_id={band_id}, sleep {sleep_time:.1f} s before loading input files",
            flush=True,
        )
        time.sleep(sleep_time)

    bands_total = int(round(360.0 / band_size_deg))
    if band_id < 0 or band_id >= bands_total:
        raise SystemExit(
            f"BAND_ID out of range: {band_id}, expected [0, {bands_total - 1}]"
        )

    lon_min = -180.0 + band_id * band_size_deg
    lon_max = lon_min + band_size_deg

    eps = 1e-9
    lon_max_adj = lon_max - eps

    # ==================================================
    # Read only coordinates to compute internal slices
    # ==================================================
    with nc.Dataset(str(STATE_PATH), "r") as ds:
        if "lat" in ds.variables:
            lat = np.asarray(ds.variables["lat"][:], dtype=float)
        elif "latitude" in ds.variables:
            lat = np.asarray(ds.variables["latitude"][:], dtype=float)
        else:
            raise ValueError("Cannot find lat / latitude in states.nc")

        if "lon" in ds.variables:
            lon = np.asarray(ds.variables["lon"][:], dtype=float)
        elif "longitude" in ds.variables:
            lon = np.asarray(ds.variables["longitude"][:], dtype=float)
        else:
            raise ValueError("Cannot find lon / longitude in states.nc")

    lat_res = abs(float(np.median(np.diff(lat))))
    lon_res = abs(float(np.median(np.diff(lon))))

    nlat = lat.size
    nlon = int(round(360.0 / lon_res))

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

    def wrap180(x):
        return ((x + 180.0) % 360.0) - 180.0

    lon_min_w = wrap180(lon_min)
    lon_max_w = wrap180(lon_max_adj)

    i0 = int(np.searchsorted(lat_internal, lat_min, side="left"))
    i1 = int(np.searchsorted(lat_internal, lat_max, side="right"))
    j0 = int(np.searchsorted(lon_internal, lon_min_w, side="left"))
    j1 = int(np.searchsorted(lon_internal, lon_max_w, side="right"))

    lat_slice = slice(i0, i1)
    lon_slice = slice(j0, j1)

    expected_lat = i1 - i0
    expected_lon = j1 - j0

    print(
        f"[MAIN] band_id={band_id}, "
        f"lon=({lon_min},{lon_max}), "
        f"lat_slice=({i0},{i1}), lon_slice=({j0},{j1}), "
        f"expected_shape=(time={expected_time}, lat={expected_lat}, lon={expected_lon})",
        flush=True,
    )

    # ==================================================
    # Output file
    # ==================================================
    out_nc = OUT_DIR / f"summary_025deg.rank{band_id:03d}.nc"

    if output_complete(out_nc, expected_time, expected_lat, expected_lon):
        print(f"[SKIP] Complete output found: {out_nc}", flush=True)
        raise SystemExit(0)

    if out_nc.exists():
        print(f"[RESTART] Incomplete or corrupted output found, removing: {out_nc}", flush=True)
        out_nc.unlink()

    # ==================================================
    # Create simulator with I/O slice
    # ==================================================
    print("[INIT] Creating LULCCSimulator...", flush=True)

    sim = LULCCSimulator(
        config_path=str(CONFIG_PATH),
        LULC_path=str(STATE_PATH),
        trans_path=str(TRANS_PATH),
        pft_path=str(PFT_PATH),
        lat_slice=lat_slice,
        lon_slice=lon_slice,
        area_unit="ha",
    )

    print("[INIT] LULCCSimulator created.", flush=True)

    # ==================================================
    # Run local band
    # The simulator has already loaded local data.
    # Therefore run_simulation_grid receives lat_slice=None/lon_slice=None.
    # ==================================================
    sim.run_simulation_grid(
        years=total_years,
        out_nc=str(out_nc),
        lat_slice=None,
        lon_slice=None,
        start_year_idx=0,
        sync_every=200,
    )

    print(f"[DONE] Finished band_id={band_id}: {out_nc}", flush=True)