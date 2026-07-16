# main.py
#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import annotations
import os
from pathlib import Path
from LULCCSimulator import LULCCSimulator

ROOT_DIR = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT_DIR / "In_ncfile"
OUT_DIR = ROOT_DIR / "Out_ncfile"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# path for input files
STATE_PATH = DATA_DIR / "states_1deg.nc"
# STATE_PATH = DATA_DIR / "synthetic_states.nc"

TRANS_PATH = DATA_DIR / "transitions_1deg.nc"
# TRANS_PATH = DATA_DIR / "synthetic_transitions.nc"

PFT_PATH = DATA_DIR / "PFTmap_orchidee_1deg.nc"
# PFT_PATH = DATA_DIR / "synthetic_pft.nc"

CONFIG_PATH = ROOT_DIR / "config" / "config.yml"

DENSITY_PATH = ROOT_DIR / "config" / "dynamic_carbon_density.parquet"

# Local output file
OUT_NC = OUT_DIR / "summary_global_win_serial.nc"


if __name__ == "__main__":
    # Local test mode:
    # - Keep lat_slice/lon_slice as None to load and run the full local grid.
    # - For quick tests, set small slices below, for example:
    #     lat_slice = slice(80, 100)
    #     lon_slice = slice(120, 140)
    lat_slice = None
    lon_slice = None

    sim = LULCCSimulator(
        config_path=str(CONFIG_PATH),
        LULC_path=str(STATE_PATH),
        trans_path=str(TRANS_PATH),
        pft_path=str(PFT_PATH),
        lat_slice=lat_slice,
        lon_slice=lon_slice,
        area_unit="ha",
    )

    # If lat_slice/lon_slice were already passed to LULCCSimulator, the loaded
    # dataset is local to that slice. Therefore run_simulation_grid usually
    # receives lat_slice=None and lon_slice=None here.
    
    # lat_slice, lon_slice = sim.indices_for_bbox(
    #     lat_min=3.0, 
    #     lat_max=54.0, 
    #     lon_min=73.0, 
    #     lon_max=136.0
    # )
    
    sim.run_simulation_grid(
        years=1173,
        out_nc=str(OUT_NC),
        lat_slice=None,
        lon_slice=None,
        start_year_idx=0,
        sync_every=200,
    )
    
    # sim.run_diagnostic_scan(lat_slice, lon_slice)
    
    # n_bands = 12
    # nlon = 360  # 1度数据通常是360；0.25度则是1440

    # for band in range(n_bands):
        # j0 = band * nlon // n_bands
        # j1 = (band + 1) * nlon // n_bands

        # lat_slice = None
        # lon_slice = slice(j0, j1)

        # sim = LULCCSimulator(
            # config_path=str(CONFIG_PATH),
            # LULC_path=str(STATE_PATH),
            # trans_path=str(TRANS_PATH),
            # pft_path=str(PFT_PATH),
            # lat_slice=lat_slice,
            # lon_slice=lon_slice,
            # area_unit="ha",
        # )

        # out_nc = OUT_DIR / f"summary_global_win_band_{band:02d}.nc"

        # sim.run_simulation_grid(
            # years=173,
            # out_nc=str(out_nc),
            # lat_slice=None,
            # lon_slice=None,
            # start_year_idx=1000,
            # sync_every=200,
        # )