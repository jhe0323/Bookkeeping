#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

from src.LULCCSimulator import LULCCSimulator

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "In_ncfile"
OUT_DIR = ROOT_DIR / "Out_ncfile"
CONFIG_DIR = ROOT_DIR / "config"
EXPERIMENT_DIR = CONFIG_DIR / "experiments"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STATE_PATH = DATA_DIR / "states_1deg.nc"
TRANS_PATH = DATA_DIR / "transitions_1deg.nc"
PFT_PATH = DATA_DIR / "PFTmap_orchidee_1deg.nc"
CONFIG_PATH = CONFIG_DIR / "config.yml"

EXPERIMENT_ALIASES = {
    "area": EXPERIMENT_DIR / "harvest_area.yml",
    "bio-strict": EXPERIMENT_DIR / "harvest_bio_strict.yml",
    "bio-forced": EXPERIMENT_DIR / "harvest_bio_forced.yml",
}


def _resolve_experiment(value: str) -> Path:
    alias_path = EXPERIMENT_ALIASES.get(value)
    path = alias_path if alias_path is not None else Path(value)
    if not path.is_absolute():
        path = ROOT_DIR / path
    path = path.resolve()
    if not path.exists():
        raise FileNotFoundError(f"Experiment configuration not found: {path}")
    return path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one Bookkeeping harvest experiment.",
    )
    parser.add_argument(
        "--experiment",
        default="area",
        help=(
            "Experiment alias (area, bio-strict, bio-forced) or a YAML path. "
            "Default: area"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    experiment_path = _resolve_experiment(args.experiment)

    lat_slice = None
    lon_slice = None

    simulator = LULCCSimulator(
        config_path=str(CONFIG_PATH),
        experiment_path=str(experiment_path),
        LULC_path=str(STATE_PATH),
        trans_path=str(TRANS_PATH),
        pft_path=str(PFT_PATH),
        lat_slice=lat_slice,
        lon_slice=lon_slice,
        area_unit="ha",
    )

    output_path = OUT_DIR / (
        f"summary_global_{simulator.params.experiment_name}.nc"
    )
    print(
        "[experiment] "
        f"name={simulator.params.experiment_name}, "
        f"mode={simulator.params.harvest_mode}, "
        f"use_bioh={simulator.params.harvest_use_bioh}, "
        f"allow_expansion={simulator.params.harvest_allow_expansion}"
    )

    simulator.run_simulation_grid(
        years=1173,
        out_nc=str(output_path),
        lat_slice=None,
        lon_slice=None,
        start_year_idx=0,
        sync_every=200,
    )


if __name__ == "__main__":
    main()
