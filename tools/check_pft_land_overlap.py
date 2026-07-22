#!/usr/bin/env python3
"""Check LUH2 land cells that have no usable PFT fractions."""
from __future__ import annotations

import argparse
from collections import Counter

import numpy as np

from src.file_loader import FileLoader, resolve_year_index
from src.parameter_loader import ParameterLoader
from src.run_config import load_run_config


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/run_1deg.yml")
    parser.add_argument("--examples", type=int, default=20)
    args = parser.parse_args()

    cfg = load_run_config(args.config)
    params = ParameterLoader(
        str(cfg.base_parameter_path),
        experiment_path=str(cfg.experiment_path),
        override_config=cfg.raw.get("model_overrides", {}) or {},
    )
    loader = FileLoader()
    info = loader.inspect_dataset(cfg.state_path)
    state_idx = resolve_year_index(
        info,
        cfg.start_year,
        input_base_year=cfg.input_base_year,
        time_encoding=cfg.time_encoding,
    )

    if cfg.input_format == "vscp":
        state_channels = [(name, name) for name in ("v", "s", "p", "c", "U") if name in info.variables]
    else:
        mapping = params.config.get("LUH2toLULC", {}) or {}
        state_channels = [
            (name, str(cover))
            for name, cover in mapping.items()
            if name in info.variables
        ]

    state_data, meta = loader.load_luh2_dataset(
        str(cfg.state_path),
        time_slice=slice(state_idx, state_idx + 1),
        variable_names=sorted({name for name, _ in state_channels}),
    )
    lon_res = abs(float(np.median(np.diff(meta.internal_lon))))
    pft_data = loader.load_pft_map(
        str(cfg.pft_path),
        pft_var=cfg.pft_variable,
        target_lat_asc=meta.internal_lat,
        target_lon_180=meta.internal_lon,
        target_lon_res=lon_res,
        expected_n_pft=params.n_pft,
        model_years=[cfg.start_year],
        pft_base_year=cfg.pft_base_year,
        time_encoding=cfg.pft_time_encoding,
        min_year_policy=cfg.pft_min_year_policy,
        max_year_policy=cfg.pft_max_year_policy,
    )

    nlat = len(meta.internal_lat)
    nlon = len(meta.internal_lon)
    land_cells = 0
    empty_pft_land = 0
    examples = []
    dominant_counts = Counter()

    for i in range(nlat):
        ii = nlat - 1 - i if meta.lat_desc else i
        for j in range(nlon):
            land_sum = 0.0
            for name, _cover in state_channels:
                value = state_data.variables[name][0, ii, j]
                value = value.filled(np.nan) if hasattr(value, "filled") else value
                value = float(value)
                if np.isfinite(value) and value > 0.0:
                    land_sum += value
            if land_sum <= 0.0:
                continue

            land_cells += 1
            raw = pft_data.get_cell(cfg.start_year, i, j)
            raw = raw.filled(np.nan) if hasattr(raw, "filled") else raw
            raw = np.asarray(raw, dtype=float).reshape(-1)
            clean = np.where(np.isfinite(raw) & (raw > 0.0), raw, 0.0)
            total = float(clean.sum())
            if total <= 0.0:
                empty_pft_land += 1
                if len(examples) < args.examples:
                    examples.append((
                        i,
                        j,
                        float(meta.internal_lat[i]),
                        float(meta.internal_lon[j]),
                        land_sum,
                    ))
            else:
                dominant_counts[int(np.argmax(clean)) + 1] += 1

    print("config:", cfg.config_path)
    print("PFT variable:", pft_data.variable_name)
    print("PFT source mode:", pft_data.source_mode)
    print("PFT dynamic:", pft_data.dynamic)
    print("LUH2 land cells:", land_cells)
    print("Land cells with empty PFT:", empty_pft_land)
    print("Mismatch fraction:", empty_pft_land / land_cells if land_cells else 0.0)
    print("Dominant PFT counts:", dict(sorted(dominant_counts.items())))
    if examples:
        print("Examples: i, j, lat, lon, LUH2 land sum")
        for item in examples:
            print("  ", item)

    return 1 if empty_pft_land else 0


if __name__ == "__main__":
    raise SystemExit(main())
