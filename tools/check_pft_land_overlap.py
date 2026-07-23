#!/usr/bin/env python3
"""Diagnose LUH2 land cells that lack usable PFT fractions.

This version distinguishes:
1. raw PFT-mask mismatch;
2. cells resolvable by the configured nearest-neighbour radius;
3. cells still unresolved after applying the configured fallback policy;
4. cell-count and land-area-weighted mismatch fractions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from src.file_loader import FileLoader, resolve_year_index
from src.parameter_loader import ParameterLoader
from src.run_config import load_run_config


def _normalize_candidate(values) -> Optional[np.ndarray]:
    values = values.filled(np.nan) if hasattr(values, "filled") else values
    array = np.asarray(values, dtype=float).reshape(-1)
    cleaned = np.where(np.isfinite(array) & (array > 0.0), array, 0.0)
    total = float(cleaned.sum())
    if total <= 0.0:
        return None
    return cleaned / total


def _nearest_valid_cell(
    pft_data,
    year: int,
    i: int,
    j: int,
    nlat: int,
    nlon: int,
    max_radius: int,
) -> Tuple[Optional[int], Optional[Tuple[int, int]], Optional[np.ndarray]]:
    """Return nearest valid PFT cell within a Chebyshev search radius.

    Longitude wraps globally so cells near -180/180 degrees are handled
    consistently. Latitude does not wrap.
    """
    for radius in range(1, max_radius + 1):
        candidates: List[Tuple[int, int, int]] = []
        for di in range(-radius, radius + 1):
            ii = i + di
            if ii < 0 or ii >= nlat:
                continue
            for dj in range(-radius, radius + 1):
                if max(abs(di), abs(dj)) != radius:
                    continue
                jj = (j + dj) % nlon
                distance2 = di * di + dj * dj
                candidates.append((distance2, ii, jj))

        for _distance2, ii, jj in sorted(candidates):
            normalized = _normalize_candidate(pft_data.get_cell(year, ii, jj))
            if normalized is not None:
                return radius, (ii, jj), normalized

    return None, None, None


def _grid_area_ha(latitudes: Sequence[float], lon_res: float) -> np.ndarray:
    latitudes = np.asarray(latitudes, dtype=float)
    if latitudes.size >= 2:
        lat_res = abs(float(np.median(np.diff(latitudes))))
    else:
        raise ValueError("At least two latitude coordinates are required.")

    radius_m = 6_371_000.0
    dphi = np.deg2rad(lat_res)
    dlambda = np.deg2rad(float(lon_res))
    row_area_m2 = radius_m ** 2 * dphi * dlambda * np.cos(np.deg2rad(latitudes))
    return row_area_m2 / 1e4


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
        state_channels = [
            (name, name)
            for name in ("v", "s", "p", "c", "U")
            if name in info.variables
        ]
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
    row_area_ha = _grid_area_ha(meta.internal_lat, lon_res)

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
    raw_empty_cells = 0
    raw_empty_area_ha = 0.0
    total_land_area_ha = 0.0
    resolvable_cells = 0
    resolvable_area_ha = 0.0
    unresolved_cells = 0
    unresolved_area_ha = 0.0

    dominant_counts = Counter()
    radius_counts: Dict[int, int] = defaultdict(int)
    radius_area_ha: Dict[int, float] = defaultdict(float)

    raw_examples = []
    unresolved_examples = []

    max_radius = int(cfg.pft_nearest_search_radius)
    policy = str(cfg.pft_missing_cell_policy).lower()
    has_default = cfg.pft_default_index is not None

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

            cell_land_area_ha = float(row_area_ha[i]) * land_sum
            land_cells += 1
            total_land_area_ha += cell_land_area_ha

            normalized = _normalize_candidate(
                pft_data.get_cell(cfg.start_year, i, j)
            )
            if normalized is not None:
                dominant_counts[int(np.argmax(normalized)) + 1] += 1
                continue

            raw_empty_cells += 1
            raw_empty_area_ha += cell_land_area_ha

            if len(raw_examples) < args.examples:
                raw_examples.append(
                    (
                        i,
                        j,
                        float(meta.internal_lat[i]),
                        float(meta.internal_lon[j]),
                        land_sum,
                    )
                )

            resolved = False
            if policy == "nearest":
                radius, source_cell, source_pft = _nearest_valid_cell(
                    pft_data=pft_data,
                    year=cfg.start_year,
                    i=i,
                    j=j,
                    nlat=nlat,
                    nlon=nlon,
                    max_radius=max_radius,
                )
                if source_pft is not None and radius is not None:
                    resolved = True
                    resolvable_cells += 1
                    resolvable_area_ha += cell_land_area_ha
                    radius_counts[radius] += 1
                    radius_area_ha[radius] += cell_land_area_ha
                elif has_default:
                    resolved = True
                    resolvable_cells += 1
                    resolvable_area_ha += cell_land_area_ha
                    radius_counts[-1] += 1
                    radius_area_ha[-1] += cell_land_area_ha

            elif policy == "default" and has_default:
                resolved = True
                resolvable_cells += 1
                resolvable_area_ha += cell_land_area_ha
                radius_counts[-1] += 1
                radius_area_ha[-1] += cell_land_area_ha

            if not resolved:
                unresolved_cells += 1
                unresolved_area_ha += cell_land_area_ha
                if len(unresolved_examples) < args.examples:
                    unresolved_examples.append(
                        (
                            i,
                            j,
                            float(meta.internal_lat[i]),
                            float(meta.internal_lon[j]),
                            land_sum,
                        )
                    )

    raw_cell_fraction = raw_empty_cells / land_cells if land_cells else 0.0
    raw_area_fraction = (
        raw_empty_area_ha / total_land_area_ha
        if total_land_area_ha > 0.0
        else 0.0
    )
    unresolved_cell_fraction = (
        unresolved_cells / land_cells if land_cells else 0.0
    )
    unresolved_area_fraction = (
        unresolved_area_ha / total_land_area_ha
        if total_land_area_ha > 0.0
        else 0.0
    )

    print("config:", cfg.config_path)
    print("PFT variable:", pft_data.variable_name)
    print("PFT source mode:", pft_data.source_mode)
    print("PFT dynamic:", pft_data.dynamic)
    print("Missing-cell policy:", policy)
    print("Nearest search radius:", max_radius)
    print("Default PFT index:", cfg.pft_default_index)
    print()

    print("LUH2 land cells:", land_cells)
    print("Total LUH2 land area (ha):", total_land_area_ha)
    print("Raw land cells with empty PFT:", raw_empty_cells)
    print("Raw mismatch fraction by cell count:", raw_cell_fraction)
    print("Raw mismatch fraction by land area:", raw_area_fraction)
    print()

    print("Resolvable by configured fallback:", resolvable_cells)
    print("Resolvable land area (ha):", resolvable_area_ha)
    print("Still unresolved after configured fallback:", unresolved_cells)
    print("Unresolved fraction by cell count:", unresolved_cell_fraction)
    print("Unresolved fraction by land area:", unresolved_area_fraction)

    if radius_counts:
        print("Resolution counts by radius:")
        for radius in sorted(radius_counts):
            label = "default" if radius == -1 else str(radius)
            area_fraction = (
                radius_area_ha[radius] / total_land_area_ha
                if total_land_area_ha > 0.0
                else 0.0
            )
            print(
                "  radius {}: cells={}, land_area_ha={}, global_land_fraction={}"
                .format(
                    label,
                    radius_counts[radius],
                    radius_area_ha[radius],
                    area_fraction,
                )
            )

    print("Dominant PFT counts on originally valid land cells:",
          dict(sorted(dominant_counts.items())))

    if raw_examples:
        print()
        print("Raw mismatch examples: i, j, lat, lon, LUH2 land sum")
        for item in raw_examples:
            print("  ", item)

    if unresolved_examples:
        print()
        print("Unresolved examples after configured fallback:")
        for item in unresolved_examples:
            print("  ", item)

    if unresolved_cells:
        print()
        print(
            "RESULT: FAIL — at least one LUH2 land cell remains without a PFT "
            "after applying the configured fallback."
        )
        return 1

    print()
    print(
        "RESULT: PASS — all LUH2 land cells have a direct PFT or can be "
        "resolved by the configured fallback."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
