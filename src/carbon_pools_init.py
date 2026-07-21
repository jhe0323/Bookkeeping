"""Carbon-pool array definitions and initialization helpers."""
from __future__ import annotations

from typing import Optional
import numpy as np

pool_index = {"B": 0, "SS": 1, "SR": 2, "P1": 3, "P10": 4, "P100": 5, "A": 6}
pool_index_Cbar = {"B": 0, "SS": 1}
pool_index_Delta = dict(pool_index)
pool_index_active = {"B": 0, "SS": 1, "SR": 2}
pool_index_release = {"B": 0, "SS": 1}

cover_index = {"v": 0, "s": 1, "p": 2, "c": 3, "U": 4}
LULC_keys = list(cover_index.keys())
history_index = {"l": 0, "h": 1, "a": 2, "g": 3}

N_POOL_CBAR = len(pool_index_Cbar)
N_POOL_DELTA = len(pool_index_Delta)
N_POOL_ACTIVE = len(pool_index_active)
N_POOL_RELEASE = len(pool_index_release)
N_COVER = len(cover_index)
N_HISTORY = len(history_index)


def normalize_pft_fractions(pft_grid, n_pft: int) -> np.ndarray:
    fractions = np.asarray(pft_grid, dtype=np.float64).reshape(-1).copy()
    if fractions.size != int(n_pft):
        raise ValueError(
            f"Expected {n_pft} PFT fractions, got {fractions.size}. "
            "The PFT map and parameter configuration must use the same order."
        )
    fractions[~np.isfinite(fractions)] = 0.0
    fractions[fractions < 0.0] = 0.0
    total = float(fractions.sum())
    if total <= 0.0:
        raise ValueError("A simulated land cell has no valid PFT fraction.")
    fractions /= total
    return fractions


def resolve_n_time(*, n_time: Optional[int] = None, years: Optional[int] = None) -> int:
    if n_time is not None:
        if n_time <= 0:
            raise ValueError(f"n_time must be positive, got {n_time}.")
        return int(n_time)
    if years is not None:
        if years < 0:
            raise ValueError(f"years must be non-negative, got {years}.")
        return int(years) + 1
    raise ValueError("Either n_time or years must be provided.")


def initialize_Cbar(
    initial_land_frac,
    cell_area,
    pft_grid,
    params,
    *,
    n_time=None,
    years=None,
) -> np.ndarray:
    n_time_resolved = resolve_n_time(n_time=n_time, years=years)
    n_pft = int(params.n_pft)
    pft_fractions = normalize_pft_fractions(pft_grid, n_pft)
    C_bar = np.zeros(
        (n_time_resolved, N_POOL_CBAR, N_COVER, n_pft),
        dtype=np.float64,
    )

    for land, frac_land in initial_land_frac.items():
        if land not in cover_index or frac_land <= 0:
            continue
        j = cover_index[land]
        if land == "U":
            area = float(frac_land) * float(cell_area)
            C_bar[0, 0, j, 0] = area * params.get_carbon_density(0, "Biomass", "U")
            C_bar[0, 1, j, 0] = area * params.get_carbon_density(0, "Soil", "U")
            continue

        for p in np.flatnonzero(pft_fractions > 0.0):
            area = float(frac_land) * pft_fractions[p] * float(cell_area)
            C_bar[0, 0, j, p] = area * params.get_carbon_density(p, "Biomass", land)
            C_bar[0, 1, j, p] = area * params.get_carbon_density(p, "Soil", land)
    return C_bar


def initialize_Delta(
    *,
    n_pft: int,
    n_time: Optional[int] = None,
    years: Optional[int] = None,
) -> np.ndarray:
    n_time_resolved = resolve_n_time(n_time=n_time, years=years)
    return np.zeros(
        (n_time_resolved, N_POOL_DELTA, N_COVER, N_HISTORY, int(n_pft)),
        dtype=np.float64,
    )


def initialize_frac_area(initial_land_frac, pft_grid, *, n_pft: int) -> np.ndarray:
    pft_fractions = normalize_pft_fractions(pft_grid, n_pft)
    frac_area = np.zeros((N_COVER, int(n_pft)), dtype=np.float64)
    for land, frac_land in initial_land_frac.items():
        if land not in cover_index:
            continue
        frac_land = float(frac_land)
        if not np.isfinite(frac_land) or frac_land <= 0.0:
            continue
        j = cover_index[land]
        if land == "U":
            frac_area[j, 0] = frac_land
        else:
            frac_area[j, :] = frac_land * pft_fractions
    return frac_area


def initialize_dC_bar(*, n_pft: int, n_time=None, years=None) -> np.ndarray:
    return np.zeros(
        (resolve_n_time(n_time=n_time, years=years), N_POOL_CBAR, N_COVER, int(n_pft)),
        dtype=np.float64,
    )


def initialize_dDelta(*, n_pft: int, n_time=None, years=None) -> np.ndarray:
    return np.zeros(
        (resolve_n_time(n_time=n_time, years=years), N_POOL_ACTIVE, N_COVER, N_HISTORY, int(n_pft)),
        dtype=np.float64,
    )


def initialize_dDelta_sum(*, n_pft: int, n_time=None, years=None) -> np.ndarray:
    return np.zeros(
        (resolve_n_time(n_time=n_time, years=years), N_POOL_RELEASE, N_COVER, int(n_pft)),
        dtype=np.float64,
    )


def initialize_dC_released(*, n_pft: int, n_time=None, years=None) -> np.ndarray:
    return np.zeros(
        (resolve_n_time(n_time=n_time, years=years), N_POOL_RELEASE, N_COVER, int(n_pft)),
        dtype=np.float64,
    )


def refresh_Cbar_from_frac_area(C_bar, t_idx, frac_area, cell_area, params) -> None:
    C_bar[t_idx, :, :, :] = 0.0
    for land, j in cover_index.items():
        if land == "U":
            continue
        for p in range(int(params.n_pft)):
            area = float(frac_area[j, p]) * float(cell_area)
            if area <= 0.0:
                continue
            C_bar[t_idx, pool_index_Cbar["B"], j, p] = (
                area * params.get_carbon_density(p, "Biomass", land)
            )
            C_bar[t_idx, pool_index_Cbar["SS"], j, p] = (
                area * params.get_carbon_density(p, "Soil", land)
            )
