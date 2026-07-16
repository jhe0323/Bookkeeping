from typing import Dict, Optional
from parameter_loader import ParameterLoader
import numpy as np

# pool types (i): vegetation biomass (B), 
# soil carbon undergoing slow relaxation processes (SS) or rapid relaxation processes (SR),
# atmospheric carbon / cumulative emissions (A), and product pools in which products 
# decompose on an average time scale of 1 year (P1), 10 years (P10), or 100 years (P100)
pool_index = {"B":0, "SS":1, "SR":2, "P1":3, "P10":4, "P100":5, "A":6} 
# Pools used by equilibrium carbon C_bar[i, j, l]
pool_index_Cbar = {"B": 0, "SS": 1}
# Pools tracked in Delta[i, j, k, l]
pool_index_Delta = {"B": 0, "SS": 1, "SR": 2, "P1": 3, "P10": 4, "P100": 5, "A": 6}
# Pools that can actively change during event bookkeeping steps
pool_index_active = {"B": 0, "SS": 1, "SR": 2}
# Pools used when computing total released carbon during transitions
pool_index_release = {"B": 0, "SS": 1}
# land cover types (j): primary / virgin land (v), secondary land (s), pasture (p), and crop (c). 
# U for urban, currently not supported
cover_index = {"v":0, "s":1, "p":2, "c":3, "U":4}
LULC_keys = list(cover_index.keys())
# history types (k): after clearing (clearing, l), harvest (harvest, h), 
# abandonment (abandon, a), and without or after other transitions (others, g) 
history_index = {"l":0, "h":1, "a":2, "g":3}
# PFT index l: 15 natural PFTs, aligned with LUCE
N_PFT = 15
pft_index = {i + 1: i for i in range(N_PFT)}
# Dimensions
N_POOL_CBAR = len(pool_index_Cbar)
N_POOL_DELTA = len(pool_index_Delta)
N_POOL_ACTIVE = len(pool_index_active)
N_POOL_RELEASE = len(pool_index_release)
N_COVER = len(cover_index)
N_HISTORY = len(history_index)

def resolve_n_time(*, n_time: Optional[int] = None, years: Optional[int] = None) -> int:
    """
    Resolve the size of the time dimension used by pool arrays.
    Parameters
    ----------
    n_time : int, optional
        Explicit number of stored model states. If provided, it is used directly.
    years : int, optional
        Number of simulated annual transitions. If ``n_time`` is not provided,
        the stored states are assumed to include the initial condition plus one
        state after each simulated year, i.e. ``years + 1``.
    Returns
    -------
    int
        Length of the time dimension.
    Notes
    -----
    ``start_year_idx`` from main.py should *not* be included in ``n_time``.
    It is only an index offset into the external forcing files. For model pool
    storage, the correct length is usually ``years + 1``.
    """
    if n_time is not None:
        if n_time <= 0:
            raise ValueError(f"n_time must be positive, got {n_time}.")
        return int(n_time)

    if years is not None:
        if years < 0:
            raise ValueError(f"years must be non-negative, got {years}.")
        return int(years) + 1

    raise ValueError("Either n_time or years must be provided.")

# Initilize carbon pools (biomass and soil per LULC, and product, atmosphere)
def initialize_Cbar(
    initial_land_frac, 
    cell_area, 
    pft_grid, 
    params,
    *, 
    n_time=None, 
    years=None
) -> np.ndarray:
    """
    Initialize equilibrium carbon pools C_bar[t, i, j, l]
    Returns: np.ndarray
        Equilibrium carbon array with shape (n_time, N_POOL_CBAR, N_COVER, N_PFT).
      Notes:
      * i: B, SS
      * j: v, s, p, c
      * l: 1, 2, ..., 11
      * At initialization, assume the system starts at equilibrium:
          Cbar filled from initial areas and carbon densities.
    """
    n_time_resolved = resolve_n_time(n_time=n_time, years=years)
    C_bar = np.zeros((n_time_resolved, N_POOL_CBAR, N_COVER, N_PFT), dtype=float)

    for land, frac_land in initial_land_frac.items():
        if land not in cover_index:
            continue
        j = cover_index[land]
        if frac_land <= 0:
            continue
        # The current code base keeps an optional urban class U for compatibility.
        # Its densities can be handled separately if the parameter file supports it.
        if land == "U":
            area = frac_land * cell_area
            try:
                rho_B = params.get_carbon_density(0, "Biomass", "U")
                rho_S = params.get_carbon_density(0, "Soil", "U")
            except:
                rho_B = rho_S = 0.0

            C_bar[0, 0, j, 0] = area * rho_B
            C_bar[0, 1, j, 0] = area * rho_S
            continue

        for p in range(N_PFT):
            frac_lp = frac_land * pft_grid[p]
            if frac_lp <= 0:
                continue

            area = frac_lp * cell_area

            rho_B = params.get_carbon_density(p, "Biomass", land)
            rho_S = params.get_carbon_density(p, "Soil", land)

            C_bar[0, 0, j, p] = area * rho_B
            C_bar[0, 1, j, p] = area * rho_S

    return C_bar
    
def initialize_Delta(*, n_time: Optional[int] = None, years: Optional[int] = None) -> np.ndarray:
    """
    Initialize excess carbon pools Delta[t, i, j, k, l]
    Returns: np.ndarray
    Excess carbon array with shape (n_time, N_POOL_DELTA, N_COVER, N_HISTORY, N_PFT)
    Notes:
      * i: B, SS, SR, P1, P10, P100, A
      * j: v, s, p, c
      * k: l, h, a, g
      * l: 1, 2, ..., 11
    """
    n_time_resolved = resolve_n_time(n_time=n_time, years=years)
    return np.zeros((n_time_resolved, N_POOL_DELTA, N_COVER, N_HISTORY, N_PFT), dtype=float)

def initialize_frac_area(initial_land_frac, pft_grid) -> np.ndarray:
    """
    Initialize fractional area array frac_area[j, l]
    Returns: np.ndarray
    Array with shape (N_COVER, N_PFT)
    Notes:
      * j: v, s, p, c
      * l: 1, 2, ..., 11
    """
    frac_area = np.zeros((N_COVER, N_PFT), dtype=float)
    for land, frac_land in initial_land_frac.items():
        if land not in cover_index:
            continue
        j = cover_index[land]

        if land == "U":
            frac_area[j, 0] = frac_land
            continue

        for p in range(N_PFT):
            frac_area[j, p] = frac_land * pft_grid[p]

    return frac_area

def initialize_dC_bar(*, n_time: Optional[int] = None, years: Optional[int] = None) -> np.ndarray:
    """
    Initialize transition-time changes in equilibrium carbon dC_bar[t, i, j, l]
    Returns: np.ndarray
    Array with shape (n_time, N_POOL_CBAR, N_COVER, N_PFT)
    Notes:
      * i: B, SS
      * j: v, s, p, c
      * l: 1, 2, ..., 11
    """
    n_time_resolved = resolve_n_time(n_time=n_time, years=years)
    return np.zeros((n_time_resolved, N_POOL_CBAR, N_COVER, N_PFT), dtype=float)

def initialize_dDelta(*, n_time: Optional[int] = None, years: Optional[int] = None) -> np.ndarray:
    """
    Initialize per-event changes in active excess pools dDelta[t, i, j, k, l]
    Returns: np.ndarray
    Array with shape (n_time, N_POOL_ACTIVE, N_COVER, N_HISTORY, N_PFT)
    Notes:
      * i: B, SS
      * j: v, s, p, c
      * k: l, h, a, g
      * l: 1, 2, ..., 11
    """
    n_time_resolved = resolve_n_time(n_time=n_time, years=years)
    return np.zeros((n_time_resolved, N_POOL_ACTIVE, N_COVER, N_HISTORY, N_PFT), dtype=float)

def initialize_dDelta_sum(*, n_time: Optional[int] = None, years: Optional[int] = None) -> np.ndarray:
    """
    Initialize history-summed excess-pool changes dDelta_sum[t, i, j, l]
    Returns: np.ndarray
    Array with shape (n_time, N_POOL_RELEASE, N_COVER, N_PFT)
    Notes:
      * i: B, SS
      * j: v, s, p, c
      * l: 1, 2, ..., 11
    """
    n_time_resolved = resolve_n_time(n_time=n_time, years=years)
    return np.zeros((n_time_resolved, N_POOL_RELEASE, N_COVER, N_PFT), dtype=float)

def initialize_dC_released(*, n_time: Optional[int] = None, years: Optional[int] = None) -> np.ndarray:
    """
    Initialize released-carbon bookkeeping array dC_released[t, i, j, l]
    Returns: np.ndarray
    Array with shape (n_time, N_POOL_RELEASE, N_COVER, N_PFT)
    Notes:
      * i: B, SS
      * j: v, s, p, c
      * l: 1, 2, ..., 11
    """
    n_time_resolved = resolve_n_time(n_time=n_time, years=years)
    return np.zeros((n_time_resolved, N_POOL_RELEASE, N_COVER, N_PFT), dtype=float)

def refresh_Cbar_from_frac_area(
    C_bar,
    t_idx,
    frac_area,
    cell_area,
    params,
):
    """
    Rebuild equilibrium C_bar at one time step using current frac_area
    and current-year carbon density from params.
    """
    C_bar[t_idx, :, :, :] = 0.0

    for land, j in cover_index.items():
        if land == "U":
            continue

        for p in range(N_PFT):
            area = float(frac_area[j, p]) * float(cell_area)
            if area <= 0.0:
                continue

            rho_B = params.get_carbon_density(p, "Biomass", land)
            rho_S = params.get_carbon_density(p, "Soil", land)

            C_bar[t_idx, pool_index_Cbar["B"], j, p] = area * rho_B
            C_bar[t_idx, pool_index_Cbar["SS"], j, p] = area * rho_S