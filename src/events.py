"""Event carbon allocation models.

This module contains small, pure functions that allocate carbon pools for
different LULCC event types. They were previously implemented as private
methods on LULCCSimulator and are extracted here to keep the simulator core
lean and easier to maintain.

All functions take a `params: ParameterLoader` instance explicitly.
"""

from __future__ import annotations

from typing import Dict, Optional

from parameter_loader import ParameterLoader
import numpy as np
from carbon_pools_init import (
    cover_index,
    history_index,
    pft_index,
    pool_index_Cbar,
    pool_index_Delta,
)

def _diag_add(diag: Optional[dict], key: str, value: float) -> None:
    """Add to diagnostics only when the key is requested by the caller."""
    if diag is not None and key in diag:
        diag[key] += float(value)

def simulate_clearing(
    loss_biomass: float,
    *,
    src: str,
    p: int,
    params: ParameterLoader
) -> Dict[str, float]:
    """Partition the removed biomass from a clearing event.

    Parameters
    ----------
    loss_biomass : float
        Biomass carbon removed from the cleared patch (tC).
    src : str
        Source land-cover type: 'v' or 's'.
    p : int
        0-based pft index.
    params : ParameterLoader

    Returns
    -------
    dict with keys:
        A, P1, P10, P100, SS, SR
    """
    if loss_biomass <= 0.0:
        return {"A": 0.0, "P1": 0.0, "P10": 0.0, "P100": 0.0, "SS": 0.0, "SR": 0.0}
    
    bio_to_soil_frac = float(params.get_clearing_param(p, "Bio_Soil"))
    prod1_frac = float(params.get_clearing_param(p, "Prod_1"))
    prod10_frac = float(params.get_clearing_param(p, "Prod_10"))
    prod100_frac = float(params.get_clearing_param(p, "Prod_100"))
    f_fast = float(params.get_clearing_param(p, f"f_{src}"))

    # all fractions are applied to total removed biomass
    to_soil_total = loss_biomass * bio_to_soil_frac
    to_soil_SR = to_soil_total * f_fast
    to_soil_SS = to_soil_total * (1.0 - f_fast)

    P1 = loss_biomass * prod1_frac
    P10 = loss_biomass * prod10_frac
    P100 = loss_biomass * prod100_frac

    allocated = to_soil_total + P1 + P10 + P100
    A = loss_biomass - allocated

    # numerical safety
    if A < 0.0 and abs(A) < 1e-12:
        A = 0.0

    return {
        "A": float(A),
        "P1": float(P1),
        "P10": float(P10),
        "P100": float(P100),
        "SS": float(to_soil_SS),
        "SR": float(to_soil_SR),
    }

def apply_clearing(
    *,
    area: float,
    src: str,
    dst: str,
    cell_area: float,
    pft_grid,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    t_idx: int = 0,
    diag: Optional[dict] = None,
) -> None:
    """
    BLUE clearing: v/s -> p/c

    Key BLUE behavior:
    - source B and SS are removed proportionally across all histories
    - source SR does NOT participate in the land transfer
    - target biomass is assumed cleared immediately:
        Delta_B(dst, clearing) -= B_eq_dst
    - released biomass is partitioned to A / P1 / P10 / P100 / SR / SS
    - soil transfer uses only source current SS, not legacy SR
    """
    if area <= 0.0:
        return
    if src not in ("v", "s"):
        return
    if dst not in ("p", "c"):
        return
    if cell_area <= 0.0:
        return

    j_src = cover_index[src]
    j_dst = cover_index[dst]
    k_clr = history_index["l"]

    src_frac_total = float(LULC_frac.get(src, 0.0))
    if src_frac_total <= 0.0:
        return

    src_area_total = src_frac_total * cell_area
    if src_area_total <= 0.0:
        return

    if area > src_area_total:
        area = src_area_total

    frac_src_total = area / src_area_total
    if frac_src_total <= 0.0:
        return

    # update whole-grid cover fractions
    LULC_frac[src] -= area / cell_area
    LULC_frac[dst] += area / cell_area
    if LULC_frac[src] < 0.0:
        LULC_frac[src] = 0.0

    for p in range(len(pft_grid)):
        a_src_p = float(frac_area[j_src, p])
        if a_src_p <= 0.0:
            continue

        src_area_p_total = a_src_p * cell_area
        if src_area_p_total <= 0.0:
            continue

        area_p = src_area_p_total * frac_src_total
        if area_p <= 0.0:
            continue
        if area_p > src_area_p_total:
            area_p = src_area_p_total

        frac_src_p = area_p / src_area_p_total
        if frac_src_p <= 0.0:
            continue

        # source current carbon before area removal
        B_src_cur_total = float(
            C_bar[t_idx, pool_index_Cbar["B"], j_src, p]
            + Delta[t_idx, pool_index_Delta["B"], j_src, :, p].sum()
        )
        SS_src_cur_total = float(
            C_bar[t_idx, pool_index_Cbar["SS"], j_src, p]
            + Delta[t_idx, pool_index_Delta["SS"], j_src, :, p].sum()
        )

        # current carbon on the cleared patch
        B_removed = B_src_cur_total * frac_src_p
        SS_removed = SS_src_cur_total * frac_src_p

        # biomass partitioning after clearing
        part = simulate_clearing(B_removed, src=src, p=p, params=params)

        # 1) source excess pools reduced proportionally across ALL histories
        #    NOTE: SR does NOT participate in clearing transfer (BLUE 3.3.2)
        Delta[t_idx, pool_index_Delta["B"], j_src, :, p] *= (1.0 - frac_src_p)
        Delta[t_idx, pool_index_Delta["SS"], j_src, :, p] *= (1.0 - frac_src_p)

        # 2) move equilibrium pools according to changed area
        B_eq_src = area_p * float(params.get_carbon_density(p, "Biomass", src))
        SS_eq_src = area_p * float(params.get_carbon_density(p, "Soil", src))

        B_eq_dst = area_p * float(params.get_carbon_density(p, "Biomass", dst))
        SS_eq_dst = area_p * float(params.get_carbon_density(p, "Soil", dst))

        C_bar[t_idx, pool_index_Cbar["B"], j_src, p] -= B_eq_src
        C_bar[t_idx, pool_index_Cbar["SS"], j_src, p] -= SS_eq_src
        C_bar[t_idx, pool_index_Cbar["B"], j_dst, p] += B_eq_dst
        C_bar[t_idx, pool_index_Cbar["SS"], j_dst, p] += SS_eq_dst

        # 3) target biomass is immediately cleared (BLUE Eq.15)
        Delta[t_idx, pool_index_Delta["B"], j_dst, k_clr, p] -= B_eq_dst

        # 4) soil transfer:
        #    only source current SS participates; the mismatch to target equilibrium
        #    is split into rapid/slow fractions using clearing fSR,j
        f_fast = float(params.get_clearing_param(p, f"f_{src}"))
        soil_mismatch = SS_removed - SS_eq_dst

        Delta[t_idx, pool_index_Delta["SR"], j_dst, k_clr, p] += soil_mismatch * f_fast
        Delta[t_idx, pool_index_Delta["SS"], j_dst, k_clr, p] += soil_mismatch * (1.0 - f_fast)

        # plus litter/dead biomass transfer from clearing
        Delta[t_idx, pool_index_Delta["SR"], j_dst, k_clr, p] += float(part["SR"])
        Delta[t_idx, pool_index_Delta["SS"], j_dst, k_clr, p] += float(part["SS"])

        # 5) products and atmosphere
        Delta[t_idx, pool_index_Delta["P1"], j_dst, k_clr, p] += float(part["P1"])
        Delta[t_idx, pool_index_Delta["P10"], j_dst, k_clr, p] += float(part["P10"])
        Delta[t_idx, pool_index_Delta["P100"], j_dst, k_clr, p] += float(part["P100"])
        Delta[t_idx, pool_index_Delta["A"], j_dst, k_clr, p] += float(part["A"])

        if diag is not None:
            prod_sum = float(part["P1"] + part["P10"] + part["P100"])
            soil_sum = float(part["SR"] + part["SS"] + soil_mismatch)

            _diag_add(diag, "area_clearing", area_p)
            _diag_add(diag, "clearing_to_A", float(part["A"]))
            _diag_add(diag, "clearing_to_products", prod_sum)
            _diag_add(diag, "clearing_to_soil", soil_sum)

            if _is_forest_pft(p, params):
                _diag_add(diag, "area_deforestation", area_p)
                _diag_add(diag, "deforestation_to_A", float(part["A"]))
                _diag_add(diag, "deforestation_to_products", prod_sum)
                _diag_add(diag, "deforestation_to_soil", soil_sum)
        # 6) update fractional area
        a_removed_p = area_p / cell_area
        frac_area[j_src, p] -= a_removed_p
        frac_area[j_dst, p] += a_removed_p
        if frac_area[j_src, p] < 0.0:
            frac_area[j_src, p] = 0.0


def apply_abandonment(
    *,
    area: float,
    src: str,
    dst: str,
    cell_area: float,
    pft_grid,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    t_idx: int = 0,
    diag: Optional[dict] = None,
) -> None:
    """
    BLUE abandonment: c/p -> s

    Key BLUE behavior:
    - no direct allocation to atmosphere/products
    - target biomass is NOT cleared
    - source current biomass is inherited into destination
    - soil transfer uses only source current SS
    - source SR does NOT participate in the land transfer
    - changes are written into abandonment history 'a'
    """
    if area <= 0.0:
        return
    if cell_area <= 0.0:
        return
    if src not in ("c", "p"):
        return
    if dst != "s":
        return

    j_src = cover_index[src]
    j_dst = cover_index[dst]
    k_abd = history_index["a"]

    src_frac_total = float(LULC_frac.get(src, 0.0))
    if src_frac_total <= 0.0:
        return

    src_area_total = src_frac_total * cell_area
    if src_area_total <= 0.0:
        return

    if area > src_area_total:
        area = src_area_total

    frac_src_total = area / src_area_total
    if frac_src_total <= 0.0:
        return

    # update whole-grid cover fractions
    LULC_frac[src] -= area / cell_area
    LULC_frac[dst] += area / cell_area
    if LULC_frac[src] < 0.0:
        LULC_frac[src] = 0.0

    for p in range(len(pft_grid)):
        a_src_p = float(frac_area[j_src, p])
        if a_src_p <= 0.0:
            continue

        src_area_p_total = a_src_p * cell_area
        if src_area_p_total <= 0.0:
            continue

        area_p = src_area_p_total * frac_src_total
        if area_p <= 0.0:
            continue
        if area_p > src_area_p_total:
            area_p = src_area_p_total

        frac_src_p = area_p / src_area_p_total
        if frac_src_p <= 0.0:
            continue

        # source current carbon before area removal
        B_src_cur_total = float(
            C_bar[t_idx, pool_index_Cbar["B"], j_src, p]
            + Delta[t_idx, pool_index_Delta["B"], j_src, :, p].sum()
        )
        SS_src_cur_total = float(
            C_bar[t_idx, pool_index_Cbar["SS"], j_src, p]
            + Delta[t_idx, pool_index_Delta["SS"], j_src, :, p].sum()
        )

        # inherited current carbon from source patch
        B_inherit = B_src_cur_total * frac_src_p
        SS_inherit = SS_src_cur_total * frac_src_p

        # remove source excess proportionally across all histories
        # NOTE: SR does NOT participate in abandonment transfer
        Delta[t_idx, pool_index_Delta["B"], j_src, :, p] *= (1.0 - frac_src_p)
        Delta[t_idx, pool_index_Delta["SS"], j_src, :, p] *= (1.0 - frac_src_p)

        # move equilibrium pools by changed area
        B_eq_src = area_p * float(params.get_carbon_density(p, "Biomass", src))
        SS_eq_src = area_p * float(params.get_carbon_density(p, "Soil", src))

        B_eq_dst = area_p * float(params.get_carbon_density(p, "Biomass", dst))
        SS_eq_dst = area_p * float(params.get_carbon_density(p, "Soil", dst))

        C_bar[t_idx, pool_index_Cbar["B"], j_src, p] -= B_eq_src
        C_bar[t_idx, pool_index_Cbar["SS"], j_src, p] -= SS_eq_src
        C_bar[t_idx, pool_index_Cbar["B"], j_dst, p] += B_eq_dst
        C_bar[t_idx, pool_index_Cbar["SS"], j_dst, p] += SS_eq_dst

        # BLUE Eq.23 / Eq.24 style
        Delta[t_idx, pool_index_Delta["B"], j_dst, k_abd, p] += (B_inherit - B_eq_dst)
        Delta[t_idx, pool_index_Delta["SS"], j_dst, k_abd, p] += (SS_inherit - SS_eq_dst)

        if diag is not None:
            _diag_add(diag, "area_abandonment", area_p)
            _diag_add(diag, "abd_delta_b", (B_inherit - B_eq_dst))
            _diag_add(diag, "abd_delta_ss", (SS_inherit - SS_eq_dst))
        # update fractional area
        a_removed_p = area_p / cell_area
        frac_area[j_src, p] -= a_removed_p
        frac_area[j_dst, p] += a_removed_p
        if frac_area[j_src, p] < 0.0:
            frac_area[j_src, p] = 0.0

def apply_others(
    *,
    area: float,
    src: str,
    dst: str,
    cell_area: float,
    pft_grid,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    t_idx: int = 0,
    diag: Optional[dict] = None,
) -> None:
    """
    BLUE other transitions: c <-> p

    Key BLUE behavior:
    - no direct allocation to atmosphere/products
    - target biomass is NOT cleared
    - source current biomass is inherited into destination
    - soil transfer uses only source current SS
    - source SR does NOT participate in the land transfer
    - changes are written into general history 'g'
    """
    if area <= 0.0:
        return
    if cell_area <= 0.0:
        return
    if src not in ("c", "p"):
        return
    if dst not in ("c", "p"):
        return
    if src == dst:
        return

    j_src = cover_index[src]
    j_dst = cover_index[dst]
    k_oth = history_index["g"]

    src_frac_total = float(LULC_frac.get(src, 0.0))
    if src_frac_total <= 0.0:
        return

    src_area_total = src_frac_total * cell_area
    if src_area_total <= 0.0:
        return

    if area > src_area_total:
        area = src_area_total

    frac_src_total = area / src_area_total
    if frac_src_total <= 0.0:
        return

    # update whole-grid cover fractions
    LULC_frac[src] -= area / cell_area
    LULC_frac[dst] += area / cell_area
    if LULC_frac[src] < 0.0:
        LULC_frac[src] = 0.0

    for p in range(len(pft_grid)):
        a_src_p = float(frac_area[j_src, p])
        if a_src_p <= 0.0:
            continue

        src_area_p_total = a_src_p * cell_area
        if src_area_p_total <= 0.0:
            continue

        area_p = src_area_p_total * frac_src_total
        if area_p <= 0.0:
            continue
        if area_p > src_area_p_total:
            area_p = src_area_p_total

        frac_src_p = area_p / src_area_p_total
        if frac_src_p <= 0.0:
            continue

        # source current carbon before area removal
        B_src_cur_total = float(
            C_bar[t_idx, pool_index_Cbar["B"], j_src, p]
            + Delta[t_idx, pool_index_Delta["B"], j_src, :, p].sum()
        )
        SS_src_cur_total = float(
            C_bar[t_idx, pool_index_Cbar["SS"], j_src, p]
            + Delta[t_idx, pool_index_Delta["SS"], j_src, :, p].sum()
        )

        B_inherit = B_src_cur_total * frac_src_p
        SS_inherit = SS_src_cur_total * frac_src_p

        # remove source excess proportionally across all histories
        # NOTE: SR does NOT participate in other-transition transfer
        Delta[t_idx, pool_index_Delta["B"], j_src, :, p] *= (1.0 - frac_src_p)
        Delta[t_idx, pool_index_Delta["SS"], j_src, :, p] *= (1.0 - frac_src_p)

        # move equilibrium pools
        B_eq_src = area_p * float(params.get_carbon_density(p, "Biomass", src))
        SS_eq_src = area_p * float(params.get_carbon_density(p, "Soil", src))

        B_eq_dst = area_p * float(params.get_carbon_density(p, "Biomass", dst))
        SS_eq_dst = area_p * float(params.get_carbon_density(p, "Soil", dst))

        C_bar[t_idx, pool_index_Cbar["B"], j_src, p] -= B_eq_src
        C_bar[t_idx, pool_index_Cbar["SS"], j_src, p] -= SS_eq_src
        C_bar[t_idx, pool_index_Cbar["B"], j_dst, p] += B_eq_dst
        C_bar[t_idx, pool_index_Cbar["SS"], j_dst, p] += SS_eq_dst

        # BLUE Eq.23 / Eq.24 style, history = g
        Delta[t_idx, pool_index_Delta["B"], j_dst, k_oth, p] += (B_inherit - B_eq_dst)
        Delta[t_idx, pool_index_Delta["SS"], j_dst, k_oth, p] += (SS_inherit - SS_eq_dst)

        if diag is not None:
            _diag_add(diag, "area_other", area_p)
        # update fractional area
        a_removed_p = area_p / cell_area
        frac_area[j_src, p] -= a_removed_p
        frac_area[j_dst, p] += a_removed_p
        if frac_area[j_src, p] < 0.0:
            frac_area[j_src, p] = 0.0

HARVEST_TYPE_MAP = {
    # LUH2 wood harvest classes
    "primf": {"src_cover": "v", "veg_group": "forest"},
    "primn": {"src_cover": "v", "veg_group": "nonforest"},
    "secmf": {"src_cover": "s", "veg_group": "forest"},
    "secyf": {"src_cover": "s", "veg_group": "forest"},
    "secnf": {"src_cover": "s", "veg_group": "nonforest"},
}

def _is_forest_pft(p: int, params) -> bool:
    """
    Forest PFTs are defined by config.yml: forest_pft_indices.
    The config uses 1-based PFT indices, while p is 0-based in code.
    """
    return p in params.get_forest_pfts()


def _pft_allowed_for_harvest_type(p: int, harvest_type: str, params) -> bool:
    if harvest_type not in HARVEST_TYPE_MAP:
        raise ValueError(f"Unknown harvest_type: {harvest_type}")

    veg_group = HARVEST_TYPE_MAP[harvest_type]["veg_group"]
    is_forest = _is_forest_pft(p, params)

    if veg_group == "forest":
        return is_forest
    if veg_group == "nonforest":
        return not is_forest

    raise ValueError(f"Unsupported veg_group: {veg_group}")


def _compute_patch_current_biomass(
    *,
    t_idx: int,
    j_src: int,
    p: int,
    C_bar,
    Delta,
) -> float:
    return float(
        C_bar[t_idx, pool_index_Cbar["B"], j_src, p]
        + Delta[t_idx, pool_index_Delta["B"], j_src, :, p].sum()
    )


def _compute_patch_current_slow_soil(
    *,
    t_idx: int,
    j_src: int,
    p: int,
    C_bar,
    Delta,
) -> float:
    return float(
        C_bar[t_idx, pool_index_Cbar["SS"], j_src, p]
        + Delta[t_idx, pool_index_Delta["SS"], j_src, :, p].sum()
    )

def simulate_harvest(
    loss_biomass: float,
    *,
    src: str,
    p: int,
    params: ParameterLoader,
) -> Dict[str, float]:
    """
    Partition harvested biomass carbon following BLUE harvest logic.

    Parameters
    ----------
    loss_biomass : float
        Harvested biomass carbon on the affected patch (tC).
    src : str
        Source cover type: 'v' or 's'.
    p : int
        0-based pft index.

    Returns
    -------
    dict with keys:
        P1, P10, P100, SR
    Notes
    -----
    BLUE 3.3.4:
    - a fraction hS,j of harvested biomass goes to rapid soil pool SR
    - the rest goes to product pools P1/P10/P100
    - no direct "event-time" atmosphere term is added here
      (product/SR decay is handled at end-of-time-step relaxation)
    """
    if loss_biomass <= 0.0:
        return {"P1": 0.0, "P10": 0.0, "P100": 0.0, "SR": 0.0}

    # Map BLUE src symbols to your config naming if needed:
    # here we assume harvest params are keyed by v/s style names consistently
    hS = float(params.get_harvest_param(p, f"Bio_Soil_{src}"))
    hP1 = float(params.get_harvest_param(p, "Prod_1"))
    hP10 = float(params.get_harvest_param(p, "Prod_10"))
    hP100 = float(params.get_harvest_param(p, "Prod_100"))

    to_soil_SR = loss_biomass * hS
    to_products = loss_biomass * (1.0 - hS)

    P1 = to_products * hP1
    P10 = to_products * hP10
    P100 = to_products * hP100

    # numerical safety: keep mass balance inside products if fractions imperfect
    s_prod = P1 + P10 + P100
    if s_prod > to_products and s_prod > 0.0:
        fac = to_products / s_prod
        P1 *= fac
        P10 *= fac
        P100 *= fac

    return {
        "P1": float(P1),
        "P10": float(P10),
        "P100": float(P100),
        "SR": float(to_soil_SR),
    }

def apply_harvest(
    *,
    area: float,
    src: str,
    cell_area: float,
    pft_grid,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    biomass_harv_tC: Optional[float] = None,
    t_idx: int = 0,
) -> None:
    """
    BLUE harvest:
      - primary harvest:   v -> s
      - secondary harvest: s -> s

    Parameters
    ----------
    area : float
        Harvested area in the current grid cell (absolute area, e.g. ha).
    src : str
        'v' or 's'
    cell_area : float
        Grid-cell area.
    pft_grid : np.ndarray
        PFT fractions of the grid cell.
    biomass_harv_tC : Optional[float]
        Total harvested biomass reported by LUH2 for this source family in this cell-year.
        If given, it is distributed across PFT patches proportionally to current patch biomass.

    Notes
    -----
    BLUE 3.3.4 behavior implemented here:
    - target cover is always secondary land ('s')
    - source excess B is NOT reduced history-by-history as in clearing Eq.11
      to avoid re-attributing abandonment regrowth to harvest
    - instead, current biomass on harvested patch is computed and added as
      a new harvest-history biomass deficit on target secondary land
    - biomass partition:
        * fraction to SR
        * rest to product pools
    - soil:
        * compute current source SS on patch
        * any amount above minimum soil carbon after harvest goes to SR
        * slow-soil harvest history gets:
              -R + Cbar_src_soil_removed - Cbar_dst_soil_added
          which is BLUE Eq.33
    - no decay is applied here; decay must happen at end of timestep
    """
    if area <= 0.0 or cell_area <= 0.0 or src not in ("v", "s"):
        return

    j_src = cover_index[src]
    j_dst = cover_index["s"]
    k_har = history_index["h"]

    src_frac_total = float(LULC_frac.get(src, 0.0))
    if src_frac_total <= 0.0:
        return

    src_area_total = src_frac_total * cell_area
    if src_area_total <= 0.0:
        return

    area = min(area, src_area_total)
    frac_src_total = area / src_area_total

    # Primary harvest (v -> s) involves a total land cover area transfer
    if src == "v":
        LULC_frac["v"] -= area / cell_area
        LULC_frac["s"] += area / cell_area
        if LULC_frac["v"] < 0.0:
            LULC_frac["v"] = 0.0

    for p in range(len(pft_grid)):
        a_src_p = float(frac_area[j_src, p])
        if a_src_p <= 0.0:
            continue

        src_area_p_total = a_src_p * cell_area
        if src_area_p_total <= 0.0:
            continue

        area_p = min(src_area_p_total * frac_src_total, src_area_p_total)
        if area_p <= 0.0:
            continue
            
        frac_src_p = area_p / src_area_p_total

        # ==========================================
        # 1. Calculate carbon corresponding to equilibrium pool area (Eq. 10 & 14)
        # ==========================================
        B_eq_src = area_p * float(params.get_carbon_density(p, "Biomass", src))
        SS_eq_src = area_p * float(params.get_carbon_density(p, "Soil", src))

        B_eq_dst = area_p * float(params.get_carbon_density(p, "Biomass", "s"))
        SS_eq_dst = area_p * float(params.get_carbon_density(p, "Soil", "s"))

        # ==========================================
        # 2. Biomass calculation and partitioning (Eq. 25-29)
        # ==========================================
        # Eq. 25: Extract current Excess Biomass for the affected area
        beta_delta_sum = float(Delta[t_idx, pool_index_Delta["B"], j_src, :, p].sum()) * frac_src_p
        
        # Eq. 27: Calculate total harvested biomass carbon (beta_j_l)
        beta_j_l = B_eq_src + beta_delta_sum

        # Eq. 26: Record the new harvest history deficit in target secondary land (s)
        # Note: We subtract the target equilibrium carbon plus the source's Excess carbon
        Delta[t_idx, pool_index_Delta["B"], j_dst, k_har, p] -= (B_eq_dst + beta_delta_sum)

        # Eq. 28 & 29: Partition the harvested biomass to product pools (P1, P10, P100) and rapid soil (SR)
        part = simulate_harvest(beta_j_l, src=src, p=p, params=params)
        Delta[t_idx, pool_index_Delta["SR"], j_dst, k_har, p] += float(part["SR"])
        Delta[t_idx, pool_index_Delta["P1"], j_dst, k_har, p] += float(part["P1"])
        Delta[t_idx, pool_index_Delta["P10"], j_dst, k_har, p] += float(part["P10"])
        Delta[t_idx, pool_index_Delta["P100"], j_dst, k_har, p] += float(part["P100"])

        # ==========================================
        # 3. Soil carbon calculation and partitioning (Eq. 30-33)
        # ==========================================
        # Eq. 30 & 31: Calculate current total soil carbon in the affected area (sigma_j_l)
        sigma_delta_sum = float(Delta[t_idx, pool_index_Delta["SS"], j_src, :, p].sum()) * frac_src_p
        sigma_j_l = SS_eq_src + sigma_delta_sum

        # Eq. 32: Calculate soil carbon above the post-harvest minimum threshold (assigned to rapid decay)
        soc_min_density = float(params.get_harvest_param(p, f"SOC_min_{src}"))
        soil_min_patch = area_p * soc_min_density
        R = max(0.0, sigma_j_l - soil_min_patch)

        # Eq. 32 (cont.): Add R to the rapid soil pool (SR)
        Delta[t_idx, pool_index_Delta["SR"], j_dst, k_har, p] += R

        # Eq. 33: Update the harvest history deficit for slow soil (SS)
        Delta[t_idx, pool_index_Delta["SS"], j_dst, k_har, p] += (-R + SS_eq_src - SS_eq_dst)

        # ==========================================
        # 4. Area and equilibrium pool transfer (v -> s only)
        # ==========================================
        # Note: s -> s harvest does not change total area; legacy carbon history is automatically retained in the 's' category of the original grid
        if src == "v":
            # Move equilibrium pool C_bar
            C_bar[t_idx, pool_index_Cbar["B"], j_src, p] -= B_eq_src
            C_bar[t_idx, pool_index_Cbar["SS"], j_src, p] -= SS_eq_src
            C_bar[t_idx, pool_index_Cbar["B"], j_dst, p] += B_eq_dst
            C_bar[t_idx, pool_index_Cbar["SS"], j_dst, p] += SS_eq_dst

            # Adjust area fractions for different PFTs
            a_removed_p = area_p / cell_area
            frac_area[j_src, p] -= a_removed_p
            frac_area[j_dst, p] += a_removed_p
            if frac_area[j_src, p] < 0.0:
                frac_area[j_src, p] = 0.0
def _allocate_by_capacity(target_tC: float, capacity_by_pft: dict) -> tuple[dict, float]:
    """
    按可用 biomass capacity 分配 target。
    返回:
      take_by_pft: 每个 PFT 实际分配的 removal
      unmet: capacity 不足时剩余未满足量
    """
    if target_tC <= 0.0:
        return {p: 0.0 for p in capacity_by_pft}, 0.0

    caps = {p: max(0.0, float(v)) for p, v in capacity_by_pft.items()}
    total_cap = sum(caps.values())

    if total_cap <= 0.0:
        return {p: 0.0 for p in caps}, target_tC

    take_total = min(target_tC, total_cap)
    take_by_pft = {
        p: take_total * cap / total_cap
        for p, cap in caps.items()
    }
    unmet = max(0.0, target_tC - take_total)
    return take_by_pft, unmet


def _apply_harvest_loss_to_patch(
    *,
    p: int,
    src: str,
    loss_biomass: float,
    effective_area_p: float,
    cell_area: float,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    t_idx: int,
    diag: Optional[dict] = None,
) -> dict:
    """
    对单个 PFT/source 执行真实 biomass removal。

    注意:
    - loss_biomass 是总采伐 biomass，不是产品量。
    - effective_area_p 是为了实现该 loss 对应的本模型有效采伐面积。
    - v -> s 时更新 C_bar 和 frac_area。
    - s -> s 时不改变总面积，只重置 secondary biomass/soil history。
    """
    if effective_area_p <= 0.0:
        return {
            "products": 0.0,
            "sr_biomass": 0.0,
            "sr_soil": 0.0,
            "delta_B_h": 0.0,
            "delta_SS_h": 0.0,
            "beta_delta_sum": 0.0,
            "sigma_delta_sum": 0.0,
        }

    # Area-driven harvest may encounter very low or zero current biomass.
    # Keep the area/soil/cover bookkeeping active, but do not allow negative removal.
    loss_biomass = max(0.0, float(loss_biomass))

    j_src = cover_index[src]
    j_dst = cover_index["s"]
    k_har = history_index["h"]

    src_area_p_total = float(frac_area[j_src, p]) * cell_area
    if src_area_p_total <= 0.0:
        return {
            "products": 0.0,
            "sr_biomass": 0.0,
            "sr_soil": 0.0,
            "delta_B_h": 0.0,
            "delta_SS_h": 0.0,
        }

    effective_area_p = min(effective_area_p, src_area_p_total)
    intensity_p = effective_area_p / src_area_p_total
    intensity_p = max(0.0, min(1.0, intensity_p))

    # ===== Biomass equilibrium pools on affected area =====
    B_eq_src = effective_area_p * float(params.get_carbon_density(p, "Biomass", src))
    B_eq_dst = effective_area_p * float(params.get_carbon_density(p, "Biomass", "s"))

    # source excess biomass on affected area
    beta_delta_sum_raw = float(
        Delta[t_idx, pool_index_Delta["B"], j_src, :, p].sum()
    ) * intensity_p

    # remove source biomass history proportionally
    Delta[t_idx, pool_index_Delta["B"], j_src, :, p] *= (1.0 - intensity_p)

    # write biomass deficit to secondary harvest history
    delta_B_h = beta_delta_sum_raw + (B_eq_src - B_eq_dst) - loss_biomass
    Delta[t_idx, pool_index_Delta["B"], j_dst, k_har, p] += delta_B_h

    # partition removed biomass to products and rapid soil
    part = simulate_harvest(loss_biomass, src=src, p=p, params=params)

    Delta[t_idx, pool_index_Delta["SR"], j_dst, k_har, p] += float(part["SR"])
    Delta[t_idx, pool_index_Delta["P1"], j_dst, k_har, p] += float(part["P1"])
    Delta[t_idx, pool_index_Delta["P10"], j_dst, k_har, p] += float(part["P10"])
    Delta[t_idx, pool_index_Delta["P100"], j_dst, k_har, p] += float(part["P100"])

    products = float(part["P1"] + part["P10"] + part["P100"])
    sr_biomass = float(part["SR"])

    # ===== Soil disturbance on affected area =====
    SS_eq_src = effective_area_p * float(params.get_carbon_density(p, "Soil", src))
    SS_eq_dst = effective_area_p * float(params.get_carbon_density(p, "Soil", "s"))

    sigma_delta_sum_raw = float(
        Delta[t_idx, pool_index_Delta["SS"], j_src, :, p].sum()
    ) * intensity_p

    sigma_j_l = SS_eq_src + sigma_delta_sum_raw

    soc_min_density = float(params.get_harvest_param(p, f"SOC_min_{src}"))
    soil_min_patch = effective_area_p * soc_min_density
    R = max(0.0, sigma_j_l - soil_min_patch)

    Delta[t_idx, pool_index_Delta["SS"], j_src, :, p] *= (1.0 - intensity_p)

    delta_SS_h = sigma_delta_sum_raw + (SS_eq_src - SS_eq_dst) - R

    Delta[t_idx, pool_index_Delta["SR"], j_dst, k_har, p] += R
    Delta[t_idx, pool_index_Delta["SS"], j_dst, k_har, p] += delta_SS_h

    # ===== v -> s area and C_bar transfer =====
    if src == "v":
        C_bar[t_idx, pool_index_Cbar["B"], j_src, p] -= B_eq_src
        C_bar[t_idx, pool_index_Cbar["SS"], j_src, p] -= SS_eq_src
        C_bar[t_idx, pool_index_Cbar["B"], j_dst, p] += B_eq_dst
        C_bar[t_idx, pool_index_Cbar["SS"], j_dst, p] += SS_eq_dst

        a_removed_p = effective_area_p / cell_area
        frac_area[j_src, p] -= a_removed_p
        frac_area[j_dst, p] += a_removed_p

        if frac_area[j_src, p] < 0.0:
            frac_area[j_src, p] = 0.0

    return {
        "products": products,
        "sr_biomass": sr_biomass,
        "sr_soil": float(R),
        "delta_B_h": float(delta_B_h),
        "delta_SS_h": float(delta_SS_h),
        "beta_delta_sum": float(beta_delta_sum_raw),
        "sigma_delta_sum": float(sigma_delta_sum_raw),
    }
    
def apply_harvest_luh2(
    *,
    harvest_type: str,
    area_frac: float,
    biomass_harv_kgC: Optional[float] = None,
    cell_area: float,
    pft_grid,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    t_idx: int = 0,
    diag: Optional[dict] = None,
):
    """
    Pure LUH2-area-driven harvest.

    Core rules
    ----------
    - Use only LUH2 *_harv area fraction as the harvest forcing.
    - Ignore LUH2 *_bioh completely. The argument is kept only for backward
      compatibility with older callers.
    - Harvested biomass is diagnosed from the model state on the LUH2 harvest
      footprint:
          removed biomass = current modeled biomass × harvested area fraction
      within the corresponding source cover / allowed PFT group.
    - No harvest expansion outside the LUH2 footprint is allowed.
      Therefore forced_biomass, extra_area and unmet_biomass are always zero
      unless no valid area exists, in which case no harvest is applied.
    """

    # ===== 0. Basic checks and reported LUH2 area =====
    if cell_area <= 0.0:
        return 0.0, 0.0

    area_frac = float(area_frac) if np.isfinite(area_frac) else 0.0
    area_frac = max(0.0, min(1.0, area_frac))

    luh2_area_reported = area_frac * float(cell_area)
    _diag_add(diag, "harvest_luh2_area_frac", area_frac)
    _diag_add(diag, "harvest_luh2_area", luh2_area_reported)

    if area_frac <= 0.0 or luh2_area_reported <= 0.0:
        return 0.0, 0.0

    if harvest_type not in HARVEST_TYPE_MAP:
        return 0.0, 0.0

    meta = HARVEST_TYPE_MAP[harvest_type]
    src = meta["src_cover"]
    j_src = cover_index[src]

    src_frac_total = float(LULC_frac.get(src, 0.0))
    if src_frac_total <= 0.0:
        return 0.0, 0.0

    # ===== 1. Allowed PFT source area =====
    allowed_pfts = [
        p for p in range(len(pft_grid))
        if _pft_allowed_for_harvest_type(p, harvest_type, params)
    ]
    if not allowed_pfts:
        return 0.0, 0.0

    src_area_p_map = {}
    B_patch_map = {}
    allowed_src_area_total = 0.0

    for p in allowed_pfts:
        a_src_p = float(frac_area[j_src, p])
        if a_src_p <= 0.0:
            continue

        src_area_p_total = a_src_p * float(cell_area)
        if src_area_p_total <= 0.0:
            continue

        B_patch_cur = _compute_patch_current_biomass(
            t_idx=t_idx,
            j_src=j_src,
            p=p,
            C_bar=C_bar,
            Delta=Delta,
        )

        src_area_p_map[p] = src_area_p_total
        B_patch_map[p] = float(B_patch_cur)
        allowed_src_area_total += src_area_p_total

    if allowed_src_area_total <= 0.0 or not src_area_p_map:
        return 0.0, 0.0

    # Apply only the portion of LUH2 area that can be located within the current
    # source cover and the harvest-type-compatible PFT group.
    target_area_total = min(luh2_area_reported, allowed_src_area_total)
    if target_area_total <= 0.0:
        return 0.0, 0.0

    # ===== 2. Apply harvest by LUH2 area, distributed across allowed PFTs =====
    total_effective_area = 0.0
    total_products = 0.0
    total_sr_biomass = 0.0
    total_sr_soil = 0.0
    total_delta_B_h = 0.0
    total_delta_SS_h = 0.0
    total_beta_delta_sum = 0.0
    total_sigma_delta_sum = 0.0
    total_R = 0.0
    total_patch_biomass_before = 0.0
    actual_removed_tC = 0.0

    for p, src_area_p_total in src_area_p_map.items():
        # LUH2 harvest area is allocated by available source area, not by bioh.
        effective_area_p = target_area_total * src_area_p_total / allowed_src_area_total
        effective_area_p = min(effective_area_p, src_area_p_total)
        if effective_area_p <= 0.0:
            continue

        area_ratio = effective_area_p / src_area_p_total
        patch_biomass_before = max(0.0, B_patch_map[p] * area_ratio)
        loss_total_p = patch_biomass_before

        result = _apply_harvest_loss_to_patch(
            p=p,
            src=src,
            loss_biomass=loss_total_p,
            effective_area_p=effective_area_p,
            cell_area=cell_area,
            params=params,
            C_bar=C_bar,
            Delta=Delta,
            LULC_frac=LULC_frac,
            frac_area=frac_area,
            t_idx=t_idx,
            diag=diag,
        )

        total_effective_area += effective_area_p
        actual_removed_tC += loss_total_p
        total_patch_biomass_before += patch_biomass_before

        total_products += result["products"]
        total_sr_biomass += result["sr_biomass"]
        total_sr_soil += result["sr_soil"]
        total_delta_B_h += result["delta_B_h"]
        total_delta_SS_h += result["delta_SS_h"]
        total_beta_delta_sum += result.get("beta_delta_sum", 0.0)
        total_sigma_delta_sum += result.get("sigma_delta_sum", 0.0)
        total_R += result.get("sr_soil", 0.0)

    # ===== 3. Update top-level LULC fractions =====
    # v -> s harvest changes primary/secondary cover area;
    # s -> s harvest does not change top-level cover area.
    if src == "v":
        j_dst = cover_index["s"]
        LULC_frac["v"] = float(frac_area[j_src, :].sum())
        LULC_frac["s"] = float(frac_area[j_dst, :].sum())

    # ===== 4. Diagnostics =====
    _diag_add(diag, "area_harvest", total_effective_area)
    _diag_add(diag, "harvest_effective_area", total_effective_area)
    _diag_add(
        diag,
        "harvest_effective_area_frac",
        total_effective_area / cell_area if cell_area > 0.0 else 0.0,
    )
    _diag_add(diag, "harvest_extra_area", 0.0)

    # In area-driven mode, this is the biomass implied by the model state on the
    # LUH2 harvest area, not LUH2 *_bioh.
    _diag_add(diag, "harvest_requested_biomass", actual_removed_tC)
    _diag_add(diag, "harvest_met_biomass", actual_removed_tC)
    _diag_add(diag, "harvest_forced_biomass", 0.0)
    _diag_add(diag, "harvest_unmet_biomass", 0.0)
    _diag_add(diag, "harvest_unmet_raw_biomass", 0.0)

    _diag_add(diag, "harvest_biomass_removed", actual_removed_tC)
    _diag_add(diag, "harvest_loss_biomass", actual_removed_tC)
    _diag_add(diag, "harvest_patch_biomass_before", total_patch_biomass_before)

    _diag_add(diag, "harvest_to_products", total_products)
    _diag_add(diag, "harvest_to_SR_from_biomass", total_sr_biomass)
    _diag_add(diag, "harvest_to_SR_from_soil", total_sr_soil)
    _diag_add(diag, "harvest_to_soil", total_sr_biomass + total_sr_soil)

    _diag_add(diag, "harvest_delta_B_h", total_delta_B_h)
    _diag_add(diag, "harvest_delta_SS_h", total_delta_SS_h)
    _diag_add(diag, "harvest_beta_delta_sum", total_beta_delta_sum)
    _diag_add(diag, "harvest_sigma_delta_sum", total_sigma_delta_sum)
    _diag_add(diag, "harvest_R", total_R)

    _diag_add(diag, "harvest_forced_to_products", 0.0)
    _diag_add(diag, "harvest_forced_to_soil", 0.0)

    return actual_removed_tC, 0.0

