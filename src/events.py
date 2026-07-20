"""Event carbon allocation models.

This module contains small, pure functions that allocate carbon pools for
different LULCC event types. They were previously implemented as private
methods on LULCCSimulator and are extracted here to keep the simulator core
lean and easier to maintain.

All functions take a `params: ParameterLoader` instance explicitly.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from src.parameter_loader import ParameterLoader
from src.carbon_pools_init import (
    cover_index,
    history_index,
    pool_index_Cbar,
    pool_index_Delta,
)

def _diag_add(diag: Optional[dict], key: str, value: float) -> None:
    """Add to diagnostics only when the key is requested by the caller."""
    if diag is not None and key in diag:
        diag[key] += float(value)


def _is_forest_pft(p: int, params: ParameterLoader) -> bool:
    """Return whether a zero-based PFT index is configured as forest."""
    return p in params.get_forest_pfts()


def _active_source_pfts(pft_grid, frac_area, j_src: int) -> np.ndarray:
    """Return PFT indices that currently occupy the source land-cover class."""
    pft_fractions = np.asarray(pft_grid, dtype=float).reshape(-1)
    area = np.asarray(frac_area)
    if area.ndim != 2 or area.shape[1] != pft_fractions.size:
        raise ValueError(
            f"frac_area shape {area.shape} is incompatible with "
            f"{pft_fractions.size} PFT classes."
        )
    source = np.asarray(area[j_src, :], dtype=float)
    return np.flatnonzero(np.isfinite(source) & (source > 0.0))

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

    for p in _active_source_pfts(pft_grid, frac_area, j_src):
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

    for p in _active_source_pfts(pft_grid, frac_area, j_src):
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

    for p in _active_source_pfts(pft_grid, frac_area, j_src):
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
