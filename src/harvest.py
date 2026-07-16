"""Wood-harvest implementations shared by all experiment configurations."""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from src.parameter_loader import ParameterLoader
from src.carbon_pools_init import (
    cover_index,
    history_index,
    pool_index_Cbar,
    pool_index_Delta,
)


HARVEST_TYPE_MAP = {
    "primf": {"src_cover": "v", "veg_group": "forest"},
    "primn": {"src_cover": "v", "veg_group": "nonforest"},
    "secmf": {"src_cover": "s", "veg_group": "forest"},
    "secyf": {"src_cover": "s", "veg_group": "forest"},
    "secnf": {"src_cover": "s", "veg_group": "nonforest"},
}


def _diag_add(diag: Optional[dict], key: str, value: float) -> None:
    """Add a value only when the caller requested that diagnostic key."""
    if diag is not None and key in diag:
        diag[key] += float(value)


def _is_forest_pft(p: int, params: ParameterLoader) -> bool:
    return p in params.get_forest_pfts()


def _pft_allowed_for_harvest_type(
    p: int,
    harvest_type: str,
    params: ParameterLoader,
) -> bool:
    if harvest_type not in HARVEST_TYPE_MAP:
        raise ValueError(f"Unknown harvest_type: {harvest_type}")

    veg_group = HARVEST_TYPE_MAP[harvest_type]["veg_group"]
    is_forest = _is_forest_pft(p, params)

    if veg_group == "forest":
        return is_forest
    if veg_group == "nonforest":
        return not is_forest
    raise ValueError(f"Unsupported harvest vegetation group: {veg_group}")


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


def simulate_harvest(
    loss_biomass: float,
    *,
    src: str,
    p: int,
    params: ParameterLoader,
) -> Dict[str, float]:
    """Partition physically removed biomass into product and rapid-soil pools."""
    if loss_biomass <= 0.0:
        return {"P1": 0.0, "P10": 0.0, "P100": 0.0, "SR": 0.0}

    h_soil = float(params.get_harvest_param(p, f"Bio_Soil_{src}"))
    h_p1 = float(params.get_harvest_param(p, "Prod_1"))
    h_p10 = float(params.get_harvest_param(p, "Prod_10"))
    h_p100 = float(params.get_harvest_param(p, "Prod_100"))

    to_soil = loss_biomass * h_soil
    to_products = loss_biomass * (1.0 - h_soil)

    p1 = to_products * h_p1
    p10 = to_products * h_p10
    p100 = to_products * h_p100

    product_sum = p1 + p10 + p100
    if product_sum > to_products and product_sum > 0.0:
        factor = to_products / product_sum
        p1 *= factor
        p10 *= factor
        p100 *= factor

    return {
        "P1": float(p1),
        "P10": float(p10),
        "P100": float(p100),
        "SR": float(to_soil),
    }


def _allocate_by_capacity(
    target_tC: float,
    capacity_by_pft: Dict[int, float],
) -> Tuple[Dict[int, float], float]:
    """Allocate a biomass target proportionally to non-negative PFT capacity."""
    capacities = {
        p: max(0.0, float(value))
        for p, value in capacity_by_pft.items()
    }

    if target_tC <= 0.0:
        return {p: 0.0 for p in capacities}, 0.0

    total_capacity = sum(capacities.values())
    if total_capacity <= 0.0:
        return {p: 0.0 for p in capacities}, float(target_tC)

    allocated_total = min(float(target_tC), total_capacity)
    allocation = {
        p: allocated_total * capacity / total_capacity
        for p, capacity in capacities.items()
    }
    unmet = max(0.0, float(target_tC) - allocated_total)
    return allocation, unmet


def _empty_patch_result() -> Dict[str, float]:
    return {
        "products": 0.0,
        "sr_biomass": 0.0,
        "sr_soil": 0.0,
        "delta_B_h": 0.0,
        "delta_SS_h": 0.0,
        "beta_delta_sum": 0.0,
        "sigma_delta_sum": 0.0,
    }


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
) -> Dict[str, float]:
    """Apply a real biomass removal and associated soil disturbance to one PFT."""
    if effective_area_p <= 0.0 or cell_area <= 0.0:
        return _empty_patch_result()

    loss_biomass = max(0.0, float(loss_biomass))

    j_src = cover_index[src]
    j_dst = cover_index["s"]
    k_har = history_index["h"]

    src_area_p_total = float(frac_area[j_src, p]) * cell_area
    if src_area_p_total <= 0.0:
        return _empty_patch_result()

    effective_area_p = min(float(effective_area_p), src_area_p_total)
    intensity_p = np.clip(effective_area_p / src_area_p_total, 0.0, 1.0)

    # Biomass pools on the affected area.
    b_eq_src = effective_area_p * float(
        params.get_carbon_density(p, "Biomass", src)
    )
    b_eq_dst = effective_area_p * float(
        params.get_carbon_density(p, "Biomass", "s")
    )
    beta_delta_sum = float(
        Delta[t_idx, pool_index_Delta["B"], j_src, :, p].sum()
    ) * intensity_p

    Delta[t_idx, pool_index_Delta["B"], j_src, :, p] *= (1.0 - intensity_p)

    delta_b_h = beta_delta_sum + (b_eq_src - b_eq_dst) - loss_biomass
    Delta[t_idx, pool_index_Delta["B"], j_dst, k_har, p] += delta_b_h

    allocation = simulate_harvest(
        loss_biomass,
        src=src,
        p=p,
        params=params,
    )
    Delta[t_idx, pool_index_Delta["SR"], j_dst, k_har, p] += allocation["SR"]
    Delta[t_idx, pool_index_Delta["P1"], j_dst, k_har, p] += allocation["P1"]
    Delta[t_idx, pool_index_Delta["P10"], j_dst, k_har, p] += allocation["P10"]
    Delta[t_idx, pool_index_Delta["P100"], j_dst, k_har, p] += allocation["P100"]

    products = allocation["P1"] + allocation["P10"] + allocation["P100"]
    sr_biomass = allocation["SR"]

    # Slow-soil disturbance on the affected area.
    ss_eq_src = effective_area_p * float(
        params.get_carbon_density(p, "Soil", src)
    )
    ss_eq_dst = effective_area_p * float(
        params.get_carbon_density(p, "Soil", "s")
    )
    sigma_delta_sum = float(
        Delta[t_idx, pool_index_Delta["SS"], j_src, :, p].sum()
    ) * intensity_p
    sigma_current = ss_eq_src + sigma_delta_sum

    soc_min_density = float(params.get_harvest_param(p, f"SOC_min_{src}"))
    soil_min_patch = effective_area_p * soc_min_density
    soil_to_rapid = max(0.0, sigma_current - soil_min_patch)

    Delta[t_idx, pool_index_Delta["SS"], j_src, :, p] *= (1.0 - intensity_p)

    delta_ss_h = sigma_delta_sum + (ss_eq_src - ss_eq_dst) - soil_to_rapid
    Delta[t_idx, pool_index_Delta["SR"], j_dst, k_har, p] += soil_to_rapid
    Delta[t_idx, pool_index_Delta["SS"], j_dst, k_har, p] += delta_ss_h

    # Primary harvest transfers affected area from primary to secondary land.
    if src == "v":
        C_bar[t_idx, pool_index_Cbar["B"], j_src, p] -= b_eq_src
        C_bar[t_idx, pool_index_Cbar["SS"], j_src, p] -= ss_eq_src
        C_bar[t_idx, pool_index_Cbar["B"], j_dst, p] += b_eq_dst
        C_bar[t_idx, pool_index_Cbar["SS"], j_dst, p] += ss_eq_dst

        area_fraction = effective_area_p / cell_area
        frac_area[j_src, p] -= area_fraction
        frac_area[j_dst, p] += area_fraction
        if frac_area[j_src, p] < 0.0:
            frac_area[j_src, p] = 0.0

    return {
        "products": float(products),
        "sr_biomass": float(sr_biomass),
        "sr_soil": float(soil_to_rapid),
        "delta_B_h": float(delta_b_h),
        "delta_SS_h": float(delta_ss_h),
        "beta_delta_sum": float(beta_delta_sum),
        "sigma_delta_sum": float(sigma_delta_sum),
    }


def _update_primary_secondary_fractions(src: str, LULC_frac, frac_area) -> None:
    if src != "v":
        return
    LULC_frac["v"] = float(frac_area[cover_index["v"], :].sum())
    LULC_frac["s"] = float(frac_area[cover_index["s"], :].sum())


def _valid_area_fraction(value: float) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return 0.0
    if not np.isfinite(value):
        return 0.0
    return float(np.clip(value, 0.0, 1.0))


def _apply_area_driven_harvest(
    *,
    harvest_type: str,
    area_frac: float,
    cell_area: float,
    pft_grid,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    t_idx: int,
    diag: Optional[dict],
) -> Tuple[float, float]:
    """Use only the LUH2 harvest footprint and diagnose removal from model stocks."""
    if cell_area <= 0.0 or harvest_type not in HARVEST_TYPE_MAP:
        return 0.0, 0.0

    area_frac = _valid_area_fraction(area_frac)
    luh2_area_reported = area_frac * float(cell_area)
    _diag_add(diag, "harvest_luh2_area_frac", area_frac)
    _diag_add(diag, "harvest_luh2_area", luh2_area_reported)

    if luh2_area_reported <= 0.0:
        return 0.0, 0.0

    src = HARVEST_TYPE_MAP[harvest_type]["src_cover"]
    j_src = cover_index[src]
    if float(LULC_frac.get(src, 0.0)) <= 0.0:
        return 0.0, 0.0

    allowed_pfts = [
        p
        for p in range(len(pft_grid))
        if _pft_allowed_for_harvest_type(p, harvest_type, params)
    ]

    source_area_by_pft: Dict[int, float] = {}
    biomass_by_pft: Dict[int, float] = {}
    allowed_source_area = 0.0

    for p in allowed_pfts:
        source_area = float(frac_area[j_src, p]) * cell_area
        if source_area <= 0.0:
            continue
        source_area_by_pft[p] = source_area
        biomass_by_pft[p] = _compute_patch_current_biomass(
            t_idx=t_idx,
            j_src=j_src,
            p=p,
            C_bar=C_bar,
            Delta=Delta,
        )
        allowed_source_area += source_area

    if allowed_source_area <= 0.0:
        return 0.0, 0.0

    target_area = min(luh2_area_reported, allowed_source_area)

    total_effective_area = 0.0
    total_products = 0.0
    total_sr_biomass = 0.0
    total_sr_soil = 0.0
    total_delta_b_h = 0.0
    total_delta_ss_h = 0.0
    total_beta_delta = 0.0
    total_sigma_delta = 0.0
    total_removed = 0.0
    total_patch_biomass_before = 0.0

    for p, source_area in source_area_by_pft.items():
        effective_area = target_area * source_area / allowed_source_area
        effective_area = min(effective_area, source_area)
        if effective_area <= 0.0:
            continue

        area_ratio = effective_area / source_area
        patch_biomass = max(0.0, biomass_by_pft[p] * area_ratio)

        result = _apply_harvest_loss_to_patch(
            p=p,
            src=src,
            loss_biomass=patch_biomass,
            effective_area_p=effective_area,
            cell_area=cell_area,
            params=params,
            C_bar=C_bar,
            Delta=Delta,
            LULC_frac=LULC_frac,
            frac_area=frac_area,
            t_idx=t_idx,
        )

        total_effective_area += effective_area
        total_removed += patch_biomass
        total_patch_biomass_before += patch_biomass
        total_products += result["products"]
        total_sr_biomass += result["sr_biomass"]
        total_sr_soil += result["sr_soil"]
        total_delta_b_h += result["delta_B_h"]
        total_delta_ss_h += result["delta_SS_h"]
        total_beta_delta += result["beta_delta_sum"]
        total_sigma_delta += result["sigma_delta_sum"]

    _update_primary_secondary_fractions(src, LULC_frac, frac_area)

    _diag_add(diag, "area_harvest", total_effective_area)
    _diag_add(diag, "harvest_effective_area", total_effective_area)
    _diag_add(
        diag,
        "harvest_effective_area_frac",
        total_effective_area / cell_area,
    )
    _diag_add(diag, "harvest_extra_area", 0.0)

    _diag_add(diag, "harvest_requested_biomass", total_removed)
    _diag_add(diag, "harvest_met_biomass", total_removed)
    _diag_add(diag, "harvest_forced_biomass", 0.0)
    _diag_add(diag, "harvest_unmet_biomass", 0.0)
    _diag_add(diag, "harvest_unmet_raw_biomass", 0.0)

    _diag_add(diag, "harvest_biomass_removed", total_removed)
    _diag_add(diag, "harvest_loss_biomass", total_removed)
    _diag_add(diag, "harvest_patch_biomass_before", total_patch_biomass_before)

    _diag_add(diag, "harvest_to_products", total_products)
    _diag_add(diag, "harvest_to_SR_from_biomass", total_sr_biomass)
    _diag_add(diag, "harvest_to_SR_from_soil", total_sr_soil)
    _diag_add(diag, "harvest_to_soil", total_sr_biomass + total_sr_soil)

    _diag_add(diag, "harvest_delta_B_h", total_delta_b_h)
    _diag_add(diag, "harvest_delta_SS_h", total_delta_ss_h)
    _diag_add(diag, "harvest_beta_delta_sum", total_beta_delta)
    _diag_add(diag, "harvest_sigma_delta_sum", total_sigma_delta)
    _diag_add(diag, "harvest_R", total_sr_soil)

    _diag_add(diag, "harvest_forced_to_products", 0.0)
    _diag_add(diag, "harvest_forced_to_soil", 0.0)
    return total_removed, 0.0


def _record_unmet(diag: Optional[dict], requested_tC: float) -> Tuple[float, float]:
    _diag_add(diag, "harvest_unmet_biomass", requested_tC)
    _diag_add(diag, "harvest_unmet_raw_biomass", requested_tC)
    return 0.0, 0.0


def _apply_biomass_demand_harvest(
    *,
    harvest_type: str,
    area_frac: float,
    biomass_harv_kgC: Optional[float],
    allow_expansion: bool,
    cell_area: float,
    pft_grid,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    t_idx: int,
    diag: Optional[dict],
) -> Tuple[float, float]:
    """Meet LUH2 biomass demand within the footprint, optionally expanding area."""
    area_frac = _valid_area_fraction(area_frac)
    luh2_area_reported = area_frac * max(0.0, float(cell_area))
    _diag_add(diag, "harvest_luh2_area_frac", area_frac)
    _diag_add(diag, "harvest_luh2_area", luh2_area_reported)

    if biomass_harv_kgC is None:
        requested_tC = 0.0
    else:
        try:
            raw_biomass = float(biomass_harv_kgC)
        except (TypeError, ValueError):
            raw_biomass = 0.0
        requested_tC = raw_biomass / 1000.0 if np.isfinite(raw_biomass) else 0.0
        requested_tC = max(0.0, requested_tC)

    _diag_add(diag, "harvest_requested_biomass", requested_tC)
    if requested_tC <= 0.0:
        return 0.0, 0.0

    if cell_area <= 0.0 or harvest_type not in HARVEST_TYPE_MAP:
        return _record_unmet(diag, requested_tC)

    src = HARVEST_TYPE_MAP[harvest_type]["src_cover"]
    j_src = cover_index[src]
    if float(LULC_frac.get(src, 0.0)) <= 0.0:
        return _record_unmet(diag, requested_tC)

    allowed_pfts = [
        p
        for p in range(len(pft_grid))
        if _pft_allowed_for_harvest_type(p, harvest_type, params)
    ]
    if not allowed_pfts:
        return _record_unmet(diag, requested_tC)

    source_area_by_pft: Dict[int, float] = {}
    biomass_by_pft: Dict[int, float] = {}
    allowed_source_area = 0.0

    for p in allowed_pfts:
        source_area = float(frac_area[j_src, p]) * cell_area
        if source_area <= 0.0:
            continue

        patch_biomass = max(
            0.0,
            _compute_patch_current_biomass(
                t_idx=t_idx,
                j_src=j_src,
                p=p,
                C_bar=C_bar,
                Delta=Delta,
            ),
        )
        if patch_biomass <= 0.0:
            continue

        source_area_by_pft[p] = source_area
        biomass_by_pft[p] = patch_biomass
        allowed_source_area += source_area

    if allowed_source_area <= 0.0:
        return _record_unmet(diag, requested_tC)

    footprint_area = min(luh2_area_reported, allowed_source_area)
    footprint_area_by_pft: Dict[int, float] = {}
    footprint_capacity_by_pft: Dict[int, float] = {}
    extra_capacity_by_pft: Dict[int, float] = {}

    for p, source_area in source_area_by_pft.items():
        footprint_area_p = footprint_area * source_area / allowed_source_area
        footprint_area_p = min(footprint_area_p, source_area)

        footprint_capacity = biomass_by_pft[p] * footprint_area_p / source_area
        footprint_area_by_pft[p] = footprint_area_p
        footprint_capacity_by_pft[p] = max(0.0, footprint_capacity)
        extra_capacity_by_pft[p] = max(0.0, biomass_by_pft[p] - footprint_capacity)

    take_footprint, unmet_after_footprint = _allocate_by_capacity(
        requested_tC,
        footprint_capacity_by_pft,
    )

    if allow_expansion:
        take_extra, unmet_final = _allocate_by_capacity(
            unmet_after_footprint,
            extra_capacity_by_pft,
        )
    else:
        take_extra = {p: 0.0 for p in source_area_by_pft}
        unmet_final = unmet_after_footprint

    met_tC = sum(take_footprint.values())
    forced_tC = sum(take_extra.values())
    actual_removed_tC = met_tC + forced_tC

    total_effective_area = 0.0
    total_extra_area = 0.0
    total_products = 0.0
    total_sr_biomass = 0.0
    total_sr_soil = 0.0
    total_delta_b_h = 0.0
    total_delta_ss_h = 0.0
    total_beta_delta = 0.0
    total_sigma_delta = 0.0
    total_patch_biomass_before = 0.0

    for p, source_area in source_area_by_pft.items():
        loss_footprint = float(take_footprint.get(p, 0.0))
        loss_extra = float(take_extra.get(p, 0.0))
        loss_total = loss_footprint + loss_extra
        if loss_total <= 0.0:
            continue

        footprint_area_p = footprint_area_by_pft[p]
        footprint_capacity = footprint_capacity_by_pft[p]
        extra_capacity = extra_capacity_by_pft[p]

        if footprint_capacity > 0.0 and loss_footprint > 0.0:
            effective_footprint_area = (
                footprint_area_p * loss_footprint / footprint_capacity
            )
        else:
            effective_footprint_area = 0.0

        outside_area = max(0.0, source_area - footprint_area_p)
        if extra_capacity > 0.0 and loss_extra > 0.0:
            effective_extra_area = outside_area * loss_extra / extra_capacity
        else:
            effective_extra_area = 0.0

        effective_area = effective_footprint_area + effective_extra_area
        result = _apply_harvest_loss_to_patch(
            p=p,
            src=src,
            loss_biomass=loss_total,
            effective_area_p=effective_area,
            cell_area=cell_area,
            params=params,
            C_bar=C_bar,
            Delta=Delta,
            LULC_frac=LULC_frac,
            frac_area=frac_area,
            t_idx=t_idx,
        )

        total_effective_area += effective_area
        total_extra_area += effective_extra_area
        total_products += result["products"]
        total_sr_biomass += result["sr_biomass"]
        total_sr_soil += result["sr_soil"]
        total_delta_b_h += result["delta_B_h"]
        total_delta_ss_h += result["delta_SS_h"]
        total_beta_delta += result["beta_delta_sum"]
        total_sigma_delta += result["sigma_delta_sum"]
        total_patch_biomass_before += biomass_by_pft[p] * effective_area / source_area

    _update_primary_secondary_fractions(src, LULC_frac, frac_area)

    _diag_add(diag, "area_harvest", total_effective_area)
    _diag_add(diag, "harvest_effective_area", total_effective_area)
    _diag_add(
        diag,
        "harvest_effective_area_frac",
        total_effective_area / cell_area,
    )
    _diag_add(diag, "harvest_extra_area", total_extra_area)

    _diag_add(diag, "harvest_met_biomass", met_tC)
    _diag_add(diag, "harvest_forced_biomass", forced_tC)
    _diag_add(diag, "harvest_unmet_biomass", unmet_final)
    _diag_add(diag, "harvest_unmet_raw_biomass", unmet_final)

    _diag_add(diag, "harvest_biomass_removed", actual_removed_tC)
    _diag_add(diag, "harvest_loss_biomass", actual_removed_tC)
    _diag_add(diag, "harvest_patch_biomass_before", total_patch_biomass_before)

    _diag_add(diag, "harvest_to_products", total_products)
    _diag_add(diag, "harvest_to_SR_from_biomass", total_sr_biomass)
    _diag_add(diag, "harvest_to_SR_from_soil", total_sr_soil)
    _diag_add(diag, "harvest_to_soil", total_sr_biomass + total_sr_soil)

    _diag_add(diag, "harvest_delta_B_h", total_delta_b_h)
    _diag_add(diag, "harvest_delta_SS_h", total_delta_ss_h)
    _diag_add(diag, "harvest_beta_delta_sum", total_beta_delta)
    _diag_add(diag, "harvest_sigma_delta_sum", total_sigma_delta)
    _diag_add(diag, "harvest_R", total_sr_soil)

    # Retained for backward-compatible output schemas. Forced removal is already
    # included in the physical product and soil totals above.
    _diag_add(diag, "harvest_forced_to_products", 0.0)
    _diag_add(diag, "harvest_forced_to_soil", 0.0)
    return met_tC, forced_tC


def apply_harvest_luh2(
    *,
    mode: str,
    allow_expansion: bool,
    harvest_type: str,
    area_frac: float,
    biomass_harv_kgC: Optional[float],
    cell_area: float,
    pft_grid,
    params: ParameterLoader,
    C_bar,
    Delta,
    LULC_frac,
    frac_area,
    t_idx: int = 0,
    diag: Optional[dict] = None,
) -> Tuple[float, float]:
    """Unified entry point for area, strict-biomass, and forced-biomass modes."""
    if mode == "area":
        if allow_expansion:
            raise ValueError("Area-driven harvest cannot enable area expansion.")
        return _apply_area_driven_harvest(
            harvest_type=harvest_type,
            area_frac=area_frac,
            cell_area=cell_area,
            pft_grid=pft_grid,
            params=params,
            C_bar=C_bar,
            Delta=Delta,
            LULC_frac=LULC_frac,
            frac_area=frac_area,
            t_idx=t_idx,
            diag=diag,
        )

    if mode == "biomass_demand":
        return _apply_biomass_demand_harvest(
            harvest_type=harvest_type,
            area_frac=area_frac,
            biomass_harv_kgC=biomass_harv_kgC,
            allow_expansion=bool(allow_expansion),
            cell_area=cell_area,
            pft_grid=pft_grid,
            params=params,
            C_bar=C_bar,
            Delta=Delta,
            LULC_frac=LULC_frac,
            frac_area=frac_area,
            t_idx=t_idx,
            diag=diag,
        )

    raise ValueError(
        f"Unsupported harvest mode {mode!r}; expected 'area' or 'biomass_demand'."
    )
