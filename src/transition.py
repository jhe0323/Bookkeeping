# relaxation_blue.py
from __future__ import annotations

import numpy as np
from typing import Dict, Any

from src.carbon_pools_init import history_index, pool_index_Delta


_EFOLD = 0.534  # keep consistent with transition.response_func / excess_remaining


def _decay_factor(dt: float, tau: float) -> float:
    """Annual decay factor to relax Delta -> 0 using BLUE response-function convention."""
    if tau <= 1e-9:
        return 0.0
    return float(np.exp(-float(dt) / (_EFOLD * float(tau))))

def _relax_slice_inplace(arr: np.ndarray, fac: float) -> float:
    old = arr.copy()
    arr *= fac
    return float((arr - old).sum())

def relax_one_year_blue(
    *,
    Delta: np.ndarray,
    Atmos: float,
    params,
    t_idx: int,
    dt: float = 1.0,
) -> float:
    """
    Apply one annual BLUE-style relaxation step to ndarray-based Delta storage.
                                                           
    Parameters
    ----------
    Delta : ndarray
        Shape (n_time, N_POOL_DELTA, N_COVER, N_HISTORY, N_PFT)
    Atmos : float
        Atmospheric/cumulative emitted carbon.

    Returns
    -------
    float
        Updated atmospheric pool.
    """
    if not (0 <= t_idx < Delta.shape[0]):
        raise IndexError(f"t_idx out of range: {t_idx}, n_time={Delta.shape[0]}")
    # history indices
    k_clr = history_index["l"]
    k_har = history_index["h"]
    k_abd = history_index["a"]
    k_oth = history_index["g"]

    iB = pool_index_Delta["B"]
    iSS = pool_index_Delta["SS"]
    iSR = pool_index_Delta["SR"]
    iP1 = pool_index_Delta["P1"]
    iP10 = pool_index_Delta["P10"]
    iP100 = pool_index_Delta["P100"]
    iA = pool_index_Delta["A"]

    dC_total = 0.0
    n_pft = Delta.shape[-1]

    # Direct atmosphere emissions from events are transferred immediately at end of time step.
    direct_A = float(Delta[t_idx, iA, :, :, :].sum())
    if direct_A != 0.0:
        Atmos += direct_A
        Delta[t_idx, iA, :, :, :] = 0.0
    
    for p in range(n_pft):
        # Clearing
        tau_clr_sr = float(params.get_clearing_param(p, "t_lapse"))
        tau_clr_ss = float(params.get_clearing_param(p, "t_rest"))
        tau_clr_b = 1.0
        # Abandonment
        tau_abd_b = float(params.get_abandonment_param(p, "t_biomass"))
        tau_abd_ss = float(params.get_abandonment_param(p, "t_soil"))
        # Harvest
        tau_har_sr = float(params.get_harvest_param(p, "t_lapse"))
        tau_har_ss = tau_abd_ss
        tau_har_b = tau_abd_b
        # Other
        # For BLUE "other transitions" (c <-> p), soil relaxation uses the
        # rapid clearing soil timescale (~tau_SR), not abandonment soil timescale.
        tau_oth_b = 1.0
        tau_oth_ss = tau_clr_sr

        if tau_clr_sr > 0.0:
            fac = _decay_factor(dt, tau_clr_sr)
            dC_total += _relax_slice_inplace(Delta[t_idx, iSR, :, k_clr, p], fac)

        if tau_har_sr > 0.0:
            fac = _decay_factor(dt, tau_har_sr)
            dC_total += _relax_slice_inplace(Delta[t_idx, iSR, :, k_har, p], fac)

        if tau_clr_ss > 0.0:
            fac = _decay_factor(dt, tau_clr_ss)
            dC_total += _relax_slice_inplace(Delta[t_idx, iSS, :, k_clr, p], fac)

        if tau_har_ss > 0.0:
            fac = _decay_factor(dt, tau_har_ss)
            dC_total += _relax_slice_inplace(Delta[t_idx, iSS, :, k_har, p], fac)

        if tau_abd_ss > 0.0:
            fac = _decay_factor(dt, tau_abd_ss)
            dC_total += _relax_slice_inplace(Delta[t_idx, iSS, :, k_abd, p], fac)

        if tau_oth_ss > 0.0:
            fac = _decay_factor(dt, tau_oth_ss)
            dC_total += _relax_slice_inplace(Delta[t_idx, iSS, :, k_oth, p], fac)

        if tau_clr_b > 0.0:
            fac = _decay_factor(dt, tau_clr_b)
            dC_total += _relax_slice_inplace(Delta[t_idx, iB, :, k_clr, p], fac)

        if tau_har_b > 0.0:
            fac = _decay_factor(dt, tau_har_b)
            dC_total += _relax_slice_inplace(Delta[t_idx, iB, :, k_har, p], fac)

        if tau_abd_b > 0.0:
            fac = _decay_factor(dt, tau_abd_b)
            dC_total += _relax_slice_inplace(Delta[t_idx, iB, :, k_abd, p], fac)

        if tau_oth_b > 0.0:
            fac = _decay_factor(dt, tau_oth_b)
            dC_total += _relax_slice_inplace(Delta[t_idx, iB, :, k_oth, p], fac)

    # Product pools decay globally irrespective of source cover/history labels.
    dC_total += _relax_slice_inplace(Delta[t_idx, iP1, :, :, :], _decay_factor(dt, 1.0))
    dC_total += _relax_slice_inplace(Delta[t_idx, iP10, :, :, :], _decay_factor(dt, 10.0))
    dC_total += _relax_slice_inplace(Delta[t_idx, iP100, :, :, :], _decay_factor(dt, 100.0))

    Atmos = float(Atmos) - dC_total
    return Atmos

def _relax_slice_inplace_with_release(arr: np.ndarray, fac: float) -> float:
    """
    Relax one Delta slice in place and return the released carbon (>0 means
    carbon left the slice and went to atmosphere).
    """
    old = arr.copy()
    arr *= fac
    return float((old - arr).sum())

def _is_forest_pft(p: int, params) -> bool:
    """
    LUCE/BLUE-compatible forest grouping for the current 15-PFT LUCE setup.
    Keep this consistent with events._is_forest_pft().
    Current convention: PFT1--PFT6 are forest PFTs, p is 0-based.
    """
    return p in params.get_forest_pfts()

def _add_flux(diag: Dict[str, float], value: float) -> None:
    """Accumulate one atmosphere flux into gross source/sink and net emissions."""
    value = float(value)
    if value > 0.0:
        diag["Gross_Sources"] += value
    elif value < 0.0:
        diag["Gross_Sinks"] += value
    diag["Net_Emissions"] += value

def relax_one_year_blue_with_diag(
    *,
    Delta: np.ndarray,
    Atmos: float,
    params,
    t_idx: int,
    dt: float = 1.0,
):
    """
    Apply one annual BLUE-style relaxation step and return updated Atmos plus
    publication-level LUCE/BLUE diagnostics.

    Positive flux means carbon released to atmosphere; negative flux means uptake.
    BLUE event terms include product emissions in their matching event category, so
    Flux_Clearing + Flux_Abandonment + Flux_Harvest_Net + Flux_Other = Net_Emissions.
    """
    if not (0 <= t_idx < Delta.shape[0]):
        raise IndexError(f"t_idx out of range: {t_idx}, n_time={Delta.shape[0]}")

    k_clr = history_index["l"]
    k_har = history_index["h"]
    k_abd = history_index["a"]
    k_oth = history_index["g"]

    iB = pool_index_Delta["B"]
    iSS = pool_index_Delta["SS"]
    iSR = pool_index_Delta["SR"]
    iP1 = pool_index_Delta["P1"]
    iP10 = pool_index_Delta["P10"]
    iP100 = pool_index_Delta["P100"]
    iA = pool_index_Delta["A"]

    diag_relax = {
        "Gross_Sources": 0.0,
        "Gross_Sinks": 0.0,
        "Net_Emissions": 0.0,
        "Flux_FD": 0.0,
        "Flux_NFC": 0.0,
        "Flux_FR": 0.0,
        "Flux_NFR": 0.0,
        "Flux_CAL": 0.0,
        "Flux_WHp": 0.0,
        "Flux_Clearing": 0.0,
        "Flux_Abandonment": 0.0,
        "Flux_Harvest_Net": 0.0,
        "Flux_Harvest_SoilSlash": 0.0,
        "Flux_Harvest_Regrowth": 0.0,
        "Flux_Other": 0.0,
        "Flux_Products_Total": 0.0,
    }

    n_pft = Delta.shape[-1]

    def account_clearing(value: float, p: int) -> None:
        nonlocal Atmos
        value = float(value)
        if value == 0.0:
            return
        Atmos += value
        _add_flux(diag_relax, value)
        diag_relax["Flux_Clearing"] += value
        if _is_forest_pft(p, params):
            diag_relax["Flux_FD"] += value
        else:
            diag_relax["Flux_NFC"] += value

    def account_harvest(value: float, p: int, *, product: bool = False) -> None:
        nonlocal Atmos
        value = float(value)
        if value == 0.0:
            return
        Atmos += value
        _add_flux(diag_relax, value)
        diag_relax["Flux_Harvest_Net"] += value
        if product:
            diag_relax["Flux_WHp"] += value
            diag_relax["Flux_Products_Total"] += value
        else:
            if value > 0.0:
                diag_relax["Flux_Harvest_SoilSlash"] += value
            else:
                diag_relax["Flux_Harvest_Regrowth"] += value
                if _is_forest_pft(p, params):
                    diag_relax["Flux_FR"] += value
                else:
                    diag_relax["Flux_NFR"] += value

    def account_abandonment(value: float, p: int) -> None:
        nonlocal Atmos
        value = float(value)
        if value == 0.0:
            return
        Atmos += value
        _add_flux(diag_relax, value)
        diag_relax["Flux_Abandonment"] += value
        if _is_forest_pft(p, params):
            diag_relax["Flux_FR"] += value
        else:
            diag_relax["Flux_NFR"] += value

    def account_other(value: float) -> None:
        nonlocal Atmos
        value = float(value)
        if value == 0.0:
            return
        Atmos += value
        _add_flux(diag_relax, value)
        diag_relax["Flux_Other"] += value
        diag_relax["Flux_CAL"] += value

    # Direct atmosphere emissions stored by event functions in A.
    for p in range(n_pft):
        rel = float(Delta[t_idx, iA, :, k_clr, p].sum())
        if rel:
            account_clearing(rel, p)
        rel = float(Delta[t_idx, iA, :, k_har, p].sum())
        if rel:
            account_harvest(rel, p)
        rel = float(Delta[t_idx, iA, :, k_abd, p].sum())
        if rel:
            account_abandonment(rel, p)
        rel = float(Delta[t_idx, iA, :, k_oth, p].sum())
        if rel:
            account_other(rel)
    Delta[t_idx, iA, :, :, :] = 0.0
    
    for p in range(n_pft):
        tau_clr_sr = float(params.get_clearing_param(p, "t_lapse"))
        tau_clr_ss = float(params.get_clearing_param(p, "t_rest"))
        tau_clr_b = 1.0

        tau_abd_b = float(params.get_abandonment_param(p, "t_biomass"))
        tau_abd_ss = float(params.get_abandonment_param(p, "t_soil"))

        tau_har_sr = float(params.get_harvest_param(p, "t_lapse"))
        tau_har_ss = tau_abd_ss
        tau_har_b = tau_abd_b

        tau_oth_b = 1.0
        tau_oth_ss = tau_clr_sr

        if tau_clr_sr > 0.0:
            account_clearing(_relax_slice_inplace_with_release(Delta[t_idx, iSR, :, k_clr, p], _decay_factor(dt, tau_clr_sr)), p)
        
        if tau_clr_ss > 0.0:
            account_clearing(_relax_slice_inplace_with_release(Delta[t_idx, iSS, :, k_clr, p], _decay_factor(dt, tau_clr_ss)), p)
        
        if tau_clr_b > 0.0:
            account_clearing(_relax_slice_inplace_with_release(Delta[t_idx, iB, :, k_clr, p], _decay_factor(dt, tau_clr_b)), p)

        if tau_har_sr > 0.0:
            account_harvest(_relax_slice_inplace_with_release(Delta[t_idx, iSR, :, k_har, p], _decay_factor(dt, tau_har_sr)), p)

        if tau_har_ss > 0.0:
            account_harvest(_relax_slice_inplace_with_release(Delta[t_idx, iSS, :, k_har, p], _decay_factor(dt, tau_har_ss)), p)

        if tau_har_b > 0.0:
            account_harvest(_relax_slice_inplace_with_release(Delta[t_idx, iB, :, k_har, p], _decay_factor(dt, tau_har_b)), p)

        if tau_abd_ss > 0.0:
            account_abandonment(_relax_slice_inplace_with_release(Delta[t_idx, iSS, :, k_abd, p], _decay_factor(dt, tau_abd_ss)), p)

        if tau_abd_b > 0.0:
            account_abandonment(_relax_slice_inplace_with_release(Delta[t_idx, iB, :, k_abd, p], _decay_factor(dt, tau_abd_b)), p)
        
        if tau_oth_ss > 0.0:
            account_other(_relax_slice_inplace_with_release(Delta[t_idx, iSS, :, k_oth, p], _decay_factor(dt, tau_oth_ss)))
        
        if tau_oth_b > 0.0:
            account_other(_relax_slice_inplace_with_release(Delta[t_idx, iB, :, k_oth, p], _decay_factor(dt, tau_oth_b)))

    # Product pools: clearing products -> FD/NFC and Flux_Clearing; harvest products -> WHp.
    for pool_i, tau in ((iP1, 1.0), (iP10, 10.0), (iP100, 100.0)):
        fac = _decay_factor(dt, tau)
        for p in range(n_pft):
            rel = _relax_slice_inplace_with_release(Delta[t_idx, pool_i, :, k_clr, p], fac)
            if rel:
                diag_relax["Flux_Products_Total"] += rel
                account_clearing(rel, p)
            rel = _relax_slice_inplace_with_release(Delta[t_idx, pool_i, :, k_har, p], fac)
            if rel:
                account_harvest(rel, p, product=True)

    diag_relax["Net_Emissions"] = diag_relax["Gross_Sources"] + diag_relax["Gross_Sinks"]
    return float(Atmos), diag_relax