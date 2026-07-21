from __future__ import annotations

from typing import Optional, Union
import numpy as np

from src.carbon_pools_init import pool_index_Cbar, pool_index_Delta


def summarize_year(
    C_bar: np.ndarray,
    Delta: np.ndarray,
    Atmos: Union[float, np.ndarray],
    t_idx: int,
    round_ndigits: Optional[int] = None,
):
    iB_c = pool_index_Cbar["B"]
    iSS_c = pool_index_Cbar["SS"]
    iB_d = pool_index_Delta["B"]
    iSS_d = pool_index_Delta["SS"]
    iSR_d = pool_index_Delta["SR"]
    iP1_d = pool_index_Delta["P1"]
    iP10_d = pool_index_Delta["P10"]
    iP100_d = pool_index_Delta["P100"]

    biomass_total = float(C_bar[t_idx, iB_c].sum() + Delta[t_idx, iB_d].sum())
    soil_total = float(
        C_bar[t_idx, iSS_c].sum()
        + Delta[t_idx, iSS_d].sum()
        + Delta[t_idx, iSR_d].sum()
    )
    P1 = float(Delta[t_idx, iP1_d].sum())
    P10 = float(Delta[t_idx, iP10_d].sum())
    P100 = float(Delta[t_idx, iP100_d].sum())
    atmosphere = float(Atmos if np.isscalar(Atmos) else np.asarray(Atmos)[t_idx])
    system_total = biomass_total + soil_total + P1 + P10 + P100 + atmosphere

    if round_ndigits is not None:
        biomass_total = round(biomass_total, round_ndigits)
        soil_total = round(soil_total, round_ndigits)
        P1 = round(P1, round_ndigits)
        P10 = round(P10, round_ndigits)
        P100 = round(P100, round_ndigits)
        atmosphere = round(atmosphere, round_ndigits)
        system_total = round(system_total, round_ndigits)

    return {
        "biomass_total": biomass_total,
        "soil_total": soil_total,
        "P1": P1,
        "P10": P10,
        "P100": P100,
        "atmosphere": atmosphere,
        "system_carbon_total": system_total,
    }
