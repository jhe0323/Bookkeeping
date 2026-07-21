"""Core LULCC bookkeeping simulator.

This version unifies original-LUH2 and pre-aggregated VSCP inputs, uses calendar
simulation years, precompiles only existing transition channels, supports static
and dynamic PFT maps, and stores carbon state in rolling one-step arrays.
"""
from __future__ import annotations

from collections import defaultdict
import json
import re
import time
from typing import Dict, Iterable, List, Optional, Set, Tuple

import netCDF4 as nc
import numpy as np

from src.parameter_loader import ParameterLoader
from src.carbon_pools_init import (
    LULC_keys,
    cover_index,
    history_index,
    pool_index_Cbar,
    pool_index_Delta,
    initialize_Cbar,
    initialize_Delta,
    initialize_frac_area,
    normalize_pft_fractions,
    refresh_Cbar_from_frac_area,
)
from src.transition import relax_one_year_blue_with_diag
from src.summary_yearly import summarize_year
from src.file_loader import FileLoader, resolve_year_index
from src.events import apply_clearing, apply_abandonment, apply_others
from src.harvest import apply_harvest_luh2


FORMAL_FLUX_KEYS = [
    "Gross_Sources", "Gross_Sinks", "Net_Emissions",
    "Closure_Error", "Atmosphere_Closure_Error",
    "Carbon_Density_Adjustment", "PFT_Remap_Adjustment",
    "Flux_FD", "Flux_NFC", "Flux_FR", "Flux_NFR", "Flux_CAL", "Flux_WHp",
    "Flux_Clearing", "Flux_Abandonment", "Flux_Harvest_Net",
    "Flux_Harvest_SoilSlash", "Flux_Harvest_Regrowth", "Flux_Other",
    "Flux_Products_Total",
]

HARVEST_DIAG_KEYS = [
    "harvest_requested_biomass", "harvest_met_biomass",
    "harvest_unmet_biomass", "harvest_unmet_raw_biomass",
    "harvest_forced_biomass", "harvest_luh2_area_frac",
    "harvest_luh2_area", "harvest_effective_area_frac",
    "harvest_effective_area", "harvest_extra_area",
    "harvest_biomass_removed", "harvest_patch_biomass_before",
    "harvest_loss_biomass", "harvest_to_products", "harvest_to_soil",
    "harvest_to_SR_from_biomass", "harvest_to_SR_from_soil",
    "harvest_forced_to_products", "harvest_forced_to_soil",
    "harvest_beta_delta_sum", "harvest_sigma_delta_sum",
    "harvest_delta_B_h", "harvest_delta_SS_h", "harvest_R",
]

TRANSITION_DIAG_KEYS = [
    "transition_requested_area",
    "transition_applied_area",
    "transition_clipped_area",
]

FORMAL_AREA_DIAG_KEYS = [
    "area_clearing", "area_abandonment", "area_other", "area_harvest",
    "area_deforestation",
] + TRANSITION_DIAG_KEYS + HARVEST_DIAG_KEYS

FORMAL_OUTPUT_KEYS = FORMAL_FLUX_KEYS + FORMAL_AREA_DIAG_KEYS
STATE_OUTPUT_KEYS = [
    "biomass_total", "soil_total", "P1", "P10", "P100",
    "atmosphere", "system_carbon_total",
]

_HARVEST_FAMILIES = ("primf", "primn", "secmf", "secyf", "secnf")
_TRANSITION_RE = re.compile(r"^(.+)_to_(.+?)(?:_frac)?$")


class LULCCSimulator:
    def __init__(
        self,
        pft_path: str,
        config_path: str,
        LULC_path: str,
        trans_path: str,
        *,
        experiment_path: Optional[str],
        input_format: str,
        input_base_year: int,
        start_year: int,
        end_year: int,
        time_encoding: str = "index",
        pft_var: Optional[str] = None,
        pft_base_year: int = 850,
        pft_time_encoding: str = "auto",
        pft_update_mode: str = "annual_conservative",
        pft_min_year_policy: str = "clip",
        pft_max_year_policy: str = "clip",
        lat_slice: Optional[slice] = None,
        lon_slice: Optional[slice] = None,
        area_unit: str = "ha",
        run_metadata: Optional[dict] = None,
        parameter_overrides: Optional[dict] = None,
        compression_level: int = 4,
    ):
        self.params = ParameterLoader(
            config_path,
            experiment_path=experiment_path,
            override_config=parameter_overrides,
        )
        self.input_format = str(input_format).lower()
        if self.input_format not in {"original_luh2", "vscp"}:
            raise ValueError("input_format must be original_luh2 or vscp")

        self.input_base_year = int(input_base_year)
        self.start_year = int(start_year)
        self.end_year = int(end_year)
        self.n_steps = self.end_year - self.start_year
        if self.n_steps <= 0:
            raise ValueError("end_year must be greater than start_year")
        self.time_encoding = str(time_encoding)
        self.pft_update_mode = str(pft_update_mode).lower()
        if self.pft_update_mode not in {"fixed_initial", "annual_conservative"}:
            raise ValueError("pft.update_mode must be fixed_initial or annual_conservative")
        self.area_unit = area_unit
        self.run_metadata = dict(run_metadata or {})
        self.compression_level = int(compression_level)

        self.loader = FileLoader()
        self.state_info = self.loader.inspect_dataset(LULC_path)
        self.trans_info = self.loader.inspect_dataset(trans_path)
        self._validate_spatial_metadata()

        state_start_idx = resolve_year_index(
            self.state_info,
            self.start_year,
            input_base_year=self.input_base_year,
            time_encoding=self.time_encoding,
        )
        transition_start_idx = resolve_year_index(
            self.trans_info,
            self.start_year,
            input_base_year=self.input_base_year,
            time_encoding=self.time_encoding,
        )
        transition_stop_idx = transition_start_idx + self.n_steps
        trans_time_len = (
            self.trans_info.dimensions[self.trans_info.time_name]
            if self.trans_info.time_name else 1
        )
        if transition_stop_idx > trans_time_len:
            raise IndexError(
                f"Transition range {self.start_year}..{self.end_year} requires "
                f"indices {transition_start_idx}:{transition_stop_idx}, but "
                f"{trans_path} has only {trans_time_len} time steps."
            )

        self.state_channels = self._compile_state_channels(self.state_info.variables)
        self.transition_channels = self._compile_transition_channels(self.trans_info.variables)
        self.harvest_channels = self._compile_harvest_channels(self.trans_info.variables)

        state_vars = sorted({name for name, _ in self.state_channels})
        transition_vars = sorted(
            {name for name, _, _ in self.transition_channels}
            | {name for pair in self.harvest_channels.values() for name in pair if name}
        )

        self.state_data, state_meta = self.loader.load_luh2_dataset(
            LULC_path,
            lat_slice=lat_slice,
            lon_slice=lon_slice,
            time_slice=slice(state_start_idx, state_start_idx + 1),
            variable_names=state_vars,
        )
        self.trans_data, trans_meta = self.loader.load_luh2_dataset(
            trans_path,
            lat_slice=lat_slice,
            lon_slice=lon_slice,
            time_slice=slice(transition_start_idx, transition_stop_idx),
            variable_names=transition_vars,
        )

        self.state_lat = state_meta.lat
        self.state_lon = state_meta.lon
        self.trans_lat = trans_meta.lat
        self.trans_lon = trans_meta.lon
        self.state_lat_name = state_meta.lat_name
        self.state_lon_name = state_meta.lon_name
        self.trans_lat_name = trans_meta.lat_name
        self.trans_lon_name = trans_meta.lon_name
        self.actual_ds_lat_is_desc = state_meta.lat_desc
        self.actual_trans_lat_is_desc = trans_meta.lat_desc

        self._build_grid_and_area(
            state_meta.internal_lat,
            state_meta.internal_lon,
            unit=area_unit,
        )
        self.pft_data = self.loader.load_pft_map(
            pft_path,
            pft_var=pft_var,
            target_lat_asc=self._latitudes,
            target_lon_180=self._longitudes,
            target_lon_res=self.lon_res,
            expected_n_pft=self.params.n_pft,
            model_years=range(self.start_year, self.end_year + 1),
            pft_base_year=int(pft_base_year),
            time_encoding=pft_time_encoding,
            min_year_policy=pft_min_year_policy,
            max_year_policy=pft_max_year_policy,
        )
        self.run_metadata.update({
            "pft_variable": self.pft_data.variable_name,
            "pft_source_mode": self.pft_data.source_mode,
            "pft_dynamic": int(self.pft_data.dynamic),
            "transition_channel_count": len(self.transition_channels),
            "state_channel_count": len(self.state_channels),
        })

    def _validate_spatial_metadata(self) -> None:
        state_lat = np.asarray(self.state_info.lat, dtype=float)
        trans_lat = np.asarray(self.trans_info.lat, dtype=float)
        state_lon = np.sort(((np.asarray(self.state_info.lon) + 180) % 360) - 180)
        trans_lon = np.sort(((np.asarray(self.trans_info.lon) + 180) % 360) - 180)
        if state_lat.size != trans_lat.size or not np.allclose(
            np.sort(state_lat), np.sort(trans_lat), atol=1e-7, rtol=0
        ):
            raise ValueError("States and transitions latitude grids differ.")
        if state_lon.size != trans_lon.size or not np.allclose(
            state_lon, trans_lon, atol=1e-7, rtol=0
        ):
            raise ValueError("States and transitions longitude grids differ.")

    def _compile_state_channels(self, variable_names: Iterable[str]) -> List[Tuple[str, str]]:
        available = set(variable_names)
        channels: List[Tuple[str, str]] = []
        if self.input_format == "vscp":
            missing = [name for name in ("v", "s", "p", "c") if name not in available]
            if missing:
                raise KeyError(f"VSCP states file is missing variables: {missing}")
            for cover in LULC_keys:
                if cover in available:
                    channels.append((cover, cover))
            return channels

        mapping = self.params.config.get("LUH2toLULC", {})
        required = list(mapping.keys())
        missing = [name for name in required if name not in available]
        if missing:
            raise KeyError(f"Original LUH2 states file is missing variables: {missing}")
        for variable_name, cover in mapping.items():
            if cover in cover_index:
                channels.append((str(variable_name), str(cover)))
        return channels

    def _compile_transition_channels(
        self,
        variable_names: Iterable[str],
    ) -> List[Tuple[str, str, str]]:
        mapping = self.params.config.get("LUH2toLULC", {})
        channels: List[Tuple[str, str, str]] = []
        seen: Set[Tuple[str, str, str]] = set()
        for name in variable_names:
            match = _TRANSITION_RE.match(str(name))
            if not match:
                continue
            raw_src, raw_dst = match.groups()
            if raw_src.endswith("_harv") or raw_src.endswith("_bioh"):
                continue
            if self.input_format == "vscp":
                src, dst = raw_src, raw_dst
                if src not in cover_index or dst not in cover_index:
                    continue
            else:
                if raw_src not in mapping or raw_dst not in mapping:
                    continue
                src, dst = mapping[raw_src], mapping[raw_dst]
            if src == dst:
                continue
            item = (str(name), str(src), str(dst))
            if item not in seen:
                channels.append(item)
                seen.add(item)
        if not channels:
            raise ValueError("No usable transition channels were found in the transition file.")
        return channels

    def _compile_harvest_channels(self, variable_names: Iterable[str]) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
        available = set(variable_names)
        channels = {}
        for family in _HARVEST_FAMILIES:
            area_name = f"{family}_harv" if f"{family}_harv" in available else None
            bio_name = f"{family}_bioh" if f"{family}_bioh" in available else None
            if self.params.harvest_mode == "area" and area_name is None:
                continue
            if self.params.harvest_use_bioh and bio_name is None and area_name is None:
                continue
            if area_name is not None or bio_name is not None:
                channels[family] = (area_name, bio_name)
        if not channels:
            raise ValueError("No usable LUH2 harvest channels were found.")
        return channels

    def _build_grid_and_area(self, internal_lat, internal_lon, *, unit: str) -> None:
        self._latitudes = np.asarray(internal_lat, dtype=float)
        self._longitudes = np.asarray(internal_lon, dtype=float)
        self.lat_res = (
            abs(float(np.median(np.diff(self._latitudes))))
            if self._latitudes.size >= 2
            else abs(float(np.median(np.diff(np.asarray(self.state_info.lat, dtype=float)))))
        )
        self.lon_res = (
            abs(float(np.median(np.diff(self._longitudes))))
            if self._longitudes.size >= 2
            else abs(float(np.median(np.diff(np.asarray(self.state_info.lon, dtype=float)))))
        )
        radius = 6_371_000.0
        dphi = np.deg2rad(self.lat_res)
        dlambda = np.deg2rad(self.lon_res)
        lat_grid = np.repeat(self._latitudes[:, None], self._longitudes.size, axis=1)
        area_m2 = radius**2 * dphi * dlambda * np.cos(np.deg2rad(lat_grid))
        normalized = unit.lower()
        if normalized in ("m2", "m^2"):
            area = area_m2
        elif normalized == "ha":
            area = area_m2 / 1e4
        elif normalized in ("km2", "km^2"):
            area = area_m2 / 1e6
        else:
            raise ValueError(f"Unsupported area unit: {unit}")
        self._area_grid = area.astype(np.float64)

    def _normalize_ij(self, lat_idx: int, lon_idx: int) -> Tuple[int, int]:
        nlat, nlon = self._area_grid.shape
        i, j = int(lat_idx), int(lon_idx)
        if not (0 <= i < nlat and 0 <= j < nlon):
            raise IndexError(f"Grid index outside local slice: i={i}, j={j}, shape={(nlat, nlon)}")
        return i, j

    def _ds_lat_idx(self, i: int, which: str) -> int:
        if which == "state":
            n = self.state_data.dimensions[self.state_lat_name].size
            descending = self.actual_ds_lat_is_desc
        else:
            n = self.trans_data.dimensions[self.trans_lat_name].size
            descending = self.actual_trans_lat_is_desc
        return n - 1 - int(i) if descending else int(i)

    @staticmethod
    def _clean_fraction(value, *, name: str, year: int, i: int, j: int) -> float:
        value = value.filled(np.nan) if hasattr(value, "filled") else value
        result = float(value)
        if not np.isfinite(result):
            return 0.0
        if result < -1e-10 or result > 1.0 + 1e-8:
            raise ValueError(
                f"Fraction outside [0,1]: variable={name}, year={year}, "
                f"cell=({i},{j}), value={result}"
            )
        return float(np.clip(result, 0.0, 1.0))

    def _get_cell_pft_grid(self, lat_idx: int, lon_idx: int, calendar_year: int) -> np.ndarray:
        i, j = self._normalize_ij(lat_idx, lon_idx)
        values = self.pft_data.get_cell(calendar_year, i, j)
        return normalize_pft_fractions(values, self.params.n_pft)

    def _get_initial_fractions(self, lat_idx: int, lon_idx: int) -> Dict[str, float]:
        i, j = self._normalize_ij(lat_idx, lon_idx)
        ii = self._ds_lat_idx(i, "state")
        fractions = {cover: 0.0 for cover in LULC_keys}
        for variable_name, cover in self.state_channels:
            variable = self.state_data.variables[variable_name]
            value = variable[0, ii, j]
            fractions[cover] += self._clean_fraction(
                value,
                name=variable_name,
                year=self.start_year,
                i=i,
                j=j,
            )
        return fractions

    def _parse_transitions(self, local_year_idx: int, lat_idx: int, lon_idx: int) -> Dict[str, float]:
        i, j = self._normalize_ij(lat_idx, lon_idx)
        ii = self._ds_lat_idx(i, "trans")
        calendar_year = self.start_year + int(local_year_idx)
        transitions: Dict[str, float] = defaultdict(float)
        for variable_name, src, dst in self.transition_channels:
            variable = self.trans_data.variables[variable_name]
            value = variable[local_year_idx, ii, j]
            fraction = self._clean_fraction(
                value,
                name=variable_name,
                year=calendar_year,
                i=i,
                j=j,
            )
            if fraction > 0.0:
                transitions[f"{src}_to_{dst}"] += fraction
        return dict(transitions)

    def _parse_wood_harvest(self, local_year_idx: int, lat_idx: int, lon_idx: int):
        i, j = self._normalize_ij(lat_idx, lon_idx)
        ii = self._ds_lat_idx(i, "trans")
        calendar_year = self.start_year + int(local_year_idx)
        result = {}
        for family, (area_name, bio_name) in self.harvest_channels.items():
            area_frac = 0.0
            biomass_kgC = 0.0
            if area_name is not None:
                area_frac = self._clean_fraction(
                    self.trans_data.variables[area_name][local_year_idx, ii, j],
                    name=area_name,
                    year=calendar_year,
                    i=i,
                    j=j,
                )
            if self.params.harvest_use_bioh and bio_name is not None:
                raw = self.trans_data.variables[bio_name][local_year_idx, ii, j]
                raw = raw.filled(np.nan) if hasattr(raw, "filled") else raw
                biomass_kgC = float(raw)
                if not np.isfinite(biomass_kgC):
                    biomass_kgC = 0.0
                if biomass_kgC < -1e-8:
                    raise ValueError(
                        f"Negative harvest biomass: {bio_name}, year={calendar_year}, "
                        f"cell=({i},{j}), value={biomass_kgC}"
                    )
                biomass_kgC = max(0.0, biomass_kgC)
            active = area_frac > 0.0 or (
                self.params.harvest_use_bioh and biomass_kgC > 0.0
            )
            if active:
                result[family] = {
                    "area_frac": area_frac,
                    "biomass_kgC": biomass_kgC,
                }
        return result

    @staticmethod
    def _system_total(C_bar, Delta, atmosphere: float) -> float:
        return float(C_bar.sum() + Delta.sum() + float(atmosphere))

    def _refresh_dynamic_density(self, C_bar, frac_area, cell_area) -> float:
        if not self.params.use_dynamic_density:
            return 0.0
        before = float(C_bar.sum())
        refresh_Cbar_from_frac_area(C_bar, 0, frac_area, cell_area, self.params)
        return float(C_bar.sum() - before)

    def _conservative_pft_remap(
        self,
        *,
        new_pft: np.ndarray,
        LULC_frac: dict,
        frac_area: np.ndarray,
        C_bar: np.ndarray,
        Delta: np.ndarray,
        cell_area: float,
    ) -> float:
        """Remap cover-level stocks to new PFT fractions without changing total carbon.

        Existing product and rapid-soil pools retain their historical PFT labels.
        Biomass and slow-soil actual stocks are redistributed within each cover;
        equilibrium mismatch is placed in the general-history (g) Delta pool.
        """
        before_total = float(C_bar.sum() + Delta.sum())
        weights = normalize_pft_fractions(new_pft, self.params.n_pft)
        k_general = history_index["g"]

        for land, j in cover_index.items():
            if land == "U":
                continue
            land_fraction = float(LULC_frac.get(land, 0.0))
            if land_fraction <= 0.0:
                frac_area[j, :] = 0.0
                C_bar[0, :, j, :] = 0.0
                continue

            old_actual = {}
            old_history_totals = {}
            for pool_name in ("B", "SS"):
                c_index = pool_index_Cbar[pool_name]
                d_index = pool_index_Delta[pool_name]
                old_actual[pool_name] = float(
                    C_bar[0, c_index, j, :].sum()
                    + Delta[0, d_index, j, :, :].sum()
                )
                old_history_totals[pool_name] = np.asarray(
                    Delta[0, d_index, j, :, :].sum(axis=1),
                    dtype=np.float64,
                )

            frac_area[j, :] = land_fraction * weights
            for p in range(self.params.n_pft):
                area = float(frac_area[j, p]) * cell_area
                C_bar[0, pool_index_Cbar["B"], j, p] = (
                    area * self.params.get_carbon_density(p, "Biomass", land)
                )
                C_bar[0, pool_index_Cbar["SS"], j, p] = (
                    area * self.params.get_carbon_density(p, "Soil", land)
                )

            for pool_name in ("B", "SS"):
                c_index = pool_index_Cbar[pool_name]
                d_index = pool_index_Delta[pool_name]
                Delta[0, d_index, j, :, :] = (
                    old_history_totals[pool_name][:, None] * weights[None, :]
                )
                current = float(
                    C_bar[0, c_index, j, :].sum()
                    + Delta[0, d_index, j, :, :].sum()
                )
                mismatch = old_actual[pool_name] - current
                Delta[0, d_index, j, k_general, :] += mismatch * weights

        after_total = float(C_bar.sum() + Delta.sum())
        return after_total - before_total

    @staticmethod
    def _scale_transitions_to_year_start(
        transitions: Dict[str, float],
        year_start_fraction: Dict[str, float],
        cell_area: float,
    ) -> Tuple[Dict[str, float], float, float, float]:
        by_source: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
        requested_total = 0.0
        for key, fraction in transitions.items():
            src, _ = key.split("_to_", 1)
            by_source[src].append((key, float(fraction)))
            requested_total += float(fraction) * cell_area

        applied: Dict[str, float] = {}
        for src, items in by_source.items():
            requested_fraction = sum(value for _, value in items)
            available_fraction = max(0.0, float(year_start_fraction.get(src, 0.0)))
            scale = 1.0 if requested_fraction <= available_fraction + 1e-15 else (
                available_fraction / requested_fraction if requested_fraction > 0 else 0.0
            )
            for key, value in items:
                scaled = value * scale
                if scaled > 0.0:
                    applied[key] = scaled

        applied_total = sum(applied.values()) * cell_area
        clipped_total = max(0.0, requested_total - applied_total)
        return applied, requested_total, applied_total, clipped_total

    def _apply_transition_batch(
        self,
        *,
        transitions: Dict[str, float],
        cell_area: float,
        pft_grid: np.ndarray,
        C_bar: np.ndarray,
        Delta: np.ndarray,
        LULC_frac: dict,
        frac_area: np.ndarray,
        diag: dict,
    ) -> None:
        """Apply all ordinary transitions against the same beginning-of-year state.

        Each event is evaluated on a private copy of the year-start state.  Its
        carbon, area, and diagnostic increments are accumulated and then applied
        together.  Reciprocal or multi-destination transitions therefore cannot
        consume area/carbon that arrived earlier in the same annual loop.
        """
        if not transitions:
            return

        base_cbar = C_bar.copy()
        base_delta = Delta.copy()
        base_lulc = dict(LULC_frac)
        base_frac = frac_area.copy()

        accumulated_cbar = np.zeros_like(C_bar)
        accumulated_delta = np.zeros_like(Delta)
        accumulated_frac = np.zeros_like(frac_area)
        accumulated_lulc = {key: 0.0 for key in LULC_frac}
        accumulated_diag = {key: 0.0 for key in diag}

        for trans_key, trans_frac in transitions.items():
            src, dst = trans_key.split("_to_", 1)
            area_abs = float(trans_frac) * float(cell_area)
            if area_abs <= 0.0:
                continue

            local_cbar = base_cbar.copy()
            local_delta = base_delta.copy()
            local_lulc = dict(base_lulc)
            local_frac = base_frac.copy()
            local_diag = {key: 0.0 for key in diag}
            common = dict(
                area=area_abs,
                src=src,
                dst=dst,
                cell_area=cell_area,
                pft_grid=pft_grid,
                params=self.params,
                C_bar=local_cbar,
                Delta=local_delta,
                LULC_frac=local_lulc,
                frac_area=local_frac,
                t_idx=0,
                diag=local_diag,
            )

            handled = True
            if src in ("v", "s") and dst in ("p", "c"):
                apply_clearing(**common)
            elif src in ("c", "p") and dst == "s":
                apply_abandonment(**common)
            elif src in ("c", "p") and dst in ("c", "p") and src != dst:
                apply_others(**common)
            else:
                handled = False

            if not handled:
                continue
            accumulated_cbar += local_cbar - base_cbar
            accumulated_delta += local_delta - base_delta
            accumulated_frac += local_frac - base_frac
            for key in accumulated_lulc:
                accumulated_lulc[key] += float(local_lulc[key] - base_lulc[key])
            for key, value in local_diag.items():
                accumulated_diag[key] += float(value)

        C_bar[:] = base_cbar + accumulated_cbar
        Delta[:] = base_delta + accumulated_delta
        frac_area[:] = base_frac + accumulated_frac
        for key in LULC_frac:
            LULC_frac[key] = max(0.0, float(base_lulc[key] + accumulated_lulc[key]))
        for key, value in accumulated_diag.items():
            diag[key] += value

    @staticmethod
    def _empty_record() -> dict:
        record = {key: 0.0 for key in STATE_OUTPUT_KEYS + FORMAL_OUTPUT_KEYS}
        return record

    def run_simulation(self, lat_idx: int, lon_idx: int):
        cell_area = float(self._area_grid[lat_idx, lon_idx])
        LULC_frac = self._get_initial_fractions(lat_idx, lon_idx)
        if sum(LULC_frac.values()) <= 0.0:
            return [self._empty_record() for _ in range(self.n_steps + 1)]

        pft_grid = self._get_cell_pft_grid(lat_idx, lon_idx, self.start_year)
        self.params.set_density_year(self.start_year)
        C_bar = initialize_Cbar(
            LULC_frac,
            cell_area,
            pft_grid,
            self.params,
            n_time=1,
        )
        Delta = initialize_Delta(n_pft=self.params.n_pft, n_time=1)
        frac_area = initialize_frac_area(
            LULC_frac,
            pft_grid,
            n_pft=self.params.n_pft,
        )
        atmosphere = 0.0

        yearly = []
        initial = summarize_year(C_bar, Delta, atmosphere, t_idx=0)
        initial.update({key: 0.0 for key in FORMAL_OUTPUT_KEYS})
        yearly.append(initial)
        total_before = self._system_total(C_bar, Delta, atmosphere)

        for local_year_idx in range(self.n_steps):
            current_year = self.start_year + local_year_idx
            next_year = current_year + 1
            atmosphere_before = atmosphere
            diag = {key: 0.0 for key in FORMAL_OUTPUT_KEYS}

            self.params.set_density_year(next_year)
            density_adjustment = self._refresh_dynamic_density(
                C_bar,
                frac_area,
                cell_area,
            )
            diag["Carbon_Density_Adjustment"] = density_adjustment

            next_pft = self._get_cell_pft_grid(lat_idx, lon_idx, next_year)
            if self.pft_data.dynamic and self.pft_update_mode == "annual_conservative":
                pft_adjustment = self._conservative_pft_remap(
                    new_pft=next_pft,
                    LULC_frac=LULC_frac,
                    frac_area=frac_area,
                    C_bar=C_bar,
                    Delta=Delta,
                    cell_area=cell_area,
                )
                pft_grid = next_pft
            else:
                pft_adjustment = 0.0
            diag["PFT_Remap_Adjustment"] = pft_adjustment

            raw_transitions = self._parse_transitions(
                local_year_idx,
                lat_idx,
                lon_idx,
            )
            transitions, requested, applied, clipped = self._scale_transitions_to_year_start(
                raw_transitions,
                dict(LULC_frac),
                cell_area,
            )
            diag["transition_requested_area"] = requested
            diag["transition_applied_area"] = applied
            diag["transition_clipped_area"] = clipped

            self._apply_transition_batch(
                transitions=transitions,
                cell_area=cell_area,
                pft_grid=pft_grid,
                C_bar=C_bar,
                Delta=Delta,
                LULC_frac=LULC_frac,
                frac_area=frac_area,
                diag=diag,
            )

            wood_harvest = self._parse_wood_harvest(
                local_year_idx,
                lat_idx,
                lon_idx,
            )
            for harvest_type, harvest_input in wood_harvest.items():
                biomass_kgC = float(harvest_input.get("biomass_kgC", 0.0))
                apply_harvest_luh2(
                    mode=self.params.harvest_mode,
                    allow_expansion=self.params.harvest_allow_expansion,
                    harvest_type=harvest_type,
                    area_frac=float(harvest_input.get("area_frac", 0.0)),
                    biomass_harv_kgC=(
                        biomass_kgC
                        if self.params.harvest_use_bioh and biomass_kgC > 0.0
                        else None
                    ),
                    cell_area=cell_area,
                    pft_grid=pft_grid,
                    params=self.params,
                    C_bar=C_bar,
                    Delta=Delta,
                    LULC_frac=LULC_frac,
                    frac_area=frac_area,
                    t_idx=0,
                    diag=diag,
                )

            atmosphere, diag_relax = relax_one_year_blue_with_diag(
                Delta=Delta,
                Atmos=atmosphere,
                params=self.params,
                t_idx=0,
                dt=1.0,
            )
            for key, value in diag_relax.items():
                diag[key] += value

            total_after = self._system_total(C_bar, Delta, atmosphere)
            external_adjustment = density_adjustment + pft_adjustment
            diag["Closure_Error"] = total_after - total_before - external_adjustment
            diag["Atmosphere_Closure_Error"] = (
                diag["Net_Emissions"] - (atmosphere - atmosphere_before)
            )

            record = summarize_year(C_bar, Delta, atmosphere, t_idx=0)
            record.update(diag)
            yearly.append(record)
            total_before = total_after
        return yearly

    def _set_global_attributes(self, ds: nc.Dataset) -> None:
        metadata = dict(self.run_metadata)
        metadata.update({
            "experiment_name": self.params.experiment_name,
            "harvest_mode": self.params.harvest_mode,
            "harvest_use_bioh": int(self.params.harvest_use_bioh),
            "harvest_allow_expansion": int(self.params.harvest_allow_expansion),
            "dynamic_carbon_density": int(self.params.use_dynamic_density),
            "pft_count": int(self.params.n_pft),
            "pft_names": ",".join(self.params.pft_names),
            "area_unit": self.area_unit,
            "carbon_unit": "t C",
            "start_year": self.start_year,
            "end_year": self.end_year,
        })
        for key, value in metadata.items():
            if value is None:
                continue
            if isinstance(value, (dict, list, tuple)):
                value = json.dumps(value, ensure_ascii=False)
            try:
                ds.setncattr(str(key), value)
            except Exception:
                ds.setncattr(str(key), str(value))

    def run_simulation_grid(
        self,
        *,
        out_nc: str,
        lat_slice: Optional[slice] = None,
        lon_slice: Optional[slice] = None,
        sync_every: int = 200,
    ) -> None:
        started = time.time()
        nlat_all, nlon_all = self._area_grid.shape
        lat_slice = lat_slice or slice(0, nlat_all)
        lon_slice = lon_slice or slice(0, nlon_all)
        i0, i1 = lat_slice.start or 0, lat_slice.stop or nlat_all
        j0, j1 = lon_slice.start or 0, lon_slice.stop or nlon_all
        ni, nj = i1 - i0, j1 - j0
        if ni <= 0 or nj <= 0:
            raise ValueError(f"Empty output slice: lat={lat_slice}, lon={lon_slice}")

        time_values = np.arange(self.start_year, self.end_year + 1, dtype=np.int32)
        nt = time_values.size
        chunks = (min(nt, 64), 1, 1)
        kwargs = {
            "zlib": True,
            "complevel": self.compression_level,
            "fill_value": np.nan,
            "chunksizes": chunks,
        }

        ds = nc.Dataset(out_nc, "w")
        try:
            self._set_global_attributes(ds)
            ds.createDimension("time", nt)
            ds.createDimension("lat", ni)
            ds.createDimension("lon", nj)
            vtime = ds.createVariable("time", "i4", ("time",))
            vlat = ds.createVariable("lat", "f4", ("lat",))
            vlon = ds.createVariable("lon", "f4", ("lon",))
            vtime[:] = time_values
            vlat[:] = self._latitudes[i0:i1].astype(np.float32)
            vlon[:] = self._longitudes[j0:j1].astype(np.float32)
            vtime.units = "calendar year"
            vtime.long_name = "calendar year"
            vlat.units = "degrees_north"
            vlon.units = "degrees_east"

            output_names = tuple(STATE_OUTPUT_KEYS + FORMAL_OUTPUT_KEYS)
            out_vars = {
                name: ds.createVariable(name, "f8", ("time", "lat", "lon"), **kwargs)
                for name in output_names
            }
            for variable in out_vars.values():
                variable.units = "t C"
            area_variables = {
                "area_clearing", "area_abandonment", "area_other", "area_harvest",
                "area_deforestation", "transition_requested_area",
                "transition_applied_area", "transition_clipped_area",
                "harvest_luh2_area", "harvest_effective_area", "harvest_extra_area",
            }
            for name in area_variables:
                out_vars[name].units = self.area_unit
            out_vars["harvest_luh2_area_frac"].units = "1"
            out_vars["harvest_effective_area_frac"].units = "1"
            out_vars["Closure_Error"].long_name = (
                "Change in total system carbon minus external carbon-density/PFT adjustments"
            )
            out_vars["Atmosphere_Closure_Error"].long_name = (
                "Net emissions minus annual atmosphere-pool increment"
            )
            out_vars["Carbon_Density_Adjustment"].long_name = (
                "External system-carbon adjustment caused by transient density refresh"
            )
            out_vars["PFT_Remap_Adjustment"].long_name = (
                "Numerical carbon adjustment during conservative PFT remapping; expected near zero"
            )

            done_var = ds.createVariable(
                "done",
                "i1",
                ("lat", "lon"),
                zlib=True,
                complevel=1,
                chunksizes=(min(ni, 64), min(nj, 64)),
            )
            done_var[:, :] = 0

            total = ni * nj
            completed = 0
            ema = None
            for ii, i in enumerate(range(i0, i1)):
                for jj, j in enumerate(range(j0, j1)):
                    t_cell = time.time()
                    initial_fractions = self._get_initial_fractions(i, j)
                    if sum(initial_fractions.values()) <= 0.0:
                        for variable in out_vars.values():
                            variable[:, ii, jj] = 0.0
                        done_var[ii, jj] = 1
                        completed += 1
                        continue

                    yearly = self.run_simulation(i, j)
                    block = np.empty((nt, len(output_names)), dtype=np.float64)
                    for t, record in enumerate(yearly):
                        block[t, :] = [float(record.get(name, 0.0)) for name in output_names]
                    for k, name in enumerate(output_names):
                        out_vars[name][:, ii, jj] = block[:, k]
                    done_var[ii, jj] = 1

                    elapsed = time.time() - t_cell
                    ema = elapsed if ema is None else 0.25 * elapsed + 0.75 * ema
                    completed += 1
                    percent = 100.0 * completed / total
                    eta = (total - completed) * (ema or 0.0)
                    print(
                        f"\r[{completed:6d}/{total:6d}] {percent:5.1f}% "
                        f"cell=({i},{j}) ETA={eta/3600:.2f}h",
                        end="",
                        flush=True,
                    )
                    if sync_every and completed % sync_every == 0:
                        ds.sync()
            print("\n[grid] DONE")
            ds.history = (
                f"Created {time.strftime('%Y-%m-%d %H:%M:%S')}; "
                f"rolling carbon state; elapsed_seconds={time.time()-started:.1f}"
            )
        finally:
            ds.close()
