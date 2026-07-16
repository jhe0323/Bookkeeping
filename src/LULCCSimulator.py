"""Core script. Simulate LULCC and carbon change"""
# simulator.py
# -*- coding: UTF-8 -*-

import copy
from typing import Dict, Optional, Tuple

import netCDF4 as nc
import numpy as np

from src.parameter_loader import ParameterLoader
from src.carbon_pools_init import (
    initialize_Cbar,
    initialize_Delta,
    initialize_frac_area,
    refresh_Cbar_from_frac_area,
)
from src.transition import (
    # relax_one_year_blue,
    relax_one_year_blue_with_diag
)
from src.summary_yearly import summarize_year
from src.file_loader import FileLoader
from src.events import (
    apply_clearing,
    apply_abandonment,
    apply_others,
)
from src.harvest import apply_harvest_luh2

FORMAL_FLUX_KEYS = [
    "Gross_Sources", "Gross_Sinks", "Net_Emissions", "Closure_Error",
    "Flux_FD", "Flux_NFC", "Flux_FR", "Flux_NFR", "Flux_CAL", "Flux_WHp",
    "Flux_Clearing", "Flux_Abandonment", "Flux_Harvest_Net",
    "Flux_Harvest_SoilSlash", "Flux_Harvest_Regrowth", "Flux_Other",
    "Flux_Products_Total",
]

HARVEST_DIAG_KEYS = [
    # Harvest biomass diagnostics
    "harvest_requested_biomass",      # Area-implied removal or LUH2 *_bioh demand, depending on mode
    "harvest_met_biomass",            # Demand met within the LUH2 harvest footprint
    "harvest_unmet_biomass",          # Demand not satisfied after the selected harvest strategy
    "harvest_unmet_raw_biomass",      # Demand not satisfied after the selected harvest strategy
    "harvest_forced_biomass",         # Removal outside the LUH2 footprint in expansion mode

    # Harvest area diagnostics
    "harvest_luh2_area_frac",         # LUH2 *_harv fraction of whole grid cell
    "harvest_luh2_area",              # LUH2 *_harv area fraction converted to absolute area
    "harvest_effective_area_frac",    # Effective modeled harvest area / cell_area
    "harvest_effective_area",         # Effective modeled harvest area inferred from biomass intensity
    "harvest_extra_area",             # Harvest area outside LUH2 footprint; zero in strict-area version

    # Physical stock-removal diagnostics
    "harvest_biomass_removed",        # Biomass actually removed from modeled current stock
    "harvest_patch_biomass_before",   # Current biomass stock available before removal
    "harvest_loss_biomass",           # Same as removed biomass, kept for debugging old output names

    # Allocation diagnostics
    "harvest_to_products",            # Product-pool input from actually removed harvest biomass
    "harvest_to_soil",                # SR input from biomass slash plus soil disturbance
    "harvest_to_SR_from_biomass",     # Biomass slash fraction sent to SR
    "harvest_to_SR_from_soil",        # Slow-soil disturbance R sent to SR
    "harvest_forced_to_products",     # Forced unmet product/removal demand sent to product pools
    "harvest_forced_to_soil",         # Forced unmet sent to SR; should remain zero in product-only forced mode

    # Equation-level diagnostics
    "harvest_beta_delta_sum",
    "harvest_sigma_delta_sum",
    "harvest_delta_B_h",
    "harvest_delta_SS_h",
    "harvest_R",
]

FORMAL_AREA_DIAG_KEYS = [
    "area_clearing", "area_abandonment", "area_other", "area_harvest",
    "area_deforestation",
] + HARVEST_DIAG_KEYS

FORMAL_OUTPUT_KEYS = FORMAL_FLUX_KEYS + FORMAL_AREA_DIAG_KEYS

class LULCCSimulator:
    def __init__(
        self,
        pft_path: str,
        config_path: str = "config.yml",
        LULC_path: str = "states.nc",
        trans_path: str = "transitions.nc",
        *,
        experiment_path: Optional[str] = None,
        lat_slice: slice = None,
        lon_slice: slice = None,
        area_unit: str = "ha", # using ha for Qing et al. (2024) data. BECAREFUL CHANGING TO OTHER UNITS
    ):
        """
        Parameters
        ----------
        config_path : str
            YAML config file for parameters.
        LULC_path : str
            LUH2 states file (e.g., states.nc).
        trans_path : str
            LUH2 transitions file (e.g., transitions.nc).
        pft_path : str
            PFT map file (required).
        experiment_path : str or None
            Optional experiment YAML overlaid on the base config.
        area_unit : str
            "m2" (default), "ha", or "km2".
        """
        self.params = ParameterLoader(
            config_path,
            experiment_path=experiment_path,
        )
        # --- Load LUH2 + PFT through a dedicated loader (I/O + format normalization) ---
        # Server version: only the requested lat/lon slice is loaded into memory.
        self.loader = FileLoader()
        self.io_lat_slice = lat_slice
        self.io_lon_slice = lon_slice

        self.state_data, state_meta = self.loader.load_luh2_dataset(
            LULC_path,
            lat_slice=lat_slice,
            lon_slice=lon_slice,
        )
        self.trans_data, trans_meta = self.loader.load_luh2_dataset(
            trans_path,
            lat_slice=lat_slice,
            lon_slice=lon_slice,
        )

        # Keep original coordinate arrays & flags (for index mapping)
        self.state_lat = state_meta.lat
        self.state_lon = state_meta.lon
        self.trans_lat = trans_meta.lat
        self.trans_lon = trans_meta.lon
        
        self.state_lat_desc = state_meta.lat_desc
        self.trans_lat_desc = trans_meta.lat_desc

        self.state_lat_name = state_meta.lat_name
        self.state_lon_name = state_meta.lon_name
        self.actual_ds_lat_is_desc = state_meta.lat_desc

        self.trans_lat_name = trans_meta.lat_name
        self.trans_lon_name = trans_meta.lon_name
        self.actual_trans_lat_is_desc = trans_meta.lat_desc
        self.state_lat = state_meta.lat
        self.state_lon = state_meta.lon

        # Keep lon conventions separately for states and transitions
        self.state_lon_0_360 = state_meta.lon_0_360
        self.trans_lon_0_360 = trans_meta.lon_0_360

        # Build local internal grid centers + area for the loaded slice.
        self._build_grid_and_area(
            state_meta.internal_lat,
            state_meta.internal_lon,
            unit=area_unit,
        ) # produce self._area_grid[i,j]
        # Backward-compat alias used in the rest of the code
        # PFT map. FileLoader will infer lat/lon conventions and align to the internal grid.
        self.pft_map = self.loader.load_pft_map(
            pft_path,
            pft_var="maxvegetfrac",
            target_lat_asc=self._latitudes,
            target_lon_180=self._longitudes,
            target_lon_res=self.lon_res,
            lat_slice=lat_slice,
            lon_slice=lon_slice,
        )

        
        self.area_unit = area_unit

    # ---------------------------------------------------------------------
    # Loading-related helpers (grid + index normalization)
    # ---------------------------------------------------------------------
    def _build_grid_and_area(self, internal_lat, internal_lon, *, unit: str = "ha") -> None:
        """
        Build local internal grid centers and per-cell area for the loaded slice.

        The server FileLoader already maps the global internal lat/lon slice to
        the source files and returns local internal coordinates. Therefore this
        function must not rebuild a global grid from source coordinate length.
        """
        self._latitudes = np.asarray(internal_lat, dtype=float)
        self._longitudes = np.asarray(internal_lon, dtype=float)

        if self._latitudes.size >= 2:
            self.lat_res = abs(float(np.median(np.diff(self._latitudes))))
        else:
            # Fallback for a one-row slice: use metadata/source spacing if available.
            self.lat_res = abs(float(np.median(np.diff(np.asarray(self.state_lat, dtype=float)))))

        if self._longitudes.size >= 2:
            self.lon_res = abs(float(np.median(np.diff(self._longitudes))))
        else:
            self.lon_res = abs(float(np.median(np.diff(np.asarray(self.state_lon, dtype=float)))))

        nlat = self._latitudes.size
        nlon = self._longitudes.size

        R = 6371000.0
        dphi = np.deg2rad(self.lat_res)
        dlmb = np.deg2rad(self.lon_res)
        LAT = np.repeat(self._latitudes[:, None], nlon, axis=1)
        area_m2 = (R**2) * dphi * dlmb * np.cos(np.deg2rad(LAT))

        u = (unit or "m2").lower()
        if u in ("m2", "m^2"):
            area = area_m2
        elif u == "ha":
            area = area_m2 / 1e4
        elif u in ("km2", "km^2"):
            area = area_m2 / 1e6
        else:
            raise ValueError(f"Unsupported area unit: {unit}")

        self._area_grid = area.astype(np.float64)

    def _normalize_ij(self, lat_idx: int, lon_idx: int) -> Tuple[int, int]:
        """
        Normalize any possibly-flattened index into a legal 2D (i, j).
        Prefer pft_map shape if available; otherwise fall back to area grid shape.
        """
        if getattr(self, "pft_map", None) is not None and getattr(self.pft_map, "ndim", 0) >= 2:
            nlat, nlon = self.pft_map.shape[:2]
        else:
            nlat, nlon = self._area_grid.shape

        i, j = int(lat_idx), int(lon_idx)
        if i >= nlat or j >= nlon:
            flat = i if i >= nlat else j
            i, j = divmod(int(flat), int(nlon))

        if not (0 <= i < nlat and 0 <= j < nlon):
            raise IndexError(f"[normalize] out of range: i={i}, j={j}, nlat={nlat}, nlon={nlon}")
        return i, j

    def _get_cell_pft_grid(self, lat_idx: int, lon_idx: int) -> np.ndarray:
        i, j = self._normalize_ij(lat_idx, lon_idx)

        n_pft = self.params.n_pft

        if self.pft_map is None:
            out = np.zeros(n_pft, dtype=np.float64)
            out[0] = 1.0
            return out

        vals = self.pft_map[i, j]

        # multi-PFT fraction map: vals.shape = (n_pft,)
        if np.ndim(vals) == 1:
            out = np.asarray(vals, dtype=np.float64).copy()
            if out.size != n_pft:
                raise ValueError(f"pft_map last dim = {out.size}, but params.n_pft = {n_pft}")

        # legacy dominant-pft map
        else:
            out = np.zeros(n_pft, dtype=np.float64)
            if not np.isfinite(vals):
                out[0] = 1.0
                return out

            p = int(vals)
            if 1 <= p <= n_pft:
                out[p - 1] = 1.0
            else:
                out[0] = 1.0

        out[~np.isfinite(out)] = 0.0
        out[out < 0.0] = 0.0

        s = out.sum()
        if s > 0.0:
            out /= s
        else:
            out[0] = 1.0

        return out

    def _ds_lat_idx(self, i: int, which: str = "state") -> int:
        if which == "state":
            n = self.state_data.dimensions[self.state_lat_name].size
            is_desc = self.actual_ds_lat_is_desc
        else:
            n = self.trans_data.dimensions[self.trans_lat_name].size
            is_desc = self.actual_trans_lat_is_desc
        
        return (n - 1 - i) if is_desc else i

    def _wrap_lon_180(self, lon):
        return ((np.asarray(lon) + 180.0) % 360.0) - 180.0

    def _get_lon_axis_from_var(self, lon_var) -> np.ndarray:
        """
        Extract a 1D longitude axis from a lon coordinate variable.
        Supports both 1D lon(lon) and 2D lon(lat,lon).
        """
        lon_arr = np.asarray(lon_var[:], dtype=float)

        if lon_arr.ndim == 1:
            return lon_arr

        if lon_arr.ndim == 2:
            row0 = lon_arr[0, :]
            if np.allclose(lon_arr, row0[None, :], equal_nan=True):
                return row0

            col0 = lon_arr[:, 0]
            if np.allclose(lon_arr, col0[:, None], equal_nan=True):
                return col0

            std0 = np.nanstd(lon_arr, axis=0).mean()
            std1 = np.nanstd(lon_arr, axis=1).mean()

            if std1 < std0:
                return lon_arr[0, :]
            else:
                return lon_arr[:, 0]

        raise ValueError(f"Unsupported lon variable ndim={lon_arr.ndim}")

    def _ds_lon_idx(self, j: int, which: str = "state") -> int:
        """
        In the server memory-loader version, longitude variables have already
        been sliced and reordered/mapped to local internal indices. Therefore the
        model-local longitude index is the dataset-local longitude index.
        """
        return int(j)

    def indices_for_bbox(self, lat_min=-90.0, lat_max=90.0, lon_min=73.0, lon_max=135.0):
        """
        Return (lat_slice, lon_slice) for a given bounding box using internal grid centers.
        """
        lat = self._latitudes
        lon = self._longitudes

        def wrap180(x):
            return ((x + 180.0) % 360.0) - 180.0

        lon_min = wrap180(lon_min)
        lon_max = wrap180(lon_max)

        i0 = int(np.searchsorted(lat, lat_min, side="left"))
        i1 = int(np.searchsorted(lat, lat_max, side="right"))
        j0 = int(np.searchsorted(lon, lon_min, side="left"))
        j1 = int(np.searchsorted(lon, lon_max, side="right"))

        i0 = max(0, min(i0, lat.size))
        i1 = max(0, min(i1, lat.size))
        j0 = max(0, min(j0, lon.size))
        j1 = max(0, min(j1, lon.size))
        if i1 <= i0 or j1 <= j0:
            raise ValueError(f"Empty slice: lat=({lat_min},{lat_max}), lon=({lon_min},{lon_max})")
        return slice(i0, i1), slice(j0, j1)

    def _get_initial_fractions(self, lat_idx: int = None, lon_idx: int = None, year_idx: int = 0) -> Dict[str, float]:
        """Get initial area for 4+1 LULCs from states.nc (returns fractional area)."""
        i, j = self._normalize_ij(lat_idx, lon_idx)
        frac_dict = {"v": 0.0, "s": 0.0, "p": 0.0, "c": 0.0, "U": 0.0}
        try:
            ii = self._ds_lat_idx(i, "state")
            jj = self._ds_lon_idx(j, "state")
            val = self.state_data.variables['pastr'][0, ii, jj]
            
            p_val = float(val) if not np.ma.is_masked(val) else np.nan
            #print(f"DEBUG (i={i}, j={j}) -> DS(ii={ii}, jj={jj}), pastr={p_val:.4f}")
        except Exception as e:
            print(f"DEBUG Error at (i={i}, j={j}): {e}")

        for LUH2_type in range(1, 13):
            lulc_type = self.params.get_LULC(LUH2_type)
            var_name = self.params.get_LUH2Type(LUH2_type)

            try:
                v = self.state_data.variables[var_name][year_idx, ii, jj]
                frac = float(v.filled(np.nan) if hasattr(v, "filled") else v)
            except Exception:
                frac = np.nan

            # avoid nan
            if not np.isfinite(frac) or frac < 0.0:
                frac = 0.0
            elif frac > 1.0:
                frac = 1.0

            if lulc_type in frac_dict:
                frac_dict[lulc_type] += frac

        return frac_dict

    def _parse_transitions(self, year_idx: int, lat_idx: int, lon_idx: int) -> Dict[str, float]:
        """
        Use transitions.nc variables *_frac (0..1) and return aggregated transition
        FRACTIONS by coarse LULC class (v/s/p/c), not absolute area.

        Returned values are fractions of the whole grid cell per year.
        Absolute area is computed later in run_simulation() as:
            area_abs = trans_frac * cell_area

        Missing channels are treated as 0.
        """
        i, j = self._normalize_ij(lat_idx, lon_idx)
        transitions: Dict[str, float] = {}

        ii = self._ds_lat_idx(i, "trans")
        jj = self._ds_lon_idx(j, "trans")

        for LUH2_type_src in range(1, 13):
            for LUH2_type_dst in range(1, 13):
                src_name = self.params.get_LUH2Type(LUH2_type_src)
                dst_name = self.params.get_LUH2Type(LUH2_type_dst)
                key_base = f"{src_name}_to_{dst_name}"
                lulc_key = f"{self.params.get_LULC(LUH2_type_src)}_to_{self.params.get_LULC(LUH2_type_dst)}"

                var_name = f"{key_base}_frac" # be careful here when using 0.25 LUH2 transition
                var = self.trans_data.variables.get(var_name)
                if var is None:
                    var_name = key_base
                    var = self.trans_data.variables.get(var_name)
                if var is None:
                    continue
                try:
                    v = var[year_idx, ii, jj]
                    frac = float(v.filled(np.nan) if hasattr(v, "filled") else v)
                except Exception:
                    frac = np.nan

                if not np.isfinite(frac) or frac < 0.0:
                    frac = 0.0
                elif frac > 1.0:
                    frac = 1.0

                frac_val = frac
                if frac_val > 0.0:
                    transitions[lulc_key] = transitions.get(lulc_key, 0.0) + frac_val
        return transitions

    def _parse_wood_harvest(
        self,
        year_idx: int,
        lat_idx: int,
        lon_idx: int,
    ) -> Dict[str, Dict[str, float]]:
        """Parse the LUH2 harvest channels required by the active experiment."""
        i, j = self._normalize_ij(lat_idx, lon_idx)
        ii = self._ds_lat_idx(i, "trans")
        jj = self._ds_lon_idx(j, "trans")

        families = ["primf", "primn", "secmf", "secyf", "secnf"]
        use_bioh = self.params.harvest_use_bioh
        out: Dict[str, Dict[str, float]] = {}

        for family in families:
            area_var = self.trans_data.variables.get(f"{family}_harv")
            biomass_var = (
                self.trans_data.variables.get(f"{family}_bioh")
                if use_bioh
                else None
            )

            area_frac = 0.0
            if area_var is not None:
                try:
                    value = area_var[year_idx, ii, jj]
                    value = value.filled(np.nan) if hasattr(value, "filled") else value
                    area_frac = float(value)
                except Exception:
                    area_frac = np.nan

            if not np.isfinite(area_frac) or area_frac < 0.0:
                area_frac = 0.0
            elif area_frac > 1.0:
                area_frac = 1.0

            biomass_kgC = 0.0
            if biomass_var is not None:
                try:
                    value = biomass_var[year_idx, ii, jj]
                    value = value.filled(np.nan) if hasattr(value, "filled") else value
                    biomass_kgC = float(value)
                except Exception:
                    biomass_kgC = np.nan

            if not np.isfinite(biomass_kgC) or biomass_kgC < 0.0:
                biomass_kgC = 0.0

            if self.params.harvest_mode == "area":
                active = area_frac > 0.0
            else:
                active = area_frac > 0.0 or biomass_kgC > 0.0

            if active:
                out[family] = {
                    "area_frac": float(area_frac),
                    "biomass_kgC": float(biomass_kgC),
                }

        return out

    def run_simulation_grid(self, years: int, out_nc: str = "summary.nc",
                        lat_slice: slice = None, lon_slice: slice = None,
                        start_year_idx: int = 0, sync_every: int = 200):
        """
        traverse grids within slices. run while writting results
        - create NixNj dimension
        - chunksizes=(T,1,1)
        - after each sync_every clear ache
        """
        import time
        t0 = time.time()
    
        # resolve slice
        nlat_all, nlon_all = self._area_grid.shape
        if lat_slice is None: lat_slice = slice(0, nlat_all, 1)
        if lon_slice is None: lon_slice = slice(0, nlon_all, 1)
    
        i0 = lat_slice.start or 0
        i1 = lat_slice.stop  or nlat_all
        j0 = lon_slice.start or 0
        j1 = lon_slice.stop  or nlon_all
    
        Ni = i1 - i0
        Nj = j1 - j0
        if Ni <= 0 or Nj <= 0:
            raise ValueError(f"Empty slice: lat_slice=({i0},{i1}), lon_slice=({j0},{j1})")
    
        T = years + 1
        print(f"[grid] Ni={Ni}, Nj={Nj}, T={T}, lat_slice=({i0},{i1}), lon_slice=({j0},{j1})")
    
        # create deminsion
        latitudes  = self._latitudes[i0:i1].astype(np.float32)
        longitudes = self._longitudes[j0:j1].astype(np.float32)
        time_vec   = (np.arange(T) + start_year_idx).astype(np.int32)
    
        # create output file and write variables
        chunks = (min(T, 64), 1, 1)
        kwargs = dict(zlib=True, complevel=4, fill_value=np.nan, chunksizes=chunks)
    
        ds = nc.Dataset(out_nc, "w")
        try:
            ds.experiment_name = self.params.experiment_name
            ds.harvest_mode = self.params.harvest_mode
            ds.harvest_use_bioh = int(self.params.harvest_use_bioh)
            ds.harvest_allow_expansion = int(
                self.params.harvest_allow_expansion
            )

            ds.createDimension("time", T)
            ds.createDimension("lat",  Ni)
            ds.createDimension("lon",  Nj)
    
            vtime = ds.createVariable("time", "i4", ("time",))
            vlat  = ds.createVariable("lat",  "f4", ("lat",))
            vlon  = ds.createVariable("lon",  "f4", ("lon",))
            vtime[:] = time_vec
            vlat[:] = latitudes
            vlon[:]  = longitudes
            vlat.units  = "degrees_north";  vlat.standard_name = "latitude"
            vlon.units  = "degrees_east";   vlon.standard_name = "longitude"
            vtime.long_name = "LUH2 index year"
    
            # State variables kept for pool/closure checking
            state_vars = ["biomass_total", "soil_total", "P1", "P10", "P100", "atmosphere"]
            out_vars = {}
            for name in state_vars + FORMAL_OUTPUT_KEYS:
                out_vars[name] = ds.createVariable(name, "f8", ("time", "lat", "lon"), **kwargs)
            # Minimal metadata for downstream scripts
            out_vars["Gross_Sources"].long_name = "Gross positive LULCC carbon flux to atmosphere"
            out_vars["Gross_Sinks"].long_name = "Gross negative LULCC carbon uptake from atmosphere"
            out_vars["Net_Emissions"].long_name = "ELUC = Gross_Sources + Gross_Sinks"
            out_vars["Closure_Error"].long_name = "Net_Emissions - annual atmosphere increment"
            out_vars["Flux_FD"].long_name = "LUCE category: deforestation"
            out_vars["Flux_NFC"].long_name = "LUCE category: non-forest conversion"
            out_vars["Flux_FR"].long_name = "LUCE category: reforestation and harvest regrowth"
            out_vars["Flux_NFR"].long_name = "LUCE category: non-forest reconstruction"
            out_vars["Flux_CAL"].long_name = "LUCE category: conversions between anthropogenic land"
            out_vars["Flux_WHp"].long_name = "LUCE category: wood harvest product emissions"
            
    
            # done flag for checking if simulation for grid is done
            vdone = ds.createVariable("done", "i1", ("lat","lon"), zlib=True, complevel=1,
                                      chunksizes=(min(Ni,64), min(Nj,64)))
            vdone[:, :] = 0
    
            # main loop
            total = Ni * Nj
            done  = 0
            ema   = None; alpha = 0.25
            def _eta(sec):
                sec = int(max(0, round(sec))); h, r = divmod(sec,3600); m,s = divmod(r,60)
                return f"{h}h{m:02d}m{s:02d}s" if h else (f"{m}m{s:02d}s" if m else f"{s}s")
    
            for ii, i in enumerate(range(i0, i1)):
                for jj, j in enumerate(range(j0, j1)):
                    # skip ocean/void grid
                    fracs0 = self._get_initial_fractions(i, j, year_idx=start_year_idx)
                    if sum(fracs0.values()) <= 0.0:
                        vdone[ii, jj] = 1
                        done += 1
                        continue
                    # skip if done flag = 1
                    if vdone[ii, jj] == 1:
                        done += 1
                        continue
    
                    t1 = time.time()
                    # TRAVERSE TIME
                    yearly = self.run_simulation(years, i, j, start_year_idx=start_year_idx)
    
                    upto = min(T, len(yearly))
                    for t in range(upto):
                        rec = yearly[t]
                        for name, var in out_vars.items():
                            var[t, ii, jj] = float(rec.get(name, 0.0))

                    vdone[ii, jj] = 1  # done flag
    
                    # progress bar
                    dt = time.time() - t1
                    ema = dt if ema is None else (alpha*dt + (1-alpha)*ema)
                    done += 1
                    pct = 100.0 * done / total
                    spent = time.time() - t0
                    eta   = ema * max(0, total - done) if ema else 0.0
                    barw  = 28; bar = int(barw * pct / 100.0)
                    print(f"\r[{done:6d}/{total:6d}] {pct:5.1f}% |{'█'*bar}{'.'*(barw-bar)}| "
                          f"spent {_eta(spent)} ETA {_eta(eta)} (i={i},j={j})",
                          end="", flush=True)
    
                    # synchronization  
                    if sync_every and (done % sync_every == 0):
                        ds.sync()
    
            print("\n[grid] DONE")
            ds.history = (getattr(ds, "history", "") + f"\nstreaming write; created at {time.ctime()}").strip()
        finally:
            ds.close()

    # CORE SCRIPT. grid scale simualtion
    def run_simulation(self, years: int, lat_idx: int, lon_idx: int, start_year_idx: int = 0):
        # Initialization
        dt = 1.0 #running at yearly timestep
        n_time = years + 1
        cell_area = float(self._area_grid[lat_idx, lon_idx]) # area for simulation grid cell
        pft_grid = self._get_cell_pft_grid(lat_idx, lon_idx) # pft fractions in a grid cell
        LULC_frac = self._get_initial_fractions(lat_idx, lon_idx, year_idx=start_year_idx) # fraction area for each LULC in a grid
        if sum(LULC_frac.values()) <= 0.0:
            empty_record = {
                "biomass_total": 0.0,
                "soil_total": 0.0,
                "P1": 0.0, "P10": 0.0, "P100": 0.0,
                "atmosphere": 0.0,
            }
            for k in FORMAL_OUTPUT_KEYS:
                empty_record[k] = 0.0
            return [empty_record.copy() for _ in range(years + 1)]
        
        base_year = int(
            self.params.config.get("dynamic_carbon_density", {}).get("base_year", 850)
        )
        calendar_year0 = base_year + start_year_idx
        self.params.set_density_year(calendar_year0)
        C_bar = initialize_Cbar(LULC_frac, cell_area, pft_grid, self.params, n_time=n_time) # equilibrium carbon content for land cover of each pft
        Delta = initialize_Delta(n_time=n_time) # excess pool
        frac_area = initialize_frac_area(LULC_frac, pft_grid) # fraction for each LULC of each pft, change with time

        yearly = [] #save results
        Atmos = np.zeros(n_time, dtype=np.float64)

        # save initial condition
        rec0 = summarize_year(C_bar, Delta, Atmos, t_idx=0)
        for k in FORMAL_OUTPUT_KEYS:
            rec0[k] = 0.0
        yearly.append(rec0)

        # traverse time 
        for year in range(years):
            t_cur = year
            t_next = year + 1
            Delta[t_next, :, :, :, :] = Delta[t_cur, :, :, :, :]
            Atmos[t_next] = Atmos[t_cur]

            calendar_year = base_year + start_year_idx + year + 1
            self.params.set_density_year(calendar_year)

            if self.params.use_dynamic_density:
                refresh_Cbar_from_frac_area(
                    C_bar=C_bar,
                    t_idx=t_next,
                    frac_area=frac_area,
                    cell_area=cell_area,
                    params=self.params,
                )
            else:
                C_bar[t_next, :, :, :] = C_bar[t_cur, :, :, :]
            # Only formal output diagnostics are requested here.
            # events.py uses _diag_add(), so old low-level diagnostics are skipped.
            diag = {k: 0.0 for k in FORMAL_OUTPUT_KEYS}
            # get LULCC from XX to YY, and area in year#
            # starting year can be set in here, e.g., 1000 year+: transitions = self._parse_transitions(year + 1000, lat_idx, lon_idx)  
            # get transition matrix
            transitions = self._parse_transitions(year + start_year_idx, lat_idx, lon_idx)  
            # traverse land use change types in year#
            for trans_key, trans_frac in transitions.items():
                src, dst = trans_key.split('_to_')
                area_abs = trans_frac * cell_area

                if area_abs <= 0.0:
                    continue
                # determine LULCC type. Clearing, Abandonment, Harvest, Others
                if src in ['v', 's'] and dst in ['p', 'c']:  # Clearing
                    apply_clearing(
                        area=area_abs,
                        src=src,
                        dst=dst,
                        cell_area=cell_area,
                        pft_grid=pft_grid,
                        params=self.params,
                        C_bar=C_bar,
                        Delta=Delta,
                        LULC_frac=LULC_frac,
                        frac_area=frac_area,
                        t_idx=t_next,
                        diag=diag
                    )
                    continue
                    
                elif src in ['c', 'p'] and dst == 's':  # Abandonment
                    apply_abandonment(
                        area=area_abs,
                        src=src,
                        dst=dst,
                        cell_area=cell_area,
                        pft_grid=pft_grid,
                        params=self.params,
                        C_bar=C_bar,
                        Delta=Delta,
                        LULC_frac=LULC_frac,
                        frac_area=frac_area,
                        t_idx=t_next,
                        diag=diag
                    )
                    continue
                
                elif src in ['c', 'p'] and dst in ['c', 'p'] and src != dst:  # Others
                    apply_others(
                        area=area_abs,
                        src=src,
                        dst=dst,
                        cell_area=cell_area,
                        pft_grid=pft_grid,
                        params=self.params,
                        C_bar=C_bar,
                        Delta=Delta,
                        LULC_frac=LULC_frac,
                        frac_area=frac_area,
                        t_idx=t_next,
                        diag=diag
                    )
                    continue

                else:
                    continue
                
            # LUH2 wood harvest is handled separately from state transitions.
            wood_harv = self._parse_wood_harvest(
                year + start_year_idx,
                lat_idx,
                lon_idx,
            )

            for harvest_type, harvest_input in wood_harv.items():
                area_frac = float(harvest_input.get("area_frac", 0.0))
                biomass_kgC = float(harvest_input.get("biomass_kgC", 0.0))

                apply_harvest_luh2(
                    mode=self.params.harvest_mode,
                    allow_expansion=self.params.harvest_allow_expansion,
                    harvest_type=harvest_type,
                    area_frac=area_frac,
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
                    t_idx=t_next,
                    diag=diag,
                )

            # End-of-time-step decay / relaxation
            Atmos[t_next], diag_relax = relax_one_year_blue_with_diag(
                Delta=Delta,
                Atmos=Atmos[t_next],
                params=self.params,
                t_idx=t_next,
                dt=dt,
            )

            # report yearly results
            for k, v in diag_relax.items():
                diag[k] += v

            # Harvest diagnostics are accumulated inside events.py via _diag_add().
            diag["Closure_Error"] = diag["Net_Emissions"] - (Atmos[t_next] - Atmos[t_cur])
            rec = summarize_year(C_bar, Delta, Atmos, t_idx=t_next)
            for k, v in diag.items():
                rec[k] = v
            yearly.append(rec)
        return yearly
    
    def run_diagnostic_scan(self, lat_slice: slice, lon_slice: slice):
        i0, i1 = lat_slice.start, lat_slice.stop
        j0, j1 = lon_slice.start, lon_slice.stop

        print(f"--- 诊断开始 ---")
        print(f"目标范围: Lat[{i0}:{i1}], Lon[{j0}:{j1}]")

        count_show = 0
        valid_land_cells = 0
        zero_cells = 0
        lon_axis = self._get_lon_axis_from_var(self.state_data.variables[self.state_lon_name])
        lat_axis = self._get_lat_axis_from_var(self.state_data.variables[self.state_lat_name])

        for i in range(i0, i1):
            for j in range(j0, j1):
                ii = self._ds_lat_idx(i, "state")
                jj = self._ds_lon_idx(j, "state")

                lon_internal = self._longitudes[j]
                lat_internal = self._latitudes[i]

                lon_ds = float(lon_axis[jj])
                lat_ds = float(lat_axis[ii])

                fracs = self._get_initial_fractions(i, j, year_idx=0)
                s = sum(fracs.values())

                if s > 0:
                    valid_land_cells += 1
                else:
                    zero_cells += 1

                if count_show < 20:
                    print(
                        f"(i={i},j={j}) "
                        f"internal=({lat_internal:.1f},{lon_internal:.1f}) -> "
                        f"ds=({lat_ds:.1f},{lon_ds:.1f}) "
                        f"LULCsum={s:.3f}"
                    )
                    count_show += 1

        print(f"有效格点: {valid_land_cells}")
        print(f"LULC为0格点: {zero_cells}")
        print(f"--- 诊断结束 ---")
    
    def _get_lat_axis_from_var(self, lat_var) -> np.ndarray:
        """
        Extract a 1D latitude axis from a lat coordinate variable.
        Supports both 1D lat(lat) and 2D lat(lat,lon).
        """
        lat_arr = np.asarray(lat_var[:], dtype=float)

        if lat_arr.ndim == 1:
            return lat_arr

        if lat_arr.ndim == 2:
            # 常见情况1：每一列都一样，纬度沿第1维变化
            col0 = lat_arr[:, 0]
            if np.allclose(lat_arr, col0[:, None], equal_nan=True):
                return col0

            # 常见情况2：每一行都一样，纬度沿第2维变化
            row0 = lat_arr[0, :]
            if np.allclose(lat_arr, row0[None, :], equal_nan=True):
                return row0

            std0 = np.nanstd(lat_arr, axis=0).mean()
            std1 = np.nanstd(lat_arr, axis=1).mean()

            if std0 < std1:
                return lat_arr[:, 0]
            else:
                return lat_arr[0, :]

        raise ValueError(f"Unsupported lat variable ndim={lat_arr.ndim}")
