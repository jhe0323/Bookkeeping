"""Parameter load and management"""
# parameter_loader.py
# -*- coding: UTF-8 -*-
from pathlib import Path
from typing import Dict, Any
import yaml
import pandas as pd
import numpy as np

class ParameterLoader:
    def __init__(self, config_path: str = "config.yml"):
        self.config_path = Path(config_path)
        if not self.config_path.is_absolute():
            self.config_path = Path.cwd() / self.config_path
        self.config = self._load_config(str(self.config_path))
        # construct pft array info. assuming 11 pfts (PFT1, PFT2, ..., PFT11)
        self.pft_names = [f"PFT{i}" for i in range(1, 16)]
        self.pft_index = {name: i for i, name in enumerate(self.pft_names)}
        self.n_pft = len(self.pft_names)
        
        self.dynamic_cfg = self.config.get("dynamic_carbon_density", {})
        self.use_dynamic_density = bool(self.dynamic_cfg.get("enabled", False))

        self.current_density_year = None
        self.dynamic_density = None
        self.dynamic_years = None

        if self.use_dynamic_density:
            self._load_dynamic_density()

    #Load YAML file
    def _load_config(self, path: str) -> Dict[str, Any]:
        try:
            full_path = Path(path)
            with open(full_path, 'r', encoding='utf-8') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            raise ValueError(f"File {path} Not Find")
        except yaml.YAMLError as e:
            raise ValueError(f"File loading error: {e}")

    def get_forest_pfts(self):
    # convert forest pft to 0-based
        return set(int(i) - 1 for i in self.config.get("forest_pft_indices", []))

    # get carbon density for specific PFT, carbon pool, and LULC type
    def get_carbon_density(self, p: int, pool_type: str, land_cover: str) -> float:
        if self.use_dynamic_density:
            return self._get_dynamic_carbon_density(p, pool_type, land_cover)

        return self._get_static_carbon_density(p, pool_type, land_cover)
    
    def _load_dynamic_density(self):
        p = Path(self.dynamic_cfg.get("path", "dynamic_carbon_density.parquet"))

        if not p.is_absolute():
            p = self.config_path.parent / p

        required = [
            "year", "pft_id",
            "B_v_tCha", "B_s_tCha", "B_c_tCha", "B_p_tCha",
            "SS_v_tCha", "SS_s_tCha", "SS_c_tCha", "SS_p_tCha",
        ]

        # Only read columns actually used by the model.
        # C_bar_tCha, Rveg and Rsoil are intentionally not loaded.
        df = pd.read_parquet(p, columns=required, engine="pyarrow")

        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Dynamic carbon density missing columns: {missing}")

        self.dynamic_density = df.set_index(["year", "pft_id"]).sort_index()
        self.dynamic_years = np.asarray(sorted(df["year"].unique()), dtype=int)

        self.current_density_year = int(self.dynamic_years[0])


    def set_density_year(self, year: int):
        """
        Set current calendar year for dynamic carbon density.
        Static mode does nothing.
        """
        if not self.use_dynamic_density:
            return

        y = int(year)
        ymin = int(self.dynamic_years.min())
        ymax = int(self.dynamic_years.max())

        min_policy = self.dynamic_cfg.get("min_year_policy", "clip")
        max_policy = self.dynamic_cfg.get("max_year_policy", "clip")

        if y < ymin:
            if min_policy == "clip":
                y = ymin
            else:
                raise ValueError(f"Density year {y} < available minimum year {ymin}")

        if y > ymax:
            if max_policy == "clip":
                y = ymax
            else:
                raise ValueError(f"Density year {y} > available maximum year {ymax}")

        self.current_density_year = y


    def _get_dynamic_carbon_density(self, p: int, pool_type: str, land_cover: str) -> float:
        if land_cover == "U":
            return 0.0

        col_map = {
            ("Biomass", "v"): "B_v_tCha",
            ("Biomass", "s"): "B_s_tCha",
            ("Biomass", "c"): "B_c_tCha",
            ("Biomass", "p"): "B_p_tCha",
            ("Soil", "v"): "SS_v_tCha",
            ("Soil", "s"): "SS_s_tCha",
            ("Soil", "c"): "SS_c_tCha",
            ("Soil", "p"): "SS_p_tCha",
        }

        key = (pool_type, land_cover)
        if key not in col_map:
            raise KeyError(f"Unsupported dynamic density key: {key}")

        year = int(self.current_density_year)
        pft_id = int(p) + 1
        col = col_map[key]

        return float(self.dynamic_density.loc[(year, pft_id), col])


    def _get_static_carbon_density(self, p: int, pool_type: str, land_cover: str) -> float:
        cd = self.config["carbon_density"]

        pool_key = pool_type
        lulc_key = land_cover

        if lulc_key == "U":
            pft_key = "PFT_U"
        else:
            pft_key = self.get_pft_name(p)

        try:
            return float(cd[pft_key][pool_key][lulc_key])
        except KeyError as e:
            raise KeyError(
                f"Missing density: pft_key={pft_key}, pool={pool_key}, lulc={lulc_key}"
            ) from e
    
    # convert LUH2 land use number e.g., 6 'c3ann', 11 'pastr' to 4+1 LULC type
    def get_LULC(self, pft: int) -> str:
        return self.config['LUH2toLULC'][self.config['numtoLUH2Type'][pft]]

    # convert LUH2 land use number e.g., 1, 2 to land use names e.g., primf, primn in LUH2
    def get_LUH2Type(self, pft: int) -> str:
        return self.config['numtoLUH2Type'][pft]

    # get parameters for land clearing
    def get_clearing_param(self, p: int, param_name: str) -> float:
        return self.config["clearing_param"][self.get_pft_name(p)][param_name]

    # get parameters for land abandonment
    def get_abandonment_param(self, p: int, param_name: str) -> float:
        return self.config["abandonment_param"][self.get_pft_name(p)][param_name]

    # get parameters for harvest
    def get_harvest_param(self, p: int, param_name: str) -> float:
        return self.config["harvest_param"][self.get_pft_name(p)][param_name]

    def get_pft_name(self, p: int) -> str:
        """0-based pft index -> config key, e.g. 0 -> 'PFT1'"""
        if not (0 <= p < self.n_pft):
            raise IndexError(f"p out of range: {p}, n_pft={self.n_pft}")
        return self.pft_names[p]

    def get_n_pft(self) -> int:
        return self.n_pft