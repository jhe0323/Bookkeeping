"""Parameter loading and experiment-configuration management."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
import yaml


class ParameterLoader:
    def __init__(
        self,
        config_path: str = "config.yml",
        experiment_path: Optional[str] = None,
    ):
        self.config_path = self._resolve_path(config_path)
        self.base_config = self._load_config(self.config_path)

        self.experiment_path = (
            self._resolve_path(experiment_path)
            if experiment_path is not None
            else None
        )
        self.experiment_config = (
            self._load_config(self.experiment_path)
            if self.experiment_path is not None
            else {}
        )
        self.config = self._deep_merge(
            self.base_config,
            self.experiment_config,
        )

        self.pft_names = [f"PFT{i}" for i in range(1, 16)]
        self.pft_index = {name: i for i, name in enumerate(self.pft_names)}
        self.n_pft = len(self.pft_names)

        self._load_experiment_settings()

        self.dynamic_cfg = self.config.get("dynamic_carbon_density", {})
        self.use_dynamic_density = bool(self.dynamic_cfg.get("enabled", False))
        self.current_density_year = None
        self.dynamic_density = None
        self.dynamic_years = None

        if self.use_dynamic_density:
            self._load_dynamic_density()

    @staticmethod
    def _resolve_path(path: str) -> Path:
        resolved = Path(path)
        if not resolved.is_absolute():
            resolved = Path.cwd() / resolved
        return resolved.resolve()

    @staticmethod
    def _load_config(path: Path) -> Dict[str, Any]:
        try:
            with path.open("r", encoding="utf-8") as handle:
                data = yaml.safe_load(handle)
        except FileNotFoundError as exc:
            raise ValueError(f"Configuration file not found: {path}") from exc
        except yaml.YAMLError as exc:
            raise ValueError(f"YAML loading error in {path}: {exc}") from exc

        if data is None:
            return {}
        if not isinstance(data, dict):
            raise ValueError(f"Top-level YAML content must be a mapping: {path}")
        return data

    @classmethod
    def _deep_merge(
        cls,
        base: Dict[str, Any],
        override: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Recursively merge override values without mutating either input."""
        merged = deepcopy(base)
        for key, value in override.items():
            if (
                key in merged
                and isinstance(merged[key], dict)
                and isinstance(value, dict)
            ):
                merged[key] = cls._deep_merge(merged[key], value)
            else:
                merged[key] = deepcopy(value)
        return merged

    def _load_experiment_settings(self) -> None:
        experiment_cfg = self.config.get("experiment", {})
        harvest_cfg = self.config.get("harvest", {})

        self.experiment_name = str(
            experiment_cfg.get("name", "harvest_area_driven")
        )
        self.harvest_mode = str(harvest_cfg.get("mode", "area"))
        self.harvest_use_bioh = bool(
            harvest_cfg.get("demand_source", {}).get(
                "use_bioh",
                self.harvest_mode == "biomass_demand",
            )
        )
        self.harvest_use_luh2_area = bool(
            harvest_cfg.get("area_constraint", {}).get(
                "use_luh2_harvest_area",
                True,
            )
        )
        self.harvest_allow_expansion = bool(
            harvest_cfg.get("allow_expansion", {}).get("enabled", False)
        )

        supported_modes = {"area", "biomass_demand"}
        if self.harvest_mode not in supported_modes:
            raise ValueError(
                f"Unsupported harvest.mode={self.harvest_mode!r}; "
                f"expected one of {sorted(supported_modes)}"
            )

        if not self.harvest_use_luh2_area:
            raise ValueError(
                "The current implementation requires "
                "harvest.area_constraint.use_luh2_harvest_area=true."
            )

        if self.harvest_mode == "area":
            if self.harvest_use_bioh:
                raise ValueError(
                    "Area-driven harvest must set demand_source.use_bioh=false."
                )
            if self.harvest_allow_expansion:
                raise ValueError(
                    "Area-driven harvest cannot enable allow_expansion."
                )

        if self.harvest_mode == "biomass_demand" and not self.harvest_use_bioh:
            raise ValueError(
                "Biomass-demand harvest must set demand_source.use_bioh=true."
            )

    def get_forest_pfts(self):
        """Return configured forest PFT indices in zero-based form."""
        return set(
            int(index) - 1
            for index in self.config.get("forest_pft_indices", [])
        )

    def get_carbon_density(
        self,
        p: int,
        pool_type: str,
        land_cover: str,
    ) -> float:
        if self.use_dynamic_density:
            return self._get_dynamic_carbon_density(p, pool_type, land_cover)
        return self._get_static_carbon_density(p, pool_type, land_cover)

    def _load_dynamic_density(self) -> None:
        path = Path(
            self.dynamic_cfg.get("path", "dynamic_carbon_density.parquet")
        )
        if not path.is_absolute():
            path = self.config_path.parent / path

        required = [
            "year",
            "pft_id",
            "B_v_tCha",
            "B_s_tCha",
            "B_c_tCha",
            "B_p_tCha",
            "SS_v_tCha",
            "SS_s_tCha",
            "SS_c_tCha",
            "SS_p_tCha",
        ]
        frame = pd.read_parquet(path, columns=required, engine="pyarrow")

        missing = [column for column in required if column not in frame.columns]
        if missing:
            raise ValueError(
                f"Dynamic carbon density missing columns: {missing}"
            )

        self.dynamic_density = frame.set_index(
            ["year", "pft_id"]
        ).sort_index()
        self.dynamic_years = np.asarray(
            sorted(frame["year"].unique()),
            dtype=int,
        )
        self.current_density_year = int(self.dynamic_years[0])

    def set_density_year(self, year: int) -> None:
        if not self.use_dynamic_density:
            return

        selected_year = int(year)
        minimum = int(self.dynamic_years.min())
        maximum = int(self.dynamic_years.max())
        min_policy = self.dynamic_cfg.get("min_year_policy", "clip")
        max_policy = self.dynamic_cfg.get("max_year_policy", "clip")

        if selected_year < minimum:
            if min_policy == "clip":
                selected_year = minimum
            else:
                raise ValueError(
                    f"Density year {selected_year} < available minimum {minimum}"
                )

        if selected_year > maximum:
            if max_policy == "clip":
                selected_year = maximum
            else:
                raise ValueError(
                    f"Density year {selected_year} > available maximum {maximum}"
                )

        self.current_density_year = selected_year

    def _get_dynamic_carbon_density(
        self,
        p: int,
        pool_type: str,
        land_cover: str,
    ) -> float:
        if land_cover == "U":
            return 0.0

        column_map = {
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
        if key not in column_map:
            raise KeyError(f"Unsupported dynamic density key: {key}")

        year = int(self.current_density_year)
        pft_id = int(p) + 1
        return float(
            self.dynamic_density.loc[(year, pft_id), column_map[key]]
        )

    def _get_static_carbon_density(
        self,
        p: int,
        pool_type: str,
        land_cover: str,
    ) -> float:
        carbon_density = self.config["carbon_density"]
        pft_key = "PFT_U" if land_cover == "U" else self.get_pft_name(p)
        try:
            return float(carbon_density[pft_key][pool_type][land_cover])
        except KeyError as exc:
            raise KeyError(
                "Missing density: "
                f"pft={pft_key}, pool={pool_type}, cover={land_cover}"
            ) from exc

    def get_LULC(self, pft: int) -> str:
        return self.config["LUH2toLULC"][
            self.config["numtoLUH2Type"][pft]
        ]

    def get_LUH2Type(self, pft: int) -> str:
        return self.config["numtoLUH2Type"][pft]

    def get_clearing_param(self, p: int, param_name: str) -> float:
        return self.config["clearing_param"][self.get_pft_name(p)][param_name]

    def get_abandonment_param(self, p: int, param_name: str) -> float:
        return self.config["abandonment_param"][self.get_pft_name(p)][param_name]

    def get_harvest_param(self, p: int, param_name: str) -> float:
        return self.config["harvest_param"][self.get_pft_name(p)][param_name]

    def get_pft_name(self, p: int) -> str:
        if not (0 <= p < self.n_pft):
            raise IndexError(f"p out of range: {p}, n_pft={self.n_pft}")
        return self.pft_names[p]

    def get_n_pft(self) -> int:
        return self.n_pft
