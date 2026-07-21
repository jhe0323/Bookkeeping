"""Parameter loading and experiment-configuration management."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import yaml


_PFT_RE = re.compile(r"^PFT(\d+)$", re.IGNORECASE)


class ParameterLoader:
    """Load the model parameter YAML and an optional experiment overlay.

    PFT names and the PFT count are configuration driven.  The preferred form is::

        pft:
          names: [PFT1, PFT2, ...]
          forest_indices: [1, 2, ...]  # one-based IDs

    For backward compatibility, ``PFT<number>`` entries are inferred from
    ``carbon_density`` when the explicit ``pft`` section is absent.
    """

    def __init__(
        self,
        config_path: str = "config.yml",
        experiment_path: Optional[str] = None,
        override_config: Optional[Dict[str, Any]] = None,
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
        self.config = self._deep_merge(self.base_config, self.experiment_config)
        if override_config:
            self.config = self._deep_merge(self.config, override_config)

        self._load_pft_settings()
        self._load_experiment_settings()

        self.dynamic_cfg = self.config.get("dynamic_carbon_density", {}) or {}
        self.use_dynamic_density = bool(self.dynamic_cfg.get("enabled", False))
        self.current_density_year: Optional[int] = None
        self.dynamic_density: Optional[pd.DataFrame] = None
        self.dynamic_years: Optional[np.ndarray] = None
        if self.use_dynamic_density:
            self._load_dynamic_density()

        self.validate_parameter_completeness()

    @staticmethod
    def _resolve_path(path: Union[str, Path]) -> Path:
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
    def _deep_merge(cls, base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
        merged = deepcopy(base)
        for key, value in override.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = cls._deep_merge(merged[key], value)
            else:
                merged[key] = deepcopy(value)
        return merged

    def _load_pft_settings(self) -> None:
        pft_cfg = self.config.get("pft", {}) or {}
        explicit_names = pft_cfg.get("names")

        if explicit_names:
            names = [str(name) for name in explicit_names]
        else:
            density = self.config.get("carbon_density", {}) or {}
            numbered: List[Tuple[int, str]] = []
            for key in density:
                match = _PFT_RE.match(str(key))
                if match:
                    numbered.append((int(match.group(1)), str(key)))
            numbered.sort()
            names = [name for _, name in numbered]

        if not names:
            raise ValueError(
                "No PFT classes configured. Add pft.names or PFT<number> entries "
                "under carbon_density."
            )
        if len(set(names)) != len(names):
            raise ValueError(f"Duplicate PFT names in configuration: {names}")

        self.pft_names = names
        self.pft_index = {name: i for i, name in enumerate(names)}
        self.n_pft = len(names)

        forest_ids = pft_cfg.get(
            "forest_indices",
            self.config.get("forest_pft_indices", []),
        )
        self._forest_pfts = set()
        for value in forest_ids:
            index = int(value) - 1
            if not (0 <= index < self.n_pft):
                raise ValueError(
                    f"Forest PFT ID {value} is outside 1..{self.n_pft}."
                )
            self._forest_pfts.add(index)

    def _load_experiment_settings(self) -> None:
        experiment_cfg = self.config.get("experiment", {}) or {}
        harvest_cfg = self.config.get("harvest", {}) or {}

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
                "use_luh2_harvest_area", True
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

    def get_forest_pfts(self) -> Set[int]:
        return set(self._forest_pfts)

    def get_carbon_density(self, p: int, pool_type: str, land_cover: str) -> float:
        if self.use_dynamic_density:
            return self._get_dynamic_carbon_density(p, pool_type, land_cover)
        return self._get_static_carbon_density(p, pool_type, land_cover)

    def _load_dynamic_density(self) -> None:
        path = Path(self.dynamic_cfg.get("path", "dynamic_carbon_density.parquet"))
        if not path.is_absolute():
            path = self.config_path.parent / path

        required = ["year", "pft_id"]
        for pool_prefix in ("B", "SS"):
            for cover in ("v", "s", "c", "p"):
                required.append(f"{pool_prefix}_{cover}_tCha")

        frame = pd.read_parquet(path, columns=required, engine="pyarrow")
        missing = [column for column in required if column not in frame.columns]
        if missing:
            raise ValueError(f"Dynamic carbon density missing columns: {missing}")

        pft_ids = set(int(value) for value in frame["pft_id"].unique())
        expected = set(range(1, self.n_pft + 1))
        missing_pfts = sorted(expected - pft_ids)
        if missing_pfts:
            raise ValueError(
                f"Dynamic carbon density is missing PFT IDs: {missing_pfts}"
            )

        self.dynamic_density = frame.set_index(["year", "pft_id"]).sort_index()
        self.dynamic_years = np.asarray(sorted(frame["year"].unique()), dtype=int)
        self.current_density_year = int(self.dynamic_years[0])

    def set_density_year(self, year: int) -> None:
        if not self.use_dynamic_density:
            return
        assert self.dynamic_years is not None
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
            ("Biomass", cover): f"B_{cover}_tCha"
            for cover in ("v", "s", "c", "p")
        }
        column_map.update({
            ("Soil", cover): f"SS_{cover}_tCha"
            for cover in ("v", "s", "c", "p")
        })
        key = (pool_type, land_cover)
        if key not in column_map:
            raise KeyError(f"Unsupported dynamic density key: {key}")
        if self.dynamic_density is None or self.current_density_year is None:
            raise RuntimeError("Dynamic carbon density has not been initialized.")
        return float(
            self.dynamic_density.loc[
                (int(self.current_density_year), int(p) + 1),
                column_map[key],
            ]
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
        return self.config["LUH2toLULC"][self.config["numtoLUH2Type"][pft]]

    def get_LUH2Type(self, pft: int) -> str:
        return self.config["numtoLUH2Type"][pft]

    def get_clearing_param(self, p: int, param_name: str) -> float:
        return float(self.config["clearing_param"][self.get_pft_name(p)][param_name])

    def get_abandonment_param(self, p: int, param_name: str) -> float:
        return float(self.config["abandonment_param"][self.get_pft_name(p)][param_name])

    def get_harvest_param(self, p: int, param_name: str) -> float:
        return float(self.config["harvest_param"][self.get_pft_name(p)][param_name])

    def get_pft_name(self, p: int) -> str:
        if not (0 <= p < self.n_pft):
            raise IndexError(f"p out of range: {p}, n_pft={self.n_pft}")
        return self.pft_names[p]

    def get_n_pft(self) -> int:
        return self.n_pft

    def validate_parameter_completeness(self) -> None:
        """Fail early when a configured PFT lacks parameters used by the model."""
        required_density = {
            "Biomass": ("v", "s", "p", "c"),
            "Soil": ("v", "s", "p", "c"),
        }
        clearing_keys = (
            "Prod_1", "Prod_10", "Prod_100", "Bio_Soil",
            "t_lapse", "t_rest", "f_v", "f_s",
        )
        abandonment_keys = ("t_biomass", "t_soil")
        harvest_keys = (
            "Prod_1", "Prod_10", "Prod_100", "Bio_Soil_v", "Bio_Soil_s",
            "SOC_min_v", "SOC_min_s", "t_lapse", "f_v", "f_s",
        )

        errors: List[str] = []
        for p, name in enumerate(self.pft_names):
            if not self.use_dynamic_density:
                for pool, covers in required_density.items():
                    for cover in covers:
                        try:
                            self._get_static_carbon_density(p, pool, cover)
                        except Exception as exc:  # collect all missing fields
                            errors.append(str(exc))
            for section, keys in (
                ("clearing_param", clearing_keys),
                ("abandonment_param", abandonment_keys),
                ("harvest_param", harvest_keys),
            ):
                values = self.config.get(section, {}).get(name, {})
                for key in keys:
                    if key not in values:
                        errors.append(f"Missing {section}.{name}.{key}")

        if errors:
            preview = "\n".join(f"  - {item}" for item in errors[:30])
            suffix = "" if len(errors) <= 30 else f"\n  ... and {len(errors)-30} more"
            raise ValueError(f"Parameter configuration is incomplete:\n{preview}{suffix}")
