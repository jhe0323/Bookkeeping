"""Monte Carlo parameter sampling utilities for the Bookkeeping model.

This module is intentionally independent of the heavy NetCDF inputs.  It only
operates on YAML-like mappings and produces nested ``parameter_overrides`` that
can be passed directly to :class:`src.parameter_loader.ParameterLoader`.

Python 3.8 compatible.
"""
from __future__ import annotations

from copy import deepcopy
import fnmatch
import json
import math
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import numpy as np


PathTuple = Tuple[str, ...]


def deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> Dict[str, Any]:
    """Recursively merge two mappings without mutating either input."""
    merged = deepcopy(dict(base))
    for key, value in override.items():
        if (
            key in merged
            and isinstance(merged[key], dict)
            and isinstance(value, Mapping)
        ):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def path_to_string(path: Sequence[str]) -> str:
    return ".".join(path)


def _split_pattern(pattern: str) -> PathTuple:
    parts = tuple(part.strip() for part in str(pattern).split(".") if part.strip())
    if not parts:
        raise ValueError("Parameter target path cannot be empty.")
    return parts


def expand_target_pattern(config: Mapping[str, Any], pattern: str) -> List[PathTuple]:
    """Expand a dotted target path that may contain shell-style wildcards.

    Example
    -------
    ``carbon_density.PFT*.Biomass.v``
    """
    parts = _split_pattern(pattern)
    out: List[PathTuple] = []

    def walk(node: Any, index: int, prefix: Tuple[str, ...]) -> None:
        if index == len(parts):
            if isinstance(node, Mapping):
                raise ValueError(
                    "Target resolves to a mapping, not a scalar: {}".format(
                        path_to_string(prefix)
                    )
                )
            out.append(prefix)
            return
        if not isinstance(node, Mapping):
            return
        token = parts[index]
        matches = [
            str(key)
            for key in node.keys()
            if fnmatch.fnmatchcase(str(key), token)
        ]
        for key in sorted(matches):
            walk(node[key], index + 1, prefix + (key,))

    walk(config, 0, tuple())
    if not out:
        raise KeyError("Target pattern matched nothing: {}".format(pattern))
    return out


def get_path(config: Mapping[str, Any], path: Sequence[str]) -> Any:
    node: Any = config
    for key in path:
        if not isinstance(node, Mapping) or key not in node:
            raise KeyError("Missing parameter path: {}".format(path_to_string(path)))
        node = node[key]
    return node


def set_path(config: MutableMapping[str, Any], path: Sequence[str], value: Any) -> None:
    if not path:
        raise ValueError("Cannot set an empty parameter path.")
    node: MutableMapping[str, Any] = config
    for key in path[:-1]:
        child = node.get(key)
        if child is None:
            child = {}
            node[key] = child
        if not isinstance(child, MutableMapping):
            raise TypeError(
                "Cannot descend through non-mapping parameter path: {}".format(
                    path_to_string(path)
                )
            )
        node = child
    node[path[-1]] = value


def _clip(value: float, lower: Any = None, upper: Any = None) -> float:
    out = float(value)
    if lower is not None:
        out = max(float(lower), out)
    if upper is not None:
        out = min(float(upper), out)
    return out


def _draw_multiplier(rng: np.random.Generator, spec: Mapping[str, Any]) -> float:
    dist = str(spec.get("distribution", "uniform_multiplier")).lower()

    if dist == "uniform_multiplier":
        low = float(spec.get("low", 0.9))
        high = float(spec.get("high", 1.1))
        if high < low:
            raise ValueError("uniform_multiplier requires high >= low")
        return float(rng.uniform(low, high))

    if dist == "normal_multiplier":
        mean = float(spec.get("mean", 1.0))
        sd = float(spec["sd"])
        return _clip(
            rng.normal(mean, sd),
            spec.get("lower"),
            spec.get("upper"),
        )

    if dist == "lognormal_multiplier":
        cv = float(spec["cv"])
        if cv < 0.0:
            raise ValueError("lognormal_multiplier cv must be >= 0")
        sigma2 = math.log1p(cv * cv)
        sigma = math.sqrt(sigma2)
        mu = -0.5 * sigma2  # E[multiplier] = 1
        return float(rng.lognormal(mean=mu, sigma=sigma))

    raise ValueError("Unsupported multiplier distribution: {}".format(dist))


def _draw_direct(
    rng: np.random.Generator,
    spec: Mapping[str, Any],
    baseline: float,
) -> float:
    dist = str(spec.get("distribution", "normal")).lower()

    if dist == "uniform":
        low = float(spec["low"])
        high = float(spec["high"])
        if high < low:
            raise ValueError("uniform requires high >= low")
        value = rng.uniform(low, high)
    elif dist == "normal":
        value = rng.normal(float(spec.get("mean", baseline)), float(spec["sd"]))
    elif dist == "lognormal":
        mean = float(spec.get("mean", baseline))
        cv = float(spec["cv"])
        if mean <= 0.0:
            raise ValueError("lognormal direct sampling requires positive mean")
        sigma2 = math.log1p(cv * cv)
        sigma = math.sqrt(sigma2)
        mu = math.log(mean) - 0.5 * sigma2
        value = rng.lognormal(mu, sigma)
    elif dist == "beta_from_baseline":
        concentration = float(spec.get("concentration", 50.0))
        if concentration <= 0.0:
            raise ValueError("beta_from_baseline concentration must be > 0")
        mean = float(baseline)
        if mean <= 0.0 or mean >= 1.0:
            # Preserve exact structural zeros/ones instead of creating tiny
            # artificial fractions where the reference configuration has none.
            value = mean
        else:
            alpha = mean * concentration
            beta = (1.0 - mean) * concentration
            value = rng.beta(alpha, beta)
    else:
        raise ValueError("Unsupported direct distribution: {}".format(dist))

    return _clip(value, spec.get("lower"), spec.get("upper"))


def _coerce_sampled_value(value: float, spec: Mapping[str, Any]) -> Any:
    if bool(spec.get("integer", False)):
        return int(round(float(value)))
    return float(value)


def _apply_scalar_rule(
    *,
    baseline_config: Mapping[str, Any],
    overrides: MutableMapping[str, Any],
    rng: np.random.Generator,
    rule: Mapping[str, Any],
    draws: MutableMapping[str, Any],
) -> None:
    name = str(rule.get("name", "unnamed_parameter"))
    targets_raw = rule.get("targets")
    if isinstance(targets_raw, str):
        targets = [targets_raw]
    else:
        targets = list(targets_raw or [])
    if not targets:
        raise ValueError("Scalar rule {!r} has no targets".format(name))

    expanded: List[PathTuple] = []
    for pattern in targets:
        expanded.extend(expand_target_pattern(baseline_config, str(pattern)))
    # Deduplicate while preserving deterministic order.
    expanded = sorted(set(expanded))

    scope = str(rule.get("scope", "independent")).lower()
    if scope not in {"independent", "shared"}:
        raise ValueError("Rule {!r}: scope must be independent or shared".format(name))

    distribution = str(rule.get("distribution", "uniform_multiplier")).lower()
    multiplier_mode = distribution.endswith("_multiplier")
    shared_draw = None
    if scope == "shared":
        if multiplier_mode:
            shared_draw = _draw_multiplier(rng, rule)
        else:
            # A shared direct draw is only sensible when all targets are on the
            # same numerical scale; support it deliberately but keep it explicit.
            first_baseline = float(get_path(baseline_config, expanded[0]))
            shared_draw = _draw_direct(rng, rule, first_baseline)

    rule_draws: Dict[str, Any] = {}
    for path in expanded:
        baseline = float(get_path(baseline_config, path))
        if multiplier_mode:
            multiplier = float(shared_draw) if shared_draw is not None else _draw_multiplier(rng, rule)
            sampled = baseline * multiplier
            sampled = _clip(sampled, rule.get("value_lower"), rule.get("value_upper"))
            draw_value: Any = multiplier
        else:
            sampled = float(shared_draw) if shared_draw is not None else _draw_direct(rng, rule, baseline)
            draw_value = sampled

        sampled = _coerce_sampled_value(sampled, rule)
        set_path(overrides, path, sampled)
        rule_draws[path_to_string(path)] = draw_value

    draws[name] = rule_draws


def _expand_parent_pattern(config: Mapping[str, Any], pattern: str) -> List[PathTuple]:
    """Expand a dotted wildcard pattern expected to resolve to mappings."""
    parts = _split_pattern(pattern)
    out: List[PathTuple] = []

    def walk(node: Any, index: int, prefix: Tuple[str, ...]) -> None:
        if index == len(parts):
            if not isinstance(node, Mapping):
                raise ValueError(
                    "Simplex parent must resolve to a mapping: {}".format(
                        path_to_string(prefix)
                    )
                )
            out.append(prefix)
            return
        if not isinstance(node, Mapping):
            return
        token = parts[index]
        for key in sorted(str(k) for k in node.keys() if fnmatch.fnmatchcase(str(k), token)):
            walk(node[key], index + 1, prefix + (key,))

    walk(config, 0, tuple())
    if not out:
        raise KeyError("Simplex parent pattern matched nothing: {}".format(pattern))
    return out


def _apply_simplex_group(
    *,
    baseline_config: Mapping[str, Any],
    overrides: MutableMapping[str, Any],
    rng: np.random.Generator,
    group: Mapping[str, Any],
    draws: MutableMapping[str, Any],
) -> None:
    name = str(group.get("name", "unnamed_simplex"))
    parent_pattern = str(group["parent_pattern"])
    fields = [str(value) for value in group.get("fields", [])]
    if len(fields) < 2:
        raise ValueError("Simplex group {!r} requires at least two fields".format(name))
    concentration = float(group.get("concentration", 100.0))
    if concentration <= 0.0:
        raise ValueError("Simplex group concentration must be > 0")
    preserve_total = bool(group.get("preserve_total", True))

    group_draws: Dict[str, Any] = {}
    for parent in _expand_parent_pattern(baseline_config, parent_pattern):
        values = np.asarray(
            [float(get_path(baseline_config, parent + (field,))) for field in fields],
            dtype=float,
        )
        if np.any(values < 0.0):
            raise ValueError(
                "Simplex baseline contains a negative fraction at {}".format(
                    path_to_string(parent)
                )
            )
        total = float(values.sum())
        if total <= 0.0:
            # Structural all-zero group: keep it all zero.
            sampled = values.copy()
        else:
            positive = values > 0.0
            if int(positive.sum()) <= 1:
                sampled = values.copy()
            else:
                proportions = values[positive] / total
                alpha = np.maximum(proportions * concentration, 1.0e-9)
                sampled_positive = rng.dirichlet(alpha)
                sampled = np.zeros_like(values)
                sampled[positive] = sampled_positive
                if preserve_total:
                    sampled *= total
        for field, value in zip(fields, sampled):
            set_path(overrides, parent + (field,), float(value))
        group_draws[path_to_string(parent)] = {
            field: float(value) for field, value in zip(fields, sampled)
        }
    draws[name] = group_draws


def validate_model_parameters(config: Mapping[str, Any], tol: float = 1.0e-10) -> None:
    """Validate MC-sensitive physical constraints in the effective parameter set."""
    errors: List[str] = []

    density = config.get("carbon_density", {}) or {}
    for pft, pft_cfg in density.items():
        if not isinstance(pft_cfg, Mapping):
            continue
        for pool_name in ("Biomass", "Soil"):
            pool = pft_cfg.get(pool_name, {}) or {}
            if isinstance(pool, Mapping):
                for cover, value in pool.items():
                    try:
                        number = float(value)
                    except Exception:
                        errors.append("{}.{}.{} is not numeric".format(pft, pool_name, cover))
                        continue
                    if number < -tol:
                        errors.append("{}.{}.{} < 0".format(pft, pool_name, cover))

    clearing = config.get("clearing_param", {}) or {}
    for pft, cfg in clearing.items():
        if not isinstance(cfg, Mapping):
            continue
        fraction_fields = ("Prod_1", "Prod_10", "Prod_100", "Bio_Soil", "f_v", "f_s")
        for field in fraction_fields:
            if field in cfg:
                value = float(cfg[field])
                if value < -tol or value > 1.0 + tol:
                    errors.append("clearing_param.{}.{} outside [0,1]".format(pft, field))
        alloc_fields = ("Prod_1", "Prod_10", "Prod_100", "Bio_Soil")
        if all(field in cfg for field in alloc_fields):
            total = sum(float(cfg[field]) for field in alloc_fields)
            if total > 1.0 + tol:
                errors.append(
                    "clearing_param.{} allocation sum {:.8f} > 1".format(pft, total)
                )
        for field in ("t_lapse", "t_rest"):
            if field in cfg and float(cfg[field]) < -tol:
                errors.append("clearing_param.{}.{} < 0".format(pft, field))

    abandonment = config.get("abandonment_param", {}) or {}
    for pft, cfg in abandonment.items():
        if not isinstance(cfg, Mapping):
            continue
        for field in ("t_biomass", "t_soil"):
            if field in cfg and float(cfg[field]) < -tol:
                errors.append("abandonment_param.{}.{} < 0".format(pft, field))

    harvest = config.get("harvest_param", {}) or {}
    for pft, cfg in harvest.items():
        if not isinstance(cfg, Mapping):
            continue
        for field in ("Prod_1", "Prod_10", "Prod_100", "Bio_Soil_v", "Bio_Soil_s", "f_v", "f_s"):
            if field in cfg:
                value = float(cfg[field])
                if value < -tol or value > 1.0 + tol:
                    errors.append("harvest_param.{}.{} outside [0,1]".format(pft, field))
        product_fields = ("Prod_1", "Prod_10", "Prod_100")
        if all(field in cfg for field in product_fields):
            total = sum(float(cfg[field]) for field in product_fields)
            if total > 1.0 + tol:
                errors.append(
                    "harvest_param.{} product fraction sum {:.8f} > 1".format(pft, total)
                )
        for field in ("SOC_min_v", "SOC_min_s", "t_lapse"):
            if field in cfg and float(cfg[field]) < -tol:
                errors.append("harvest_param.{}.{} < 0".format(pft, field))

    if errors:
        preview = "\n  - " + "\n  - ".join(errors[:30])
        if len(errors) > 30:
            preview += "\n  - ... {} more".format(len(errors) - 30)
        raise ValueError("Invalid Monte Carlo parameter set:" + preview)


def sample_overrides(
    *,
    baseline_config: Mapping[str, Any],
    rng: np.random.Generator,
    spec: Mapping[str, Any],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Draw one parameter realization and return ``(overrides, draws)``."""
    overrides: Dict[str, Any] = {}
    draws: Dict[str, Any] = {}

    for rule in spec.get("parameters", []) or []:
        if not bool(rule.get("enabled", True)):
            continue
        _apply_scalar_rule(
            baseline_config=baseline_config,
            overrides=overrides,
            rng=rng,
            rule=rule,
            draws=draws,
        )

    for group in spec.get("simplex_groups", []) or []:
        if not bool(group.get("enabled", True)):
            continue
        _apply_simplex_group(
            baseline_config=baseline_config,
            overrides=overrides,
            rng=rng,
            group=group,
            draws=draws,
        )

    effective = deep_merge(baseline_config, overrides)
    validate_model_parameters(effective)
    return overrides, draws


def canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
