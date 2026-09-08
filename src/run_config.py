"""Configuration-file based run management."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import os
import uuid
from typing import Any, Dict, Optional, Union

import yaml


# Files whose contents can change deterministic numerical results.
# Deliberately excludes README, scheduler scripts, MC-only code, and run
# orchestration/manifest code so non-numerical maintenance does not invalidate
# already-computed deterministic outputs.
_MODEL_RESULT_FILES = (
    "src/LULCCSimulator.py",
    "src/parameter_loader.py",
    "src/carbon_pools_init.py",
    "src/events.py",
    "src/harvest.py",
    "src/transition.py",
    "src/summary_yearly.py",
    "src/file_loader.py",
)


# Numerical setup logic that lives in orchestration/configuration modules.
# Only these AST nodes/calls are hashed, so edits to logs/manifests/comments do
# not invalidate numerical results, while changes to config interpretation,
# band slicing, simulator construction, or the grid-run call do.
_NUMERICAL_SETUP_AST_NODES = {
    "src/run_config.py": ("ResolvedRunConfig", "load_run_config"),
    "src/run_manager.py": ("_band_slice",),
}
_NUMERICAL_SETUP_CALLS = {
    "src/run_manager.py": (
        ("run_from_config", "LULCCSimulator"),
        ("run_from_config", "run_simulation_grid"),
    ),
}


@dataclass
class ResolvedRunConfig:
    config_path: Path
    repo_root: Path
    raw: Dict[str, Any]
    raw_text: str
    run_name: str
    resolution: str
    input_format: str
    input_base_year: int
    start_year: int
    end_year: int
    time_encoding: str
    area_unit: str
    state_path: Path
    transition_path: Path
    pft_path: Path
    pft_variable: Optional[str]
    pft_base_year: int
    pft_time_encoding: str
    pft_update_mode: str
    pft_min_year_policy: str
    pft_max_year_policy: str
    pft_missing_cell_policy: str
    pft_nearest_search_radius: int
    pft_default_index: Optional[int]
    experiment_path: Path
    base_parameter_path: Path
    output_root: Path
    output_prefix: str
    sync_every: int
    compression_level: int
    validate_before_run: bool
    validation_full_scan: bool
    band_size_deg: float
    stagger_max: float

    @property
    def years(self) -> int:
        return self.end_year - self.start_year

    @property
    def calendar_years(self):
        return range(self.start_year, self.end_year + 1)


def _load_yaml(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Top-level YAML must be a mapping: {path}")
    return data


def _load_yaml_text(text: str) -> Dict[str, Any]:
    data = yaml.safe_load(text) or {}
    if not isinstance(data, dict):
        raise ValueError("Top-level YAML text must be a mapping")
    return data


def _resolve(repo_root: Path, value: Union[str, Path]) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = repo_root / path
    return path.resolve()


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )


def semantic_yaml_sha256(text: str) -> str:
    """Hash YAML meaning rather than comments/spacing/key order."""
    parsed = _load_yaml_text(text)
    return hashlib.sha256(_canonical_json(parsed).encode("utf-8")).hexdigest()


def load_run_config(
    path: Union[str, Path],
    *,
    expected_resolution: Optional[str] = None,
) -> ResolvedRunConfig:
    config_path = Path(path).resolve()
    raw_text = config_path.read_text(encoding="utf-8")
    raw = _load_yaml(config_path)
    repo_root = config_path.parent.parent.resolve()

    run = raw.get("run", {}) or {}
    paths = raw.get("paths", {}) or {}
    inputs = raw.get("inputs", {}) or {}
    output = raw.get("output", {}) or {}
    validation = raw.get("validation", {}) or {}
    server = raw.get("server", {}) or {}
    pft = raw.get("pft", {}) or {}

    resolution = str(run.get("resolution", "025deg"))
    if expected_resolution is not None and resolution != expected_resolution:
        raise ValueError(
            f"Run config resolution={resolution!r}, expected {expected_resolution!r}."
        )

    input_format = str(run.get("input_format", "original_luh2")).lower()
    if input_format not in {"original_luh2", "vscp"}:
        raise ValueError("run.input_format must be original_luh2 or vscp")

    start_year = int(run.get("start_year", 850))
    end_year = int(run.get("end_year", 2020))
    if end_year <= start_year:
        raise ValueError("run.end_year must be greater than run.start_year")

    data_dir = _resolve(repo_root, paths.get("data_dir", "../In_ncfile"))
    output_root = _resolve(repo_root, paths.get("output_dir", "../Out_ncfile"))

    state_path = _resolve(data_dir, inputs.get("states", "states.nc"))
    transition_path = _resolve(data_dir, inputs.get("transitions", "transitions.nc"))
    pft_path = _resolve(data_dir, inputs.get("pft", "IBIS_PFT_dominant_0.25deg.nc"))

    experiment_path = _resolve(
        repo_root,
        raw.get("experiment", {}).get("path", "config/experiments/harvest_area.yml"),
    )
    base_parameter_path = _resolve(
        repo_root,
        raw.get("parameters", {}).get("path", "config/config.yml"),
    )

    resolved = ResolvedRunConfig(
        config_path=config_path,
        repo_root=repo_root,
        raw=deepcopy(raw),
        raw_text=raw_text,
        run_name=str(run.get("name", f"{resolution}_{input_format}")),
        resolution=resolution,
        input_format=input_format,
        input_base_year=int(run.get("input_base_year", 850)),
        start_year=start_year,
        end_year=end_year,
        time_encoding=str(inputs.get("time_encoding", "index")),
        area_unit=str(run.get("area_unit", "ha")),
        state_path=state_path,
        transition_path=transition_path,
        pft_path=pft_path,
        pft_variable=inputs.get("pft_variable"),
        pft_base_year=int(pft.get("base_year", run.get("input_base_year", 850))),
        pft_time_encoding=str(pft.get("time_encoding", "auto")),
        pft_update_mode=str(pft.get("update_mode", "annual_conservative")),
        pft_min_year_policy=str(pft.get("min_year_policy", "clip")),
        pft_max_year_policy=str(pft.get("max_year_policy", "clip")),
        pft_missing_cell_policy=str(pft.get("missing_cell_policy", "error")),
        pft_nearest_search_radius=int(pft.get("nearest_search_radius", 8)),
        pft_default_index=(
            None if pft.get("default_index") is None
            else int(pft.get("default_index"))
        ),
        experiment_path=experiment_path,
        base_parameter_path=base_parameter_path,
        output_root=output_root,
        output_prefix=str(output.get("filename_prefix", f"summary_{resolution}")),
        sync_every=int(output.get("sync_every", 200)),
        compression_level=int(output.get("compression_level", 4)),
        validate_before_run=bool(validation.get("run_before_simulation", False)),
        validation_full_scan=bool(validation.get("full_scan", True)),
        band_size_deg=float(server.get("band_size_deg", 3 if resolution == "025deg" else 10)),
        stagger_max=float(server.get("stagger_max", 0)),
    )

    for required in (
        resolved.state_path,
        resolved.transition_path,
        resolved.pft_path,
        resolved.experiment_path,
        resolved.base_parameter_path,
    ):
        if not required.is_file():
            raise FileNotFoundError(f"Configured input file not found: {required}")
    return resolved


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_fingerprint(path: Path, sample_bytes: int = 1024 * 1024) -> Dict[str, Any]:
    """Return a fast fingerprint for a potentially large input file."""
    resolved = path.resolve()
    stat = resolved.stat()
    size_bytes = int(stat.st_size)
    sample_bytes = max(1, int(sample_bytes))

    digest = hashlib.sha256()
    digest.update(str(size_bytes).encode("ascii"))
    with resolved.open("rb") as handle:
        digest.update(handle.read(sample_bytes))
        if size_bytes > sample_bytes:
            handle.seek(max(0, size_bytes - sample_bytes))
            digest.update(handle.read(sample_bytes))

    return {
        "path": str(resolved),
        "size_bytes": size_bytes,
        "mtime_ns": int(stat.st_mtime_ns),
        "sample_bytes": sample_bytes,
        "sample_sha256": digest.hexdigest(),
    }


def git_commit(repo_root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def model_code_sha256(repo_root: Path, ref: Optional[str] = None) -> str:
    """Hash result-affecting deterministic numerical source code.

    If ``ref`` is supplied, hash the files as stored at that Git commit. This is
    used only to migrate/validate legacy manifests that predate direct code hashes.
    """
    digest = hashlib.sha256()
    found = 0
    for rel in _MODEL_RESULT_FILES:
        if ref is None:
            path = repo_root / rel
            if not path.is_file():
                raise FileNotFoundError(f"Model source file not found: {path}")
            file_sha = sha256_file(path)
        else:
            try:
                content = subprocess.check_output(
                    ["git", "-C", str(repo_root), "show", f"{ref}:{rel}"],
                    stderr=subprocess.DEVNULL,
                )
            except Exception as exc:
                raise RuntimeError(
                    f"Cannot read {rel} at git ref {ref!r}"
                ) from exc
            file_sha = hashlib.sha256(content).hexdigest()
        found += 1
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_sha.encode("ascii"))
        digest.update(b"\0")
    if found == 0:
        raise RuntimeError(f"Could not hash model code under {repo_root}")
    return digest.hexdigest()



def _source_text(repo_root: Path, rel: str, ref: Optional[str] = None) -> str:
    """Read one repository source file from the working tree or a Git ref."""
    if ref is None:
        path = repo_root / rel
        if not path.is_file():
            raise FileNotFoundError(f"Model source file not found: {path}")
        return path.read_text(encoding="utf-8")
    try:
        content = subprocess.check_output(
            ["git", "-C", str(repo_root), "show", f"{ref}:{rel}"],
            stderr=subprocess.DEVNULL,
        )
    except Exception as exc:
        raise RuntimeError(f"Cannot read {rel} at git ref {ref!r}") from exc
    return content.decode("utf-8")


def _find_top_level_ast_node(tree: ast.AST, name: str) -> ast.AST:
    for node in getattr(tree, "body", []):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.name == name:
                return node
    raise RuntimeError(f"Could not find AST node {name!r}")


def _find_function_call_ast(tree: ast.AST, function_name: str, call_name: str) -> ast.Call:
    function = _find_top_level_ast_node(tree, function_name)
    for node in ast.walk(function):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if isinstance(func, ast.Name) and func.id == call_name:
            return node
        if isinstance(func, ast.Attribute) and func.attr == call_name:
            return node
    raise RuntimeError(
        f"Could not find call {call_name!r} inside function {function_name!r}"
    )


def numerical_setup_sha256(repo_root: Path, ref: Optional[str] = None) -> str:
    """Hash result-affecting setup logic without hashing provenance plumbing.

    This supplements ``model_code_sha256``. It catches changes in:
    - run-YAML interpretation and resolved numerical fields;
    - longitude-band slicing;
    - arguments passed to ``LULCCSimulator``;
    - the grid simulation call.

    When ``ref`` is supplied, the same AST fragments are reconstructed from the
    historical Git commit. This allows existing manifests to be migrated
    conservatively without forcing a rerun after provenance-only maintenance.
    """
    digest = hashlib.sha256()
    parsed: Dict[str, ast.AST] = {}

    files = sorted(
        set(_NUMERICAL_SETUP_AST_NODES) | set(_NUMERICAL_SETUP_CALLS)
    )
    for rel in files:
        source = _source_text(repo_root, rel, ref=ref)
        parsed[rel] = ast.parse(source, filename=rel)

    for rel in sorted(_NUMERICAL_SETUP_AST_NODES):
        tree = parsed[rel]
        for node_name in _NUMERICAL_SETUP_AST_NODES[rel]:
            node = _find_top_level_ast_node(tree, node_name)
            payload = ast.dump(node, annotate_fields=True, include_attributes=False)
            digest.update(rel.encode("utf-8"))
            digest.update(b"\0node\0")
            digest.update(node_name.encode("utf-8"))
            digest.update(b"\0")
            digest.update(payload.encode("utf-8"))
            digest.update(b"\0")

    for rel in sorted(_NUMERICAL_SETUP_CALLS):
        tree = parsed[rel]
        for function_name, call_name in _NUMERICAL_SETUP_CALLS[rel]:
            node = _find_function_call_ast(tree, function_name, call_name)
            payload = ast.dump(node, annotate_fields=True, include_attributes=False)
            digest.update(rel.encode("utf-8"))
            digest.update(b"\0call\0")
            digest.update(function_name.encode("utf-8"))
            digest.update(b"\0")
            digest.update(call_name.encode("utf-8"))
            digest.update(b"\0")
            digest.update(payload.encode("utf-8"))
            digest.update(b"\0")

    return digest.hexdigest()

def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    merged = deepcopy(base)
    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _effective_parameter_config(config: ResolvedRunConfig) -> Dict[str, Any]:
    effective = _load_yaml(config.base_parameter_path)
    if config.experiment_path.is_file():
        effective = _deep_merge(effective, _load_yaml(config.experiment_path))
    effective = _deep_merge(effective, config.raw.get("model_overrides", {}) or {})
    return effective


def _dynamic_density_metadata(config: ResolvedRunConfig) -> Dict[str, Any]:
    effective = _effective_parameter_config(config)
    dynamic = effective.get("dynamic_carbon_density", {}) or {}
    if not bool(dynamic.get("enabled", False)):
        return {
            "dynamic_density_file": "disabled",
            "dynamic_density_fingerprint": "disabled",
        }

    density_path = Path(dynamic.get("path", "dynamic_carbon_density.parquet"))
    if not density_path.is_absolute():
        density_path = config.base_parameter_path.parent / density_path
    density_path = density_path.resolve()
    if not density_path.is_file():
        raise FileNotFoundError(
            f"Dynamic carbon-density file not found: {density_path}"
        )
    return {
        "dynamic_density_file": str(density_path),
        "dynamic_density_fingerprint": file_fingerprint(density_path),
    }


def build_run_metadata(config: ResolvedRunConfig) -> Dict[str, Any]:
    metadata = {
        "run_name": config.run_name,
        "resolution": config.resolution,
        "input_format": config.input_format,
        "input_base_year": config.input_base_year,
        "start_year": config.start_year,
        "end_year": config.end_year,
        "area_unit": config.area_unit,
        "state_file": str(config.state_path),
        "transition_file": str(config.transition_path),
        "pft_file": str(config.pft_path),
        "state_fingerprint": file_fingerprint(config.state_path),
        "transition_fingerprint": file_fingerprint(config.transition_path),
        "pft_fingerprint": file_fingerprint(config.pft_path),
        "pft_variable_requested": config.pft_variable or "auto",
        "pft_update_mode": config.pft_update_mode,
        "pft_missing_cell_policy": config.pft_missing_cell_policy,
        "pft_nearest_search_radius": config.pft_nearest_search_radius,
        "pft_default_index": (
            config.pft_default_index if config.pft_default_index is not None else "none"
        ),
        "parameter_file": str(config.base_parameter_path),
        "experiment_file": str(config.experiment_path),
        "run_config_file": str(config.config_path),
        # Semantic hash: YAML comments/formatting do not change result identity.
        "run_config_sha256": semantic_yaml_sha256(config.raw_text),
        "parameter_sha256": sha256_file(config.base_parameter_path),
        "experiment_sha256": sha256_file(config.experiment_path),
        "model_code_sha256": model_code_sha256(config.repo_root),
        "numerical_setup_sha256": numerical_setup_sha256(config.repo_root),
        # Git commit remains provenance only; it is intentionally not an identity key.
        "git_commit": git_commit(config.repo_root),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_config_yaml": config.raw_text,
    }
    metadata.update(_dynamic_density_metadata(config))
    return metadata


def _atomic_write_text(path: Path, text: str) -> None:
    temporary = path.with_name(
        ".{}.{}.{}.tmp".format(path.name, os.getpid(), uuid.uuid4().hex)
    )
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(str(temporary), str(path))
    finally:
        temporary.unlink(missing_ok=True)


def manifest_identity(metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Fields that uniquely identify one deterministic numerical result."""
    keys = (
        "run_name",
        "resolution",
        "input_format",
        "input_base_year",
        "start_year",
        "end_year",
        "state_file",
        "transition_file",
        "pft_file",
        "state_fingerprint",
        "transition_fingerprint",
        "pft_fingerprint",
        "pft_update_mode",
        "pft_missing_cell_policy",
        "pft_nearest_search_radius",
        "pft_default_index",
        "parameter_file",
        "experiment_file",
        "run_config_sha256",
        "parameter_sha256",
        "experiment_sha256",
        "model_code_sha256",
        "numerical_setup_sha256",
        "dynamic_density_file",
        "dynamic_density_fingerprint",
    )
    return {key: metadata.get(key) for key in keys}


def _normalize_metadata_value(value: Any) -> str:
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.startswith("{") or stripped.startswith("["):
            try:
                value = json.loads(stripped)
            except json.JSONDecodeError:
                pass
    return _canonical_json(value)


def metadata_identity_matches(
    existing: Dict[str, Any],
    expected: Dict[str, Any],
    *,
    repo_root: Optional[Path] = None,
) -> bool:
    """Compare identities with compatibility for pre-code-hash manifests.

    Compatibility rules are intentionally conservative:
    - legacy raw-text run-config hashes may be upgraded using stored run_config_yaml;
    - a missing model-code hash may be reconstructed from the legacy git_commit;
    - missing dynamic-density metadata is accepted only when dynamic density is disabled.
    """
    expected_identity = manifest_identity(expected)

    for key, expected_value in expected_identity.items():
        existing_value = existing.get(key)

        if key == "run_config_sha256":
            if _normalize_metadata_value(existing_value) == _normalize_metadata_value(expected_value):
                continue
            legacy_yaml = existing.get("run_config_yaml")
            if legacy_yaml is not None:
                try:
                    if semantic_yaml_sha256(str(legacy_yaml)) == str(expected_value):
                        continue
                except Exception:
                    pass
            return False

        if key == "model_code_sha256" and existing_value is None:
            if repo_root is None:
                return False
            legacy_commit = str(existing.get("git_commit", "")).strip()
            if not legacy_commit or legacy_commit == "unknown":
                return False
            try:
                legacy_hash = model_code_sha256(repo_root, ref=legacy_commit)
            except Exception:
                return False
            if legacy_hash == str(expected_value):
                continue
            return False

        if key == "numerical_setup_sha256" and existing_value is None:
            if repo_root is None:
                return False
            legacy_commit = str(existing.get("git_commit", "")).strip()
            if not legacy_commit or legacy_commit == "unknown":
                return False
            try:
                legacy_hash = numerical_setup_sha256(repo_root, ref=legacy_commit)
            except Exception:
                return False
            if legacy_hash == str(expected_value):
                continue
            return False

        if key in {"dynamic_density_file", "dynamic_density_fingerprint"}:
            if existing_value is None and expected_value == "disabled":
                continue

        if _normalize_metadata_value(existing_value) != _normalize_metadata_value(expected_value):
            return False

    return True


def write_run_manifest(
    config: ResolvedRunConfig,
    output_dir: Path,
    metadata: Dict[str, Any],
    *,
    manifest_name: str = "run_manifest.json",
    overwrite_manifest: bool = True,
) -> None:
    """Write copied run YAML and one run/band manifest safely."""
    output_dir.mkdir(parents=True, exist_ok=True)
    config_copy = output_dir / "run_config.yml"

    try:
        with config_copy.open("x", encoding="utf-8") as handle:
            handle.write(config.raw_text)
    except FileExistsError:
        existing_config = config_copy.read_text(encoding="utf-8")
        if existing_config != config.raw_text:
            try:
                same_semantics = (
                    semantic_yaml_sha256(existing_config)
                    == semantic_yaml_sha256(config.raw_text)
                )
            except Exception:
                same_semantics = False
            if not same_semantics:
                raise RuntimeError(
                    "The output directory already contains a semantically different "
                    "run_config.yml. Use a new run.name or remove the old run directory."
                )

    manifest_path = output_dir / manifest_name

    if manifest_path.exists() and not overwrite_manifest:
        try:
            existing_metadata = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
        except Exception as exc:
            raise RuntimeError(
                f"Existing manifest cannot be read: {manifest_path}"
            ) from exc

        if not metadata_identity_matches(
            existing_metadata,
            metadata,
            repo_root=config.repo_root,
        ):
            raise RuntimeError(
                "The output directory belongs to a different numerical model build "
                "or configuration. Use a new run.name. Existing manifest: {}"
                .format(manifest_path)
            )
        return

    manifest_text = json.dumps(
        metadata,
        indent=2,
        ensure_ascii=False,
        sort_keys=True,
    )
    _atomic_write_text(manifest_path, manifest_text)
