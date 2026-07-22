"""Configuration-file based run management."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import os
import uuid
from typing import Any, Dict, Optional, Union

import yaml


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


def _resolve(repo_root: Path, value: Union[str, Path]) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = repo_root / path
    return path.resolve()


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
    """Return a fast fingerprint for a potentially large input file.

    The fingerprint combines file size, nanosecond modification time, and a
    SHA-256 digest of the first and last sample blocks. It is intentionally
    much cheaper than hashing an entire multi-gigabyte NetCDF file while still
    detecting normal same-name input replacements.
    """

    resolved = path.resolve()
    stat = resolved.stat()
    size_bytes = int(stat.st_size)
    sample_bytes = max(1, int(sample_bytes))

    digest = hashlib.sha256()
    digest.update(str(size_bytes).encode("ascii"))
    with resolved.open("rb") as handle:
        first = handle.read(sample_bytes)
        digest.update(first)
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


def build_run_metadata(config: ResolvedRunConfig) -> Dict[str, Any]:
    return {
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
        "parameter_file": str(config.base_parameter_path),
        "experiment_file": str(config.experiment_path),
        "run_config_file": str(config.config_path),
        "run_config_sha256": hashlib.sha256(config.raw_text.encode("utf-8")).hexdigest(),
        "parameter_sha256": sha256_file(config.base_parameter_path),
        "experiment_sha256": sha256_file(config.experiment_path),
        "git_commit": git_commit(config.repo_root),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_config_yaml": config.raw_text,
    }


def _atomic_write_text(path: Path, text: str) -> None:
    """Atomically replace a UTF-8 text file."""

    temporary = path.with_name(
        ".{}.{}.{}.tmp".format(path.name, os.getpid(), uuid.uuid4().hex)
    )
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(str(temporary), str(path))
    finally:
        temporary.unlink(missing_ok=True)


def manifest_identity(metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Fields that uniquely identify one model configuration."""

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
        "parameter_file",
        "experiment_file",
        "run_config_sha256",
        "parameter_sha256",
        "experiment_sha256",
        "git_commit",
    )
    return {
        key: metadata.get(key)
        for key in keys
    }


def write_run_manifest(
    config: ResolvedRunConfig,
    output_dir: Path,
    metadata: Dict[str, Any],
    *,
    manifest_name: str = "run_manifest.json",
    overwrite_manifest: bool = True,
) -> None:
    """Write the copied run YAML and one run or band manifest.

    A run directory cannot be reused with a different run YAML. This prevents
    output bands produced from different model configurations from being mixed.
    """

    output_dir.mkdir(parents=True, exist_ok=True)

    config_copy = output_dir / "run_config.yml"

    try:
        # Exclusive creation is safe when multiple SLURM tasks start together.
        with config_copy.open("x", encoding="utf-8") as handle:
            handle.write(config.raw_text)
    except FileExistsError:
        existing_config = config_copy.read_text(encoding="utf-8")
        if existing_config != config.raw_text:
            raise RuntimeError(
                "The output directory already contains a different "
                "run_config.yml. Use a new run.name or remove the old "
                "run directory before starting this simulation."
            )

    manifest_path = output_dir / manifest_name

    if manifest_path.exists() and not overwrite_manifest:
        try:
            existing_metadata = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
        except Exception as exc:
            raise RuntimeError(
                "Existing manifest cannot be read: {}".format(manifest_path)
            ) from exc

        if manifest_identity(existing_metadata) != manifest_identity(metadata):
            raise RuntimeError(
                "The output directory belongs to a different model build or "
                "configuration. Use a new run.name. Existing manifest: {}"
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
