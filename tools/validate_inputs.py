#!/usr/bin/env python3
"""Run the one-time input validator.

Edit CONFIG_PATH below when validating another run configuration.
"""
from pathlib import Path
import json

from src.input_validator import validate_run_config_file

CONFIG_PATH = Path("config/run_025deg.yml")


if __name__ == "__main__":
    report = validate_run_config_file(CONFIG_PATH)
    print(json.dumps({
        "status": report["status"],
        "errors": report["errors"],
        "warnings": report["warnings"],
        "summary": report["summary"],
    }, indent=2, ensure_ascii=False))
    raise SystemExit(0 if report["status"] == "PASS" else 1)
