#!/usr/bin/env python3
"""Local entry point. Edit config/run_local.yml; no CLI options are required."""
from pathlib import Path

from src.run_manager import run_from_config

RUN_CONFIG = Path(__file__).resolve().parent.parent / "config" / "run_local.yml"


if __name__ == "__main__":
    result = run_from_config(RUN_CONFIG, server_mode=False)
    print(f"[DONE] {result}")
