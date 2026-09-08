#!/usr/bin/env python3
"""0.25-degree Monte Carlo summary-only server entry point."""
from pathlib import Path
import os
import sys

CODE_DIR = Path(__file__).resolve().parent.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from src.mc_runner import run_mc_band  # noqa: E402

RUN_CONFIG = Path(
    os.environ.get(
        "RUN_CONFIG",
        str(CODE_DIR / "config" / "run_mc_025deg.yml"),
    )
).resolve()

if __name__ == "__main__":
    result = run_mc_band(
        RUN_CONFIG,
        expected_resolution="025deg",
    )
    print("[DONE] {}".format(result))
