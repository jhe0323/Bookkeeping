#!/usr/bin/env python3
"""1 degree server entry point. Configuration: config/run_1deg.yml."""
from pathlib import Path
import sys

CODE_DIR = Path(__file__).resolve().parent.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from src.run_manager import run_from_config  # noqa: E402

RUN_CONFIG = CODE_DIR / "config" / "run_1deg.yml"


if __name__ == "__main__":
    result = run_from_config(
        RUN_CONFIG,
        expected_resolution="1deg",
        server_mode=True,
    )
    print(f"[DONE] {result}")
