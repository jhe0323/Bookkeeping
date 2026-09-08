#!/usr/bin/env python3
"""Merge all completed Monte Carlo samples listed in the sample table."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

CODE_DIR = Path(__file__).resolve().parent.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from src.mc_runner import mc_sample_table_path  # noqa: E402
from src.run_config import load_run_config  # noqa: E402
from tools.merge_mc_bands import merge_one  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/run_mc_1deg.yml")
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--skip-incomplete",
        action="store_true",
        help="Continue past samples whose bands are not all complete.",
    )
    args = parser.parse_args()

    config = load_run_config(args.config, expected_resolution="1deg")
    sample_table = mc_sample_table_path(config)
    samples = pd.read_parquet(sample_table, columns=["sample_id"], engine="pyarrow")
    ids = sorted(int(value) for value in samples["sample_id"].tolist())
    if args.start is not None:
        ids = [value for value in ids if value >= args.start]
    if args.end is not None:
        ids = [value for value in ids if value <= args.end]

    failed = []
    for sample_id in ids:
        try:
            merge_one(args.config, sample_id, force=args.force)
        except Exception as exc:
            if not args.skip_incomplete:
                raise
            print("[INCOMPLETE] sample {}: {}".format(sample_id, exc))
            failed.append(sample_id)

    print("[DONE] requested={}, incomplete={}".format(len(ids), len(failed)))
    if failed:
        print("[INFO] incomplete sample IDs: {}".format(
            ",".join(str(value) for value in failed)
        ))


if __name__ == "__main__":
    main()
