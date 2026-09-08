#!/usr/bin/env python3
"""Merge all completed MC samples listed in the selected sample table."""
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
        help="Continue past missing or stale samples.",
    )
    args = parser.parse_args()

    config = load_run_config(args.config)
    samples = pd.read_parquet(
        mc_sample_table_path(config),
        columns=["sample_id"],
        engine="pyarrow",
    )
    ids = sorted(int(v) for v in samples["sample_id"].tolist())

    if args.start is not None:
        ids = [v for v in ids if v >= args.start]
    if args.end is not None:
        ids = [v for v in ids if v <= args.end]

    failed = []
    for sample_id in ids:
        try:
            merge_one(args.config, sample_id, force=args.force)
        except Exception as exc:
            if not args.skip_incomplete:
                raise
            print("[INCOMPLETE/STALE] sample {}: {}".format(sample_id, exc))
            failed.append(sample_id)

    print(
        "[DONE] requested={}, incomplete_or_stale={}".format(
            len(ids), len(failed)
        )
    )
    if failed:
        print(
            "[INFO] affected sample IDs: {}".format(
                ",".join(str(v) for v in failed)
            )
        )


if __name__ == "__main__":
    main()
