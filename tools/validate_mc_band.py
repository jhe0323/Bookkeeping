#!/usr/bin/env python3
"""Compare one baseline summary-only MC band against deterministic NetCDF.

Every MC output variable is required to exist in the deterministic reference
unless it is explicitly excluded with ``--allow-missing-variable``. This keeps
the baseline-equivalence check from silently passing after comparing only a
small subset of variables.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mc-band", required=True, help="MC band_XXX.parquet")
    parser.add_argument("--grid-nc", required=True, help="Reference deterministic rankXXX.nc")
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--atol", type=float, default=1.0e-5)
    parser.add_argument(
        "--allow-missing-variable",
        action="append",
        default=[],
        help=(
            "Variable allowed to be absent from the deterministic NetCDF. "
            "May be repeated. Default: none."
        ),
    )
    args = parser.parse_args()

    mc_path = Path(args.mc_band).resolve()
    nc_path = Path(args.grid_nc).resolve()
    mc = pd.read_parquet(mc_path, engine="pyarrow")
    allowed_missing = set(str(v) for v in args.allow_missing_variable)

    failures = []
    missing = []
    checked = 0
    with nc.Dataset(nc_path, "r") as ds:
        if "time" not in ds.variables:
            raise ValueError("Reference NetCDF has no time variable")
        years = np.asarray(ds.variables["time"][:], dtype=int)
        if not np.array_equal(years, mc["year"].to_numpy(dtype=int)):
            raise RuntimeError("Year coordinates differ between MC band and NetCDF")

        for variable in mc.columns[1:]:
            if variable not in ds.variables:
                if variable not in allowed_missing:
                    missing.append(variable)
                continue

            raw = np.ma.asarray(ds.variables[variable][:])
            data = np.asarray(raw.filled(0.0), dtype=np.float64)
            if data.ndim != 3:
                failures.append((variable, "reference variable is not time x lat x lon"))
                continue

            reference = data.sum(axis=(1, 2), dtype=np.float64)
            candidate = mc[variable].to_numpy(dtype=np.float64)
            difference = candidate - reference
            max_abs = float(np.max(np.abs(difference)))
            scale = float(np.max(np.abs(reference)))
            ok = np.allclose(candidate, reference, rtol=args.rtol, atol=args.atol)
            checked += 1
            print(
                "{:<32s} {} max_abs={:.6g} scale={:.6g}".format(
                    variable,
                    "PASS" if ok else "FAIL",
                    max_abs,
                    scale,
                )
            )
            if not ok:
                failures.append((variable, "max_abs={:.6g}".format(max_abs)))

    if missing:
        raise RuntimeError(
            "Baseline equivalence failed: deterministic NetCDF is missing {} "
            "required MC variables: {}".format(
                len(missing), ", ".join(missing)
            )
        )
    if checked == 0:
        raise RuntimeError("No common time x lat x lon variables were found to compare")
    if failures:
        raise RuntimeError(
            "Baseline equivalence failed for {} variables: {}".format(
                len(failures),
                ", ".join("{} ({})".format(v, reason) for v, reason in failures),
            )
        )
    print("[PASS] baseline MC summary is equivalent for {} variables".format(checked))


if __name__ == "__main__":
    main()
