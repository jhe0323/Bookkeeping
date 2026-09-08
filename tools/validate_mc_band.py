#!/usr/bin/env python3
"""Compare one summary-only MC band against an existing gridded NetCDF band.

Use this on the baseline realization before launching a large ensemble.  The
script spatially sums each selected NetCDF variable and checks that it matches
the MC band Parquet time series.
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
    args = parser.parse_args()

    mc_path = Path(args.mc_band).resolve()
    nc_path = Path(args.grid_nc).resolve()
    mc = pd.read_parquet(mc_path, engine="pyarrow")

    failures = []
    checked = 0
    with nc.Dataset(nc_path, "r") as ds:
        if "time" not in ds.variables:
            raise ValueError("Reference NetCDF has no time variable")
        years = np.asarray(ds.variables["time"][:], dtype=int)
        if not np.array_equal(years, mc["year"].to_numpy(dtype=int)):
            raise RuntimeError("Year coordinates differ between MC band and NetCDF")

        for variable in mc.columns[1:]:
            if variable not in ds.variables:
                continue
            raw = np.ma.asarray(ds.variables[variable][:])
            data = np.asarray(raw.filled(0.0), dtype=np.float64)
            if data.ndim != 3:
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
                failures.append((variable, max_abs, scale))

    if checked == 0:
        raise RuntimeError("No common time x lat x lon variables were found to compare")
    if failures:
        raise RuntimeError(
            "Baseline equivalence failed for {} variables: {}".format(
                len(failures), ", ".join(item[0] for item in failures)
            )
        )
    print("[PASS] baseline MC summary is equivalent for {} variables".format(checked))


if __name__ == "__main__":
    main()
