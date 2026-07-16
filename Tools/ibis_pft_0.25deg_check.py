#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Check LUCE-compatible PFT file.

Checks:
1. required dimensions and variables exist
2. maxvegetfrac shape is correct
3. PFT fractions are finite and non-negative
4. each valid land grid has at least one PFT
5. PFT fractions sum to approximately 1 over veget
6. each grid is one-hot or fractional
7. reports which PFT types exist globally
"""

from __future__ import annotations

import argparse
import numpy as np
from netCDF4 import Dataset


def check_pft_file(
    nc_path: str,
    varname: str = "maxvegetfrac",
    tol: float = 1.0e-5,
    expected_nveg: int = 15,
):
    with Dataset(nc_path, "r") as ds:
        print("File:", nc_path)

        required_vars = ["lat", "lon", "veget", varname]
        for v in required_vars:
            if v not in ds.variables:
                raise ValueError(f"Missing variable: {v}")

        lat = ds.variables["lat"][:]
        lon = ds.variables["lon"][:]
        veget = ds.variables["veget"][:]
        pft = ds.variables[varname][:]

        if np.ma.isMaskedArray(pft):
            pft = pft.filled(np.nan)

        pft = np.asarray(pft, dtype=np.float64)

        print("lat shape:", lat.shape)
        print("lon shape:", lon.shape)
        print("veget:", veget)
        print("PFT variable shape:", pft.shape)

        if pft.ndim != 4:
            raise ValueError(f"{varname} should be 4D: time, veget, lat, lon")

        T, C, Y, X = pft.shape

        if C != expected_nveg:
            raise ValueError(f"Expected {expected_nveg} PFTs, got {C}")

        if Y != len(lat) or X != len(lon):
            raise ValueError("PFT lat/lon dimensions do not match coordinate lengths")

        print("\n--- Coordinate check ---")
        print("lat first/last:", float(lat[0]), float(lat[-1]))
        print("lon first/last:", float(lon[0]), float(lon[-1]))

        lat_ascending = lat[0] < lat[-1]
        lon_ascending = lon[0] < lon[-1]
        print("lat ascending:", bool(lat_ascending))
        print("lon ascending:", bool(lon_ascending))

        if len(lat) > 1:
            print("lat resolution:", float(np.median(np.diff(lat))))
        if len(lon) > 1:
            print("lon resolution:", float(np.median(np.diff(lon))))

        arr = pft[0, :, :, :]

        print("\n--- Value check ---")
        finite_mask = np.isfinite(arr)
        print("finite values:", int(np.sum(finite_mask)), "/", arr.size)

        if np.any(arr[finite_mask] < -tol):
            bad = np.sum(arr[finite_mask] < -tol)
            print("WARNING: negative PFT fractions found:", int(bad))
        else:
            print("negative values: none")

        if np.any(arr[finite_mask] > 1.0 + tol):
            bad = np.sum(arr[finite_mask] > 1.0 + tol)
            print("WARNING: PFT fractions > 1 found:", int(bad))
        else:
            print("values > 1: none")

        pft_sum = np.nansum(arr, axis=0)

        land_mask = pft_sum > tol
        empty_mask = pft_sum <= tol
        bad_sum_mask = land_mask & (np.abs(pft_sum - 1.0) > tol)

        print("\n--- Sum-over-PFT check ---")
        print("total grid cells:", Y * X)
        print("valid land cells:", int(np.sum(land_mask)))
        print("empty / ocean / no-PFT cells:", int(np.sum(empty_mask)))
        print("land fraction:", float(np.sum(land_mask) / (Y * X)))

        if np.any(land_mask):
            print("land-cell PFT sum min:", float(np.nanmin(pft_sum[land_mask])))
            print("land-cell PFT sum max:", float(np.nanmax(pft_sum[land_mask])))
            print("land-cell PFT sum mean:", float(np.nanmean(pft_sum[land_mask])))

        print("bad land cells where sum != 1:", int(np.sum(bad_sum_mask)))

        if np.sum(bad_sum_mask) > 0:
            iy, ix = np.where(bad_sum_mask)
            print("\nExamples of bad cells:")
            for k in range(min(10, len(iy))):
                print(
                    f"  lat={float(lat[iy[k]])}, lon={float(lon[ix[k]])}, "
                    f"sum={float(pft_sum[iy[k], ix[k]])}"
                )

        print("\n--- PFT type check ---")
        pft_area_cells = np.nansum(arr > tol, axis=(1, 2))

        existing_pfts = []
        missing_pfts = []

        for i in range(C):
            pft_id = int(veget[i]) if i < len(veget) else i + 1
            count = int(pft_area_cells[i])
            if count > 0:
                existing_pfts.append(pft_id)
            else:
                missing_pfts.append(pft_id)
            print(f"PFT {pft_id:2d}: cells = {count}")

        print("Existing PFTs:", existing_pfts)
        print("Missing PFTs:", missing_pfts)

        print("\n--- One-hot check ---")
        nonzero_count = np.sum(arr > tol, axis=0)

        onehot_land = land_mask & (nonzero_count == 1)
        fractional_land = land_mask & (nonzero_count > 1)
        strange_land = land_mask & (nonzero_count == 0)

        print("one-hot land cells:", int(np.sum(onehot_land)))
        print("fractional/mixed land cells:", int(np.sum(fractional_land)))
        print("strange land cells with sum>0 but no PFT>tol:", int(np.sum(strange_land)))

        if np.sum(fractional_land) == 0 and np.sum(bad_sum_mask) == 0:
            print("\nRESULT: PASS. File looks correct for dominant one-hot PFT input.")
        elif np.sum(bad_sum_mask) == 0:
            print("\nRESULT: PASS for fractional PFT input, but not strictly one-hot.")
        else:
            print("\nRESULT: WARNING. Some land cells do not sum to 1.")


def main():
    parser = argparse.ArgumentParser(description="Check LUCE-compatible PFT NetCDF file.")
    parser.add_argument("--input", required=True, help="Path to PFT NetCDF file")
    parser.add_argument("--var", default="maxvegetfrac", help="PFT variable name")
    parser.add_argument("--tol", type=float, default=1.0e-5, help="Tolerance for sum check")
    parser.add_argument("--nveg", type=int, default=15, help="Expected number of PFT classes")
    args = parser.parse_args()

    check_pft_file(
        nc_path=args.input,
        varname=args.var,
        tol=args.tol,
        expected_nveg=args.nveg,
    )


if __name__ == "__main__":
    main()