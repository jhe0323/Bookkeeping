#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert IBIS dominant PFT map (vegtype) from 0.5° to 0.25° LUCE-compatible one-hot PFT file.

Input
-----
- vegtype(time=1, level=1, latitude=360, longitude=720)
- dominant PFT code per 0.5° grid cell
- assumes vegtype codes already match your LUCE / yml PFT numbering (1..15)

Output
------
- maxvegetfrac(time_counter=1, veget=15, lat=720, lon=1440)
- one-hot dominant PFT fractions at 0.25°
- each 0.5° cell is split into 2x2 0.25° cells by nearest-neighbor replication
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


N_PFT = 15
FILL_VALUE = np.float32(1.0e20)


def _as_ndarray(var) -> np.ndarray:
    arr = var[:]
    if np.ma.isMaskedArray(arr):
        arr = arr.filled(np.nan)
    return np.asarray(arr)


def wrap_lon_to_180(lon: np.ndarray) -> np.ndarray:
    """Convert lon to [-180, 180)."""
    return ((lon + 180.0) % 360.0) - 180.0


def reorder_longitude_if_needed(data: np.ndarray, lon: np.ndarray):
    """
    Reorder longitude axis from [0,360) to [-180,180).
    data shape: (..., lon)
    """
    lon180 = wrap_lon_to_180(lon)
    order = np.argsort(lon180)
    return data[..., order], lon180[order]


def validate_input(lat: np.ndarray, lon: np.ndarray, veg: np.ndarray):
    if lat.ndim != 1 or lon.ndim != 1:
        raise ValueError("latitude and longitude must be 1D.")

    if veg.ndim != 4:
        raise ValueError(f"Expected vegtype ndim=4, got {veg.ndim}")

    if veg.shape[0] < 1 or veg.shape[1] < 1:
        raise ValueError(f"Unexpected vegtype shape: {veg.shape}")

    if lat.size != 360 or lon.size != 720:
        raise ValueError(
            f"Expected 0.5° global grid (lat=360, lon=720), "
            f"got lat={lat.size}, lon={lon.size}"
        )


def expand_dominant_to_025deg(
    veg2d: np.ndarray,
    lat_05: np.ndarray,
    lon_05: np.ndarray,
    *,
    n_pft: int = N_PFT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Expand 0.5° dominant PFT map to 0.25° by nearest-neighbor 2x2 replication.

    Parameters
    ----------
    veg2d : ndarray
        Shape (Y=360, X=720), dominant PFT codes.
    lat_05 : ndarray
        Shape (360,), latitude center coordinates.
    lon_05 : ndarray
        Shape (720,), longitude center coordinates.

    Returns
    -------
    out_onehot : ndarray
        Shape (1, n_pft, 720, 1440), one-hot dominant PFT fractions.
    lat_025 : ndarray
        Shape (720,), 0.25° latitude center coordinates.
    lon_025 : ndarray
        Shape (1440,), 0.25° longitude center coordinates.
    """
    Y, X = veg2d.shape
    if Y != 360 or X != 720:
        raise ValueError(f"Expected (360,720), got {veg2d.shape}")

    # 0.5° dominant class -> 0.25° dominant class
    veg025 = np.repeat(np.repeat(veg2d, 2, axis=0), 2, axis=1)

    # Generate 0.25° center coordinates from 0.5° centers
    dlat = float(np.median(np.diff(lat_05)))
    dlon = float(np.median(np.diff(lon_05)))

    lat_025 = np.empty(Y * 2, dtype=np.float32)
    lon_025 = np.empty(X * 2, dtype=np.float32)

    lat_025[0::2] = lat_05 - dlat / 4.0
    lat_025[1::2] = lat_05 + dlat / 4.0

    lon_025[0::2] = lon_05 - dlon / 4.0
    lon_025[1::2] = lon_05 + dlon / 4.0

    # dominant class -> one-hot maxvegetfrac
    out_onehot = np.zeros((1, n_pft, Y * 2, X * 2), dtype=np.float32)

    valid = np.isfinite(veg025) & (veg025 >= 1) & (veg025 <= n_pft)

    for p in range(1, n_pft + 1):
        out_onehot[0, p - 1, :, :] = ((veg025 == p) & valid).astype(np.float32)

    return out_onehot, lat_025, lon_025


def write_output(
    out_path: str | Path,
    onehot: np.ndarray,
    lat_025: np.ndarray,
    lon_025: np.ndarray,
):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    T, C, Y, X = onehot.shape

    with Dataset(out_path, "w", format="NETCDF4") as ds:
        ds.createDimension("time_counter", T)
        ds.createDimension("veget", C)
        ds.createDimension("lat", Y)
        ds.createDimension("lon", X)

        vtime = ds.createVariable("time_counter", "f4", ("time_counter",))
        vveg = ds.createVariable("veget", "i4", ("veget",))
        vlat = ds.createVariable("lat", "f4", ("lat",))
        vlon = ds.createVariable("lon", "f4", ("lon",))

        vpft = ds.createVariable(
            "maxvegetfrac",
            "f4",
            ("time_counter", "veget", "lat", "lon"),
            zlib=True,
            complevel=4,
            fill_value=FILL_VALUE,
            chunksizes=(1, min(C, 15), min(Y, 180), min(X, 360)),
        )

        vtime[:] = np.array([0], dtype=np.float32)
        vveg[:] = np.arange(1, C + 1, dtype=np.int32)
        vlat[:] = lat_025
        vlon[:] = lon_025
        vpft[:] = onehot

        vlat.units = "degrees_north"
        vlon.units = "degrees_east"
        vtime.long_name = "dummy time index"
        vveg.long_name = "LUCE PFT classes"

        vpft.long_name = (
            "0.25-degree dominant-PFT one-hot fractions converted from IBIS 0.5-degree vegtype"
        )

        ds.description = (
            "LUCE-compatible 0.25-degree dominant PFT map converted from IBIS 0.5-degree vegtype.nc"
        )
        ds.note = (
            "Each 0.5-degree grid cell is split into 2x2 0.25-degree cells. "
            "The dominant PFT code is copied directly without remapping, assuming vegtype codes "
            "already match the LUCE / yml PFT numbering."
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert IBIS dominant vegtype map from 0.5° to LUCE-compatible "
            "0.25° one-hot PFT file."
        )
    )
    parser.add_argument("--input", required=True, help="Path to input PFTmap_IBIS_0.5deg.nc")
    parser.add_argument("--output", required=True, help="Path to output 0.25° netCDF file")
    parser.add_argument("--var", default="vegtype", help="Input variable name, default: vegtype")
    args = parser.parse_args()

    with Dataset(args.input, "r") as ds:
        lon = _as_ndarray(ds.variables["longitude"]).astype(np.float64)
        lat = _as_ndarray(ds.variables["latitude"]).astype(np.float64)
        veg = _as_ndarray(ds.variables[args.var]).astype(np.float64)

    validate_input(lat, lon, veg)

    # take first time, first level
    veg2d = veg[0, 0, :, :]

    # reorder lon to [-180, 180), matching your model convention
    veg2d, lon180 = reorder_longitude_if_needed(veg2d, lon)

    # ensure latitude ascending: south -> north
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        veg2d = veg2d[::-1, :]

    onehot, lat_025, lon_025 = expand_dominant_to_025deg(
        veg2d,
        lat,
        lon180,
        n_pft=N_PFT,
    )

    # quick checks
    s = np.sum(onehot, axis=1)

    print("Output shape:", onehot.shape)
    print("lat first/last:", float(lat_025[0]), float(lat_025[-1]))
    print("lon first/last:", float(lon_025[0]), float(lon_025[-1]))
    print("Global mean sum over veget:", float(np.nanmean(s)))

    if np.any(s > 0):
        print("Land-only mean sum over veget:", float(np.nanmean(s[s > 0])))
        print(
            "Land-only min/max sum over veget:",
            float(np.nanmin(s[s > 0])),
            float(np.nanmax(s[s > 0])),
        )

    print("Cells with no valid PFT:", int(np.sum(s[0] == 0)))

    write_output(args.output, onehot, lat_025, lon_025)
    print(f"Written: {args.output}")


if __name__ == "__main__":
    main()