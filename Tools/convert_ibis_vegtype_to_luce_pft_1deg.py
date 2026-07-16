#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert IBIS dominant PFT map (vegtype) to a 1-degree LUCE-compatible one-hot PFT file.

Input
-----
- vegtype(time=1, level=1, latitude=360, longitude=720)
- dominant PFT code per 0.5° grid cell
- assumes vegtype codes already match your LUCE / yml PFT numbering (1..15)

Output
------
- maxvegetfrac(time_counter=1, veget=15, lat=180, lon=360)
- one-hot dominant PFT fractions at 1°
"""

from __future__ import annotations

import argparse
from pathlib import Path
from collections import defaultdict

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
    Reorder longitude axis from [0,360) to [-180,180) if needed.
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
            f"Expected 0.5° global grid (lat=360, lon=720), got lat={lat.size}, lon={lon.size}"
        )


def aggregate_dominant_to_1deg(
    veg2d: np.ndarray,
    lat_05: np.ndarray,
    lon_05: np.ndarray,
    *,
    n_pft: int = N_PFT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Aggregate 0.5° dominant PFT map to 1° dominant PFT map using majority vote.
    Ties are broken by cosine-latitude weighted count.

    Parameters
    ----------
    veg2d : ndarray
        Shape (Y=360, X=720), dominant PFT codes
    lat_05 : ndarray
        Shape (360,)
    lon_05 : ndarray
        Shape (720,)

    Returns
    -------
    out_onehot : ndarray
        Shape (1, n_pft, 180, 360), one-hot dominant PFT fractions
    lat_1deg : ndarray
        Shape (180,)
    lon_1deg : ndarray
        Shape (360,)
    """
    Y, X = veg2d.shape
    if Y != 360 or X != 720:
        raise ValueError(f"Expected (360,720), got {veg2d.shape}")

    Y1, X1 = Y // 2, X // 2

    out_dom = np.zeros((Y1, X1), dtype=np.int16)

    # block-center coordinates
    lat_1deg = lat_05.reshape(Y1, 2).mean(axis=1).astype(np.float32)
    lon_1deg = lon_05.reshape(X1, 2).mean(axis=1).astype(np.float32)

    for iy in range(Y1):
        for ix in range(X1):
            block = veg2d[iy*2:(iy+1)*2, ix*2:(ix+1)*2]
            block_lat = lat_05[iy*2:(iy+1)*2]

            # flatten
            vals = block.reshape(-1)

            # valid codes = 1..n_pft
            valid = np.isfinite(vals) & (vals >= 1) & (vals <= n_pft)
            vals = vals[valid].astype(int)

            if vals.size == 0:
                out_dom[iy, ix] = 0
                continue

            # majority count
            uniq, cnt = np.unique(vals, return_counts=True)
            maxcnt = cnt.max()
            candidates = uniq[cnt == maxcnt]

            if candidates.size == 1:
                out_dom[iy, ix] = int(candidates[0])
                continue

            # tie-break by latitude area weights
            # 2x2 block: weights depend only on row latitude
            row_weights = np.cos(np.deg2rad(block_lat))
            row_weights = np.clip(row_weights, 0.0, None)
            wblock = np.repeat(row_weights[:, None], 2, axis=1).reshape(-1)[valid]

            score = defaultdict(float)
            for v, w in zip(vals, wblock):
                if v in candidates:
                    score[int(v)] += float(w)

            best = max(candidates, key=lambda k: score[int(k)])
            out_dom[iy, ix] = int(best)

    # convert dominant -> one-hot
    out_onehot = np.zeros((1, n_pft, Y1, X1), dtype=np.float32)
    for p in range(1, n_pft + 1):
        out_onehot[0, p - 1, :, :] = (out_dom == p).astype(np.float32)

    return out_onehot, lat_1deg, lon_1deg


def write_output(
    out_path: str | Path,
    onehot: np.ndarray,
    lat_1deg: np.ndarray,
    lon_1deg: np.ndarray,
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
            chunksizes=(1, min(C, 15), min(Y, 90), min(X, 180)),
        )

        vtime[:] = np.array([0], dtype=np.float32)
        vveg[:] = np.arange(1, C + 1, dtype=np.int32)
        vlat[:] = lat_1deg
        vlon[:] = lon_1deg
        vpft[:] = onehot

        vlat.units = "degrees_north"
        vlon.units = "degrees_east"
        vtime.long_name = "dummy time index"
        vveg.long_name = "LUCE PFT classes"
        vpft.long_name = "1-degree dominant-PFT one-hot fractions converted from IBIS vegtype"

        ds.description = "LUCE-compatible 1-degree dominant PFT map converted from IBIS vegtype.nc"
        ds.note = (
            "Each 1-degree grid cell has a single dominant PFT encoded as one-hot maxvegetfrac. "
            "0.5-degree source dominant classes are aggregated by 2x2 majority vote, "
            "with cosine-latitude weighting used for tie-breaking."
        )


def main():
    parser = argparse.ArgumentParser(
        description="Convert IBIS dominant vegtype map to LUCE-compatible 1° one-hot PFT file."
    )
    parser.add_argument("--input", required=True, help="Path to vegtype.nc")
    parser.add_argument("--output", required=True, help="Path to output netCDF file")
    parser.add_argument("--var", default="vegtype", help="Input variable name (default: vegtype)")
    args = parser.parse_args()

    with Dataset(args.input, "r") as ds:
        lon = _as_ndarray(ds.variables["longitude"]).astype(np.float64)
        lat = _as_ndarray(ds.variables["latitude"]).astype(np.float64)
        veg = _as_ndarray(ds.variables[args.var]).astype(np.float64)

    validate_input(lat, lon, veg)

    # take first time, first level
    veg2d = veg[0, 0, :, :]

    # reorder lon to [-180, 180)
    veg2d, lon180 = reorder_longitude_if_needed(veg2d, lon)

    # ensure latitude ascending (south -> north), matching your current internal convention
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        veg2d = veg2d[::-1, :]

    onehot, lat_1deg, lon_1deg = aggregate_dominant_to_1deg(veg2d, lat, lon180, n_pft=N_PFT)

    # quick checks
    s = np.sum(onehot, axis=1)
    print("Output shape:", onehot.shape)
    print("Global mean sum over veget:", float(np.nanmean(s)))
    if np.any(s > 0):
        print("Land-only mean sum over veget:", float(np.nanmean(s[s > 0])))
        print("Land-only min/max sum over veget:", float(np.nanmin(s[s > 0])), float(np.nanmax(s[s > 0])))

    write_output(args.output, onehot, lat_1deg, lon_1deg)
    print(f"Written: {args.output}")


if __name__ == "__main__":
    main()