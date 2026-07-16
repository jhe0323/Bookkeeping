#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert ORCHIDEE 0.25° PFT fraction map to LUCE-coded 1° PFT fraction map.

Input assumptions
-----------------
- Variable name: maxvegetfrac
- Dimensions: (time_counter, veget=15, lat, lon)
- ORCHIDEE PFT codes:
    1  Bare soil
    2  Tropical Broadleaf Evergreen
    3  Tropical Broadleaf Rain green
    4  Temperate Needleleaf Evergreen
    5  Temperate Broadleaf Evergreen
    6  Temperate Broadleaf Summer green
    7  Boreal Needleleaf Evergreen
    8  Boreal Broadleaf Summergreen
    9  Boreal Needleleaf Deciduous
    10 Temperate Natural Grassland (C3)
    11 Natural Grassland (C4)
    12 Crops (C3)
    13 Crops (C4)
    14 Tropical Natural Grassland (C3)
    15 Boreal Natural Grassland (C3)

Mapping to LUCE classes
-----------------------
    1  -> 14
    2  -> 1
    3  -> 2
    4  -> 4
    5  -> 3
    6  -> 5
    7  -> 6
    8  -> 8
    9  -> 7
    10 -> 10
    11 -> 10
    12 -> 10
    13 -> 10
    14 -> 10
    15 -> 13

Notes
-----
- LUCE classes 11 and 12 are kept as valid output classes but may remain zero
  because the ORCHIDEE PFT map does not provide explicit shrubland types.
- LUCE classes 14 and 15 are placeholders in your current setup; 14 receives
  ORCHIDEE bare soil according to the agreed mapping, while 15 may remain zero.
- Aggregation from 0.25° to 1° is performed as an area-weighted mean using
  cos(latitude) weights.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


ORCHIDEE_TO_LUCE = {
    1: 14,
    2: 1,
    3: 2,
    4: 4,
    5: 3,
    6: 5,
    7: 6,
    8: 8,
    9: 7,
    10: 10,
    11: 10,
    12: 10,
    13: 10,
    14: 10,
    15: 13,
}

LUCE_CLASS_NAMES = {
    1: "Tropical Evergreen Forest/Woodland",
    2: "Tropical Deciduous Forest/Woodland",
    3: "Temperate Broadleaf Evergreen Forest/Woodland",
    4: "Temperate Needleleaf Evergreen Forest/Woodland",
    5: "Temperate Deciduous Forest/Woodland",
    6: "Boreal Evergreen Forest/Woodland",
    7: "Boreal Deciduous Forest/Woodland",
    8: "Evergreen/Deciduous Mixed Forest/Woodland",
    9: "Savanna",
    10: "Grassland/Steppe",
    11: "Dense Shrubland",
    12: "Open Shrubland",
    13: "Tundra",
    14: "Desert",
    15: "Polar Desert/Rock/Ice",
}


def _as_ndarray(var) -> np.ndarray:
    """Read a netCDF variable as ndarray with masked values converted to NaN."""
    arr = var[:]
    if np.ma.isMaskedArray(arr):
        arr = arr.filled(np.nan)
    return np.asarray(arr)


def validate_grid(lat: np.ndarray, lon: np.ndarray) -> None:
    """Check that the source grid looks like a regular 0.25° global grid."""
    if lat.ndim != 1 or lon.ndim != 1:
        raise ValueError("lat and lon must be 1D arrays.")
    if lat.size % 4 != 0 or lon.size % 4 != 0:
        raise ValueError(
            f"Expected source lat/lon lengths divisible by 4 for 1° aggregation, "
            f"got lat={lat.size}, lon={lon.size}."
        )

    if lat.size >= 2:
        dlat = np.diff(lat)
        if not np.allclose(np.abs(dlat), 0.25, atol=1e-6):
            raise ValueError("Latitude spacing is not 0.25°.")
    if lon.size >= 2:
        dlon = np.diff(lon)
        # Some files may be 0.25 increments without wrapping at 360 exactly.
        if not np.allclose(dlon, 0.25, atol=1e-6):
            raise ValueError("Longitude spacing is not 0.25°.")


def convert_orchidee_to_luce(frac: np.ndarray) -> np.ndarray:
    """
    Convert ORCHIDEE fractions to LUCE fractions.

    Parameters
    ----------
    frac : np.ndarray
        Shape (T, 15, Y, X)

    Returns
    -------
    np.ndarray
        Shape (T, 15, Y, X), LUCE-coded fractions.
    """
    if frac.ndim != 4 or frac.shape[1] != 15:
        raise ValueError(f"Expected input shape (T, 15, Y, X), got {frac.shape}")

    T, _, Y, X = frac.shape
    luce = np.zeros((T, 15, Y, X), dtype=np.float32)

    for orchidee_pft in range(1, 16):
        luce_pft = ORCHIDEE_TO_LUCE[orchidee_pft]
        luce[:, luce_pft - 1, :, :] += frac[:, orchidee_pft - 1, :, :].astype(np.float32)

    return luce


def aggregate_to_1deg_area_weighted(
    data: np.ndarray,
    lat_025: np.ndarray,
    lon_025: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Aggregate 0.25° data to 1° using area-weighted means.

    Parameters
    ----------
    data : np.ndarray
        Shape (T, C, Y, X)
    lat_025 : np.ndarray
        Shape (Y,)
    lon_025 : np.ndarray
        Shape (X,)

    Returns
    -------
    data_1deg : np.ndarray
        Shape (T, C, Y1, X1)
    lat_1deg : np.ndarray
        Shape (Y1,)
    lon_1deg : np.ndarray
        Shape (X1,)
    """
    T, C, Y, X = data.shape
    if Y % 4 != 0 or X % 4 != 0:
        raise ValueError(f"Cannot aggregate to 1° because Y={Y}, X={X} are not multiples of 4.")

    Y1 = Y // 4
    X1 = X // 4

    # Area weights vary with latitude only.
    lat_weights = np.cos(np.deg2rad(lat_025)).astype(np.float64)
    lat_weights = np.clip(lat_weights, 0.0, None)

    # Build a 2D weight field. All longitudes in a latitude ring share the same weight.
    w2d = lat_weights[:, None] * np.ones((1, X), dtype=np.float64)

    # Ocean/invalid cells in the source map are all-zero fractions. We keep them in the average,
    # which naturally yields zero for fully oceanic 1° cells.
    data64 = data.astype(np.float64)

    # Reshape into 4x4 source blocks.
    data_rs = data64.reshape(T, C, Y1, 4, X1, 4)
    w_rs = w2d.reshape(Y1, 4, X1, 4)

    weighted_sum = np.nansum(data_rs * w_rs[None, None, :, :, :, :], axis=(3, 5))
    weight_sum = np.sum(w_rs, axis=(1, 3))  # shape (Y1, X1)

    out = weighted_sum / weight_sum[None, None, :, :]
    out = out.astype(np.float32)

    # 1° cell-center coordinates from grouped 0.25° centers.
    lat_1deg = lat_025.reshape(Y1, 4).mean(axis=1).astype(np.float32)
    lon_1deg = lon_025.reshape(X1, 4).mean(axis=1).astype(np.float32)

    return out, lat_1deg, lon_1deg


def write_output_nc(
    out_path: str | Path,
    data_1deg: np.ndarray,
    lat_1deg: np.ndarray,
    lon_1deg: np.ndarray,
    time_vals: np.ndarray,
    src_description: str = "",
) -> None:
    """Write the LUCE-coded 1° map to netCDF."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    T, C, Y, X = data_1deg.shape
    fill_value = np.float32(1.0e20)

    with Dataset(out_path, "w", format="NETCDF4") as ds:
        ds.createDimension("time_counter", T)
        ds.createDimension("veget", C)
        ds.createDimension("lat", Y)
        ds.createDimension("lon", X)

        lat_var = ds.createVariable("lat", "f4", ("lat",))
        lon_var = ds.createVariable("lon", "f4", ("lon",))
        time_var = ds.createVariable("time_counter", "f4", ("time_counter",))
        veget_var = ds.createVariable("veget", "i4", ("veget",))
        frac_var = ds.createVariable(
            "maxvegetfrac",
            "f4",
            ("time_counter", "veget", "lat", "lon"),
            zlib=True,
            complevel=4,
            fill_value=fill_value,
            chunksizes=(1, min(15, C), min(90, Y), min(180, X)),
        )

        lat_var[:] = lat_1deg
        lon_var[:] = lon_1deg
        time_var[:] = time_vals
        veget_var[:] = np.arange(1, C + 1, dtype=np.int32)

        frac_to_write = np.where(np.isfinite(data_1deg), data_1deg, fill_value)
        frac_var[:] = frac_to_write

        lat_var.long_name = "Latitude"
        lat_var.units = "degrees_north"
        lat_var.axis = "Y"

        lon_var.long_name = "Longitude"
        lon_var.units = "degrees_east"
        lon_var.axis = "X"

        time_var.axis = "T"
        time_var.long_name = "Time"

        veget_var.long_name = "LUCE vegetation classes"
        veget_var.units = "-"
        veget_var.axis = "Z"
        veget_var.class_names = " | ".join(
            f"{k}:{v}" for k, v in LUCE_CLASS_NAMES.items()
        )

        frac_var.long_name = "LUCE vegetation class fractions aggregated to 1 degree"
        frac_var.units = "-"
        frac_var.mapping_note = (
            "Converted from ORCHIDEE PFT fractions using a user-defined ORCHIDEE->LUCE mapping."
        )

        ds.description = (
            "1-degree LUCE-coded PFT fraction map converted from ORCHIDEE maxvegetfrac."
        )
        ds.source_description = src_description
        ds.mapping = "; ".join(f"{k}->{v}" for k, v in ORCHIDEE_TO_LUCE.items())
        ds.note = (
            "LUCE classes 11 and 12 may remain zero because the source ORCHIDEE map "
            "does not provide explicit shrubland classes. Class 14 receives ORCHIDEE bare soil "
            "according to the agreed mapping."
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert ORCHIDEE 0.25° PFT map to LUCE-coded 1° PFT map."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input ORCHIDEE netCDF file, e.g. PFTmap_2000_ipsl.nc",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output LUCE-coded 1° netCDF file",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    with Dataset(input_path, "r") as ds:
        lat = _as_ndarray(ds.variables["lat"]).astype(np.float64)
        lon = _as_ndarray(ds.variables["lon"]).astype(np.float64)
        time_vals = _as_ndarray(ds.variables["time_counter"]).astype(np.float32)
        frac = _as_ndarray(ds.variables["maxvegetfrac"]).astype(np.float32)
        src_description = getattr(ds, "description", "")

    validate_grid(lat, lon)

    # Source shape should be (T, 15, 720, 1440)
    luce_025 = convert_orchidee_to_luce(frac)

    # Optional sanity checks
    sum_src = np.nansum(frac, axis=1)
    sum_luce_025 = np.nansum(luce_025, axis=1)

    print("Input shape:", frac.shape)
    print("LUCE 0.25° shape:", luce_025.shape)
    print("Source global mean sum over veget:", float(np.nanmean(sum_src)))
    print("Source land-only mean sum over veget:", float(np.nanmean(sum_src[sum_src > 0])))
    print("Converted 0.25° land-only mean sum over veget:", float(np.nanmean(sum_luce_025[sum_luce_025 > 0])))

    luce_1deg, lat_1deg, lon_1deg = aggregate_to_1deg_area_weighted(luce_025, lat, lon)
    sum_luce_1deg = np.nansum(luce_1deg, axis=1)

    print("LUCE 1° shape:", luce_1deg.shape)
    print("1° global mean sum over veget:", float(np.nanmean(sum_luce_1deg)))
    if np.any(sum_luce_1deg > 0):
        print("1° land-only mean sum over veget:", float(np.nanmean(sum_luce_1deg[sum_luce_1deg > 0])))
        print("1° land-only min/max sum over veget:",
              float(np.nanmin(sum_luce_1deg[sum_luce_1deg > 0])),
              float(np.nanmax(sum_luce_1deg[sum_luce_1deg > 0])))

    write_output_nc(
        out_path=output_path,
        data_1deg=luce_1deg,
        lat_1deg=lat_1deg,
        lon_1deg=lon_1deg,
        time_vals=time_vals,
        src_description=src_description,
    )

    print(f"Written: {output_path}")


if __name__ == "__main__":
    main()
