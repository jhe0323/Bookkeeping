#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nc_viewer_updated.py

适配当前 bookkeeping 输出的 NetCDF 快速查看脚本。

支持：
1. summary_*.nc: 3D(time, lat, lon)
2. 差值切片文件: 2D(lat, lon)，例如 atmosphere_1960_2020.nc
3. 自动识别变量、坐标、年份
4. 差值图按年份差自动年均化
5. 碳变量自动 /1e9 显示为 Gt C 或 Gt C/yr
6. 0 值可显示为白色，适合看稀疏区域

示例：
python nc_viewer_updated.py summary_1deg.global.nc --list
python nc_viewer_updated.py summary_1deg.global.nc --var atmosphere --year 2020
python nc_viewer_updated.py summary_1deg.global.nc --var atmosphere --diff-year 2020 1960 --save atmosphere_1960_2020.png
python nc_viewer_updated.py atmosphere_1960_2020.nc --var atmosphere --bbox china
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm

try:
    from netCDF4 import Dataset
except Exception as e:
    raise SystemExit("ERROR: netCDF4 is required. Install with `pip install netCDF4`.\n" + str(e))


LAT_CANDIDATES = ["lat", "latitude", "LAT", "nav_lat", "y", "Y"]
LON_CANDIDATES = ["lon", "longitude", "LON", "nav_lon", "x", "X"]
TIME_CANDIDATES = ["time", "Time", "TIME", "t"]
COORD_NAMES = set(LAT_CANDIDATES + LON_CANDIDATES + TIME_CANDIDATES)

BBOX_PRESETS = {
    "global": (-180.0, 180.0, -90.0, 90.0),
    "china": (73.0, 135.0, 18.0, 54.0),
    "loess": (100.0, 114.0, 33.0, 41.5),
    "xinjiang": (73.0, 97.5, 34.0, 50.5),
}

NON_CARBON_KEYWORDS = (
    "area", "frac", "fraction", "mask", "done", "land", "cell",
)


def _find_first(ds: Dataset, names: Sequence[str]) -> Optional[str]:
    for n in names:
        if n in ds.variables:
            return n
    return None


def _find_by_attr(ds: Dataset, attr: str, value: str) -> Optional[str]:
    for k, v in ds.variables.items():
        if getattr(v, attr, "").lower() == value.lower():
            return k
    return None


def find_lat_lon(ds: Dataset) -> Tuple[str, str]:
    lat_name = _find_first(ds, LAT_CANDIDATES) or _find_by_attr(ds, "standard_name", "latitude")
    lon_name = _find_first(ds, LON_CANDIDATES) or _find_by_attr(ds, "standard_name", "longitude")
    if lat_name is None or lon_name is None:
        raise ValueError("Could not find lat/lon variables.")
    return lat_name, lon_name


def find_time(ds: Dataset) -> Optional[str]:
    return _find_first(ds, TIME_CANDIDATES) or _find_by_attr(ds, "standard_name", "time")


def is_data_var(ds: Dataset, name: str) -> bool:
    if name in COORD_NAMES:
        return False
    var = ds.variables[name]
    dims = tuple(var.dimensions)
    if len(dims) == 0:
        return False
    lat_name, lon_name = find_lat_lon(ds)
    time_name = find_time(ds)
    if len(dims) == 2 and lat_name in dims and lon_name in dims:
        return True
    if len(dims) == 3 and time_name in dims and lat_name in dims and lon_name in dims:
        return True
    return False


def list_vars(ds: Dataset) -> List[str]:
    return sorted([name for name in ds.variables if is_data_var(ds, name)])


def get_time_values(ds: Dataset) -> Optional[np.ndarray]:
    time_name = find_time(ds)
    if time_name is None or time_name not in ds.variables:
        return None
    return np.asarray(ds.variables[time_name][:]).astype(int)


def year_to_index(ds: Dataset, year: int) -> int:
    tvals = get_time_values(ds)
    if tvals is None:
        raise ValueError("No time variable found; cannot use --year or --diff-year.")
    matches = np.where(tvals == year)[0]
    if len(matches) == 0:
        raise ValueError(f"Year {year} not found. Available range: {tvals.min()} .. {tvals.max()}")
    return int(matches[0])


def get_2d_field(ds: Dataset, varname: str, tidx: Optional[int] = None) -> np.ndarray:
    if varname not in ds.variables:
        raise ValueError(f"Variable '{varname}' not found. Use --list to see available variables.")

    v = ds.variables[varname]
    dims = list(v.dimensions)
    lat_name, lon_name = find_lat_lon(ds)
    time_name = find_time(ds)

    if len(dims) == 2:
        axes = [dims.index(lat_name), dims.index(lon_name)]
        arr = np.asarray(v[:], dtype=float)
        arr = np.transpose(arr, axes=axes)
        return arr

    if len(dims) == 3:
        if tidx is None:
            tidx = 0
        axes = [dims.index(time_name), dims.index(lat_name), dims.index(lon_name)]
        arr3 = np.asarray(v[:], dtype=float)
        arr3 = np.transpose(arr3, axes=axes)
        if tidx < 0 or tidx >= arr3.shape[0]:
            raise IndexError(f"tidx={tidx} out of range 0..{arr3.shape[0]-1}")
        return arr3[tidx, :, :]

    raise ValueError(f"Variable '{varname}' has unsupported dims={dims}; only 2D or 3D supported.")


def subset_by_bbox(lat: np.ndarray, lon: np.ndarray, field: np.ndarray,
                   bbox: Tuple[float, float, float, float]):
    lonmin, lonmax, latmin, latmax = bbox
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)

    # 兼容 0..360 经度输入的 bbox。
    lon_work = lon.copy()
    if lon_work.min() >= 0 and lonmin < 0:
        lon_work = ((lon_work + 180) % 360) - 180

    lat_mask = (lat >= latmin) & (lat <= latmax)
    lon_mask = (lon_work >= lonmin) & (lon_work <= lonmax)

    if not lat_mask.any() or not lon_mask.any():
        raise ValueError(
            f"No points in bbox lon[{lonmin},{lonmax}] lat[{latmin},{latmax}]. "
            f"Available lon[{lon.min()},{lon.max()}], lat[{lat.min()},{lat.max()}]."
        )

    return lat[lat_mask], lon[lon_mask], field[np.ix_(lat_mask, lon_mask)]


def mask_zero(field: np.ndarray, eps: float) -> np.ndarray:
    arr = np.asarray(field, dtype=float).copy()
    arr[np.abs(arr) <= eps] = np.nan
    return arr


def choose_norm_and_cmap(field: np.ndarray, cmap_arg: Optional[str] = None, symmetric: bool = False):
    data = np.asarray(field, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return Normalize(vmin=-1, vmax=1), (cmap_arg or "viridis")

    vmin = float(np.nanmin(data))
    vmax = float(np.nanmax(data))
    
    if symmetric:
        m = max(abs(vmin), abs(vmax))
        if m == 0:
            m = 1e-12
        return TwoSlopeNorm(vmin=-m, vcenter=0.0, vmax=m), (cmap_arg or "RdYlGn_r")
    
    if vmin < 0 < vmax:
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax), (cmap_arg or "RdYlGn_r")

    if vmax <= 0:
        return Normalize(vmin=vmin, vmax=vmax), (cmap_arg or "Reds_r")
    
    return Normalize(vmin=vmin, vmax=vmax), (cmap_arg or "Greens")


def is_carbon_like(varname: str) -> bool:
    lname = varname.lower()
    return not any(k in lname for k in NON_CARBON_KEYWORDS)


def maybe_convert_units(field: np.ndarray, varname: str, is_diff: bool, no_scale: bool):
    if no_scale or not is_carbon_like(varname):
        return field, "native units"
    return field / 1e6, ("Tg C / yr" if is_diff else "Tg C")


def plot_field(lat, lon, field, title: str, unit: str, use_cartopy: bool,
               cmap: Optional[str], symmetric: bool):
    norm, cmap_name = choose_norm_and_cmap(field, cmap_arg=cmap, symmetric=symmetric)
    field = np.where(np.isfinite(field), field, np.nan)
    extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]
    lat_desc = bool(lat[0] > lat[-1]) if len(lat) > 1 else False

    if use_cartopy:
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature

            fig = plt.figure(figsize=(11, 5.5))
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6)
            mesh = ax.pcolormesh(lon, lat, field, transform=ccrs.PlateCarree(),
                                 cmap=cmap_name, norm=norm, shading="auto")
            cb = plt.colorbar(mesh, ax=ax, orientation="vertical", pad=0.02, fraction=0.04)
            cb.set_label(unit)

            if isinstance(norm, TwoSlopeNorm):
                m = max(abs(norm.vmin), abs(norm.vmax))
                cb.set_ticks([-m, -m/2, 0, m/2, m])
            ax.set_title(title)
            plt.tight_layout()
            return fig, ax
        except Exception as e:
            print(f"[WARN] Cartopy unavailable or failed; fallback to matplotlib. Reason: {e}", file=sys.stderr)

    fig = plt.figure(figsize=(11, 5.5))
    ax = plt.gca()
    origin = "upper" if lat_desc else "lower"
    im = ax.imshow(field, extent=extent, origin=origin, cmap=cmap_name, norm=norm, aspect="auto")
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(unit)

    if isinstance(norm, TwoSlopeNorm):
        m = max(abs(norm.vmin), abs(norm.vmax))
        cb.set_ticks([-m, -m/2, 0, m/2, m])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)
    plt.tight_layout()
    return fig, ax


def parse_bbox_arg(bbox_args):
    if bbox_args is None or len(bbox_args) == 0:
        return None
    if len(bbox_args) == 1:
        key = bbox_args[0].lower()
        if key not in BBOX_PRESETS:
            raise ValueError(f"Unknown bbox preset '{key}'. Available: {list(BBOX_PRESETS.keys())}")
        return BBOX_PRESETS[key]
    if len(bbox_args) == 4:
        return tuple(map(float, bbox_args))
    raise ValueError("--bbox must be a preset name or 4 numbers: lonmin lonmax latmin latmax")


def describe_var(ds: Dataset, varname: str):
    v = ds.variables[varname]
    print(f"\nVariable: {varname}")
    print(f"  dims  = {v.dimensions}")
    print(f"  shape = {v.shape}")
    for attr in ["units", "long_name", "standard_name"]:
        if hasattr(v, attr):
            print(f"  {attr} = {getattr(v, attr)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Viewer for bookkeeping NetCDF outputs.")
    parser.add_argument("ncfile", help="Path to NetCDF file")
    parser.add_argument("--list", action="store_true", help="List variables and exit")
    parser.add_argument("--var", dest="varname", default=None, help="Variable to plot")
    parser.add_argument("--tidx", type=int, default=None, help="Time index to plot")
    parser.add_argument("--year", type=int, default=None, help="Actual year to plot, e.g. 1900")
    parser.add_argument("--diff-index", nargs=2, type=int, metavar=("N", "M"),
                        help="Plot annualized difference: (var[index N] - var[index M]) / year_diff")
    parser.add_argument("--diff-year", nargs=2, type=int, metavar=("YEAR_N", "YEAR_M"),
                        help="Plot annualized difference: (var[year N] - var[year M]) / year_diff")
    parser.add_argument("--bbox", nargs="*", default=None,
                        help="Preset bbox name or 4 numbers lonmin lonmax latmin latmax")
    parser.add_argument("--save", default=None, help="Save figure path")
    parser.add_argument("--dpi", type=int, default=200, help="Save dpi")
    parser.add_argument("--no-cartopy", action="store_true", help="Disable cartopy")
    parser.add_argument("--no-scale", action="store_true", help="Do not divide carbon-like variables by 1e9")
    parser.add_argument("--mask-zero", action="store_true", help="Display near-zero values as white/transparent")
    parser.add_argument("--zero-eps", type=float, default=0.0, help="Threshold for --mask-zero")
    parser.add_argument("--cmap", default=None, help="Override colormap name")
    parser.add_argument("--symmetric", action="store_true", help="Use symmetric colorbar around zero")
    args = parser.parse_args()

    with Dataset(args.ncfile, "r") as ds:
        if args.list:
            print("Variables:")
            for name in list_vars(ds):
                v = ds.variables[name]
                units = getattr(v, "units", "")
                print(f"  - {name:32s} dims={v.dimensions} shape={v.shape} units={units}")
            tvals = get_time_values(ds)
            if tvals is not None:
                print(f"\nTime range: {tvals.min()} .. {tvals.max()}  (n={len(tvals)})")
            return 0

        lat_name, lon_name = find_lat_lon(ds)
        lat = np.asarray(ds.variables[lat_name][:], dtype=float).squeeze()
        lon = np.asarray(ds.variables[lon_name][:], dtype=float).squeeze()
        if lat.ndim != 1 or lon.ndim != 1:
            raise ValueError(f"Expected 1D lat/lon. Got lat.ndim={lat.ndim}, lon.ndim={lon.ndim}")

        if args.varname is None:
            candidates = list_vars(ds)
            if not candidates:
                raise ValueError("No plottable variables found.")
            args.varname = candidates[0]
            print(f"[INFO] --var not provided, using '{args.varname}'")

        describe_var(ds, args.varname)
        is_diff = False

        if args.diff_year is not None:
            year_n, year_m = args.diff_year
            idx_n = year_to_index(ds, year_n)
            idx_m = year_to_index(ds, year_m)
            year_diff = year_n - year_m
            if year_diff == 0:
                raise ValueError("Year difference cannot be zero.")
            field = (get_2d_field(ds, args.varname, idx_n) - get_2d_field(ds, args.varname, idx_m)) / abs(year_diff)
            title = f"{args.varname}: ({year_n} - {year_m}) / {abs(year_diff)} yr"
            is_diff = True

        elif args.diff_index is not None:
            n, m = args.diff_index
            tvals = get_time_values(ds)
            if tvals is not None and n < len(tvals) and m < len(tvals):
                year_n, year_m = int(tvals[n]), int(tvals[m])
                year_diff = year_n - year_m
                title = f"{args.varname}: ({year_n} - {year_m}) / {abs(year_diff)} yr [index {n}-{m}]"
            else:
                year_diff = n - m
                title = f"{args.varname}: (index {n} - {m}) / {abs(year_diff)} step"
            if year_diff == 0:
                raise ValueError("Index difference leads to zero interval.")
            field = (get_2d_field(ds, args.varname, n) - get_2d_field(ds, args.varname, m)) / abs(year_diff)
            is_diff = True

        elif args.year is not None:
            idx = year_to_index(ds, args.year)
            field = get_2d_field(ds, args.varname, idx)
            title = f"{args.varname}: year {args.year} (index {idx})"

        else:
            tidx = args.tidx if args.tidx is not None else 0
            field = get_2d_field(ds, args.varname, tidx)
            tvals = get_time_values(ds)
            if tvals is not None and tidx < len(tvals):
                title = f"{args.varname}: year {int(tvals[tidx])} (index {tidx})"
            else:
                title = f"{args.varname}: 2D field" if find_time(ds) is None else f"{args.varname}: index {tidx}"

        field, unit = maybe_convert_units(field, args.varname, is_diff, args.no_scale)

        if args.mask_zero:
            field = mask_zero(field, args.zero_eps)
            title += " [zero masked]"

        bbox = parse_bbox_arg(args.bbox)
        if bbox is not None:
            lat, lon, field = subset_by_bbox(lat, lon, field, bbox)
            title += f" [bbox={args.bbox}]"

        title += f" [{unit}]"
        fig, _ = plot_field(lat, lon, field, title, unit, use_cartopy=(not args.no_cartopy),
                            cmap=args.cmap, symmetric=args.symmetric)

        if args.save:
            fig.savefig(args.save, dpi=args.dpi, bbox_inches="tight")
            print(f"[OK] Saved to: {args.save}")
        else:
            plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
