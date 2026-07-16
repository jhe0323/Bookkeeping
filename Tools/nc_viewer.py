#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nc_viewer.py
适配当前 bookkeeping summary 输出文件结构的 NetCDF 查看脚本

当前支持的典型变量：
- biomass_total
- soil_total
- P1
- P10
- P100
- atmosphere
- unmet_harvest
- done

功能：
1. 查看某变量在某个时间步/年份的空间分布
2. 查看两个时间步/年份之间的差值图
3. 支持 bbox 子区
4. 支持保存图片
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Tuple, Optional

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

BBOX_PRESETS = {
    "global": (-180.0, 180.0, -90.0, 90.0),
    "china": (73.0, 135.0, 18.0, 54.0),
    "loess": (100.0, 114.0, 33.0, 41.5),
    "xinjiang": (73.0, 97.5, 34.0, 50.5),
}


def _find_first_var(ds: Dataset, names: List[str]) -> Optional[str]:
    for n in names:
        if n in ds.variables:
            return n
    return None


def _find_by_attr(ds: Dataset, attr: str, value: str) -> Optional[str]:
    for k, v in ds.variables.items():
        if hasattr(v, attr) and str(getattr(v, attr)).lower() == value.lower():
            return k
    return None


def find_lat_lon(ds: Dataset) -> Tuple[str, str]:
    lat_name = _find_first_var(ds, LAT_CANDIDATES)
    lon_name = _find_first_var(ds, LON_CANDIDATES)

    if lat_name is None:
        lat_name = _find_by_attr(ds, "standard_name", "latitude")
    if lon_name is None:
        lon_name = _find_by_attr(ds, "standard_name", "longitude")

    if lat_name is None or lon_name is None:
        raise ValueError("Could not find lat/lon variables.")
    return lat_name, lon_name


def find_time(ds: Dataset) -> Optional[str]:
    t = _find_first_var(ds, TIME_CANDIDATES)
    if t is None:
        t = _find_by_attr(ds, "standard_name", "time")
    return t


def list_vars(ds: Dataset) -> List[str]:
    # 不排除 done，因为它也可能想看
    coord_like = set(LAT_CANDIDATES + LON_CANDIDATES + TIME_CANDIDATES)
    out = []
    for name in ds.variables:
        if name in coord_like:
            continue
        out.append(name)
    return sorted(out)


def get_time_values(ds: Dataset) -> Optional[np.ndarray]:
    time_name = find_time(ds)
    if time_name is None or time_name not in ds.variables:
        return None
    return np.asarray(ds.variables[time_name][:]).astype(int)


def year_to_index(ds: Dataset, year: int) -> int:
    tvals = get_time_values(ds)
    if tvals is None:
        raise ValueError("No time variable found in file, cannot use --year or --diff-year.")
    matches = np.where(tvals == year)[0]
    if len(matches) == 0:
        raise ValueError(
            f"Year {year} not found in time variable. "
            f"Available range: {tvals.min()} .. {tvals.max()}"
        )
    return int(matches[0])


def get_2d_field(ds: Dataset, varname: str, tidx: Optional[int] = None) -> np.ndarray:
    if varname not in ds.variables:
        raise ValueError(f"Variable '{varname}' not found.")

    v = ds.variables[varname]
    dims = list(v.dimensions)

    lat_name, lon_name = find_lat_lon(ds)
    time_name = find_time(ds)

    # 典型情况：
    # 3D: (time, lat, lon)
    # 2D: (lat, lon)
    if len(dims) == 3:
        if time_name is None or dims[0] != time_name:
            raise ValueError(
                f"Variable '{varname}' is 3D but not in expected (time, lat, lon) format. dims={dims}"
            )
        if tidx is None:
            tidx = 0
        arr = v[tidx, :, :]

    elif len(dims) == 2:
        arr = v[:, :]

    else:
        raise ValueError(
            f"Variable '{varname}' has unsupported dims={dims}. "
            "This viewer currently supports 2D (lat,lon) or 3D (time,lat,lon)."
        )

    arr = np.array(arr, dtype=float)
    return arr


def subset_by_bbox(lat: np.ndarray, lon: np.ndarray, field: np.ndarray,
                   bbox: Tuple[float, float, float, float]):
    lonmin, lonmax, latmin, latmax = bbox

    lat = np.asarray(lat).astype(float)
    lon = np.asarray(lon).astype(float)

    lat_mask = (lat >= latmin) & (lat <= latmax)
    lon_mask = (lon >= lonmin) & (lon <= lonmax)

    if not lat_mask.any() or not lon_mask.any():
        raise ValueError(
            f"No points in bbox lon[{lonmin},{lonmax}] lat[{latmin},{latmax}]. "
            f"Available lon[{lon.min()},{lon.max()}], lat[{lat.min()},{lat.max()}]."
        )

    lat_sub = lat[lat_mask]
    lon_sub = lon[lon_mask]
    field_sub = field[np.ix_(lat_mask, lon_mask)]
    return lat_sub, lon_sub, field_sub


def choose_norm(field: np.ndarray):
    data = np.asarray(field, dtype=float)
    data = data[np.isfinite(data)]

    if data.size == 0:
        return Normalize(vmin=-1, vmax=1), "viridis"

    vmin = float(np.nanmin(data))
    vmax = float(np.nanmax(data))

    if vmin < 0 < vmax:
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax), "RdYlGn_r"
    elif vmax <= 0:
        if vmin == vmax:
            vmax = vmin + 1e-12
        return Normalize(vmin=vmin, vmax=vmax), "Reds_r"
    else:
        if vmin == vmax:
            vmax = vmin + 1e-12
        return Normalize(vmin=vmin, vmax=vmax), "Greens"


def plot_field(lat, lon, field, title: str, use_cartopy: bool = True):
    norm, cmap = choose_norm(field)

    field = np.where(np.isfinite(field), field, np.nan)
    extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]
    lat_desc = bool(lat[0] > lat[-1])

    if use_cartopy:
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature

            fig = plt.figure(figsize=(11, 5.5))
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6)
            #ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.4)

            mesh = ax.pcolormesh(
                lon, lat, field,
                transform=ccrs.PlateCarree(),
                cmap=cmap, norm=norm, shading="auto"
            )
            plt.colorbar(mesh, ax=ax, orientation="vertical", pad=0.02, fraction=0.04)
            ax.set_title(title)
            plt.tight_layout()
            return fig, ax

        except Exception as e:
            print(f"[WARN] Cartopy unavailable or failed, fallback to matplotlib. Reason: {e}",
                  file=sys.stderr)

    fig = plt.figure(figsize=(11, 5.5))
    ax = plt.gca()
    origin = "upper" if lat_desc else "lower"
    im = ax.imshow(field, extent=extent, origin=origin, cmap=cmap, norm=norm, aspect="auto")
    plt.colorbar(im, ax=ax)
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


def main():
    parser = argparse.ArgumentParser(description="Viewer for current bookkeeping summary NetCDF outputs.")
    parser.add_argument("ncfile", help="Path to NetCDF file")
    parser.add_argument("--list", action="store_true", help="List variables and exit")

    parser.add_argument("--var", dest="varname", default=None, help="Variable to plot")
    parser.add_argument("--tidx", type=int, default=None, help="Time index to plot")
    parser.add_argument("--year", type=int, default=None, help="Actual year to plot, e.g. 1900")

    parser.add_argument("--diff-index", nargs=2, type=int, metavar=("N", "M"),
                        help="Plot difference var[index N] - var[index M]")
    parser.add_argument("--diff-year", nargs=2, type=int, metavar=("YEAR_N", "YEAR_M"),
                        help="Plot difference var[year N] - var[year M]")

    parser.add_argument("--bbox", nargs="*", default=None,
                        help="Preset bbox name or 4 numbers lonmin lonmax latmin latmax")
    parser.add_argument("--save", default=None, help="Save figure path")
    parser.add_argument("--dpi", type=int, default=200, help="Save dpi")
    parser.add_argument("--no-cartopy", action="store_true", help="Disable cartopy")
    args = parser.parse_args()

    with Dataset(args.ncfile, "r") as ds:
        if args.list:
            print("Variables:")
            for v in list_vars(ds):
                vv = ds.variables[v]
                print(f"  - {v:20s} dims={vv.dimensions} shape={vv.shape}")
            tvals = get_time_values(ds)
            if tvals is not None:
                print(f"\nTime range: {tvals.min()} .. {tvals.max()}  (n={len(tvals)})")
            return 0

        lat_name, lon_name = find_lat_lon(ds)
        lat = np.asarray(ds.variables[lat_name][:]).astype(float).squeeze()
        lon = np.asarray(ds.variables[lon_name][:]).astype(float).squeeze()

        if lat.ndim != 1 or lon.ndim != 1:
            raise ValueError(f"Expected 1D lat/lon. Got lat.ndim={lat.ndim}, lon.ndim={lon.ndim}")

        if args.varname is None:
            candidates = list_vars(ds)
            if not candidates:
                raise ValueError("No variables found.")
            args.varname = candidates[0]
            print(f"[INFO] --var not provided, using '{args.varname}'")

        describe_var(ds, args.varname)

        # -------------------------
        # 1) 差值图：按年份
        # -------------------------
        if args.diff_year is not None:
            year_n, year_m = args.diff_year
            idx_n = year_to_index(ds, year_n)
            idx_m = year_to_index(ds, year_m)

            f_n = get_2d_field(ds, args.varname, tidx=idx_n)
            f_m = get_2d_field(ds, args.varname, tidx=idx_m)
            year_diff = abs(year_n - year_m)
            if year_diff == 0:
                raise ValueError("Year difference cannot be zero.")

            field = (f_n - f_m) / year_diff
            title = f"{args.varname}: ({year_n} - {year_m}) / {year_diff} yr"

        # -------------------------
        # 2) 差值图：按索引
        # -------------------------
        elif args.diff_index is not None:
            n, m = args.diff_index
            f_n = get_2d_field(ds, args.varname, tidx=n)
            f_m = get_2d_field(ds, args.varname, tidx=m)
            tvals = get_time_values(ds)

            if tvals is not None and n < len(tvals) and m < len(tvals):
                year_diff = int(tvals[n] - tvals[m])
            else:
                year_diff = n - m

            if year_diff == 0:
                raise ValueError("Index difference leads to zero time interval.")

            field = (f_n - f_m) / year_diff

            tvals = get_time_values(ds)
            if tvals is not None and n < len(tvals) and m < len(tvals):
                title = f"{args.varname}: {tvals[n]} - {tvals[m]} (index {n} - {m})"
            else:
                title = f"{args.varname}: index {n} - index {m}"

            tvals = get_time_values(ds)
            if tvals is not None and n < len(tvals) and m < len(tvals):
                title = f"{args.varname}: {tvals[n]} - {tvals[m]} (index {n} - {m})"
            else:
                title = f"{args.varname}: index {n} - index {m}"

        # -------------------------
        # 3) 单时刻图：按年份
        # -------------------------
        elif args.year is not None:
            idx = year_to_index(ds, args.year)
            field = get_2d_field(ds, args.varname, tidx=idx)
            title = f"{args.varname}: year {args.year} (index {idx})"

        # -------------------------
        # 4) 单时刻图：按索引
        # -------------------------
        else:
            tidx = args.tidx if args.tidx is not None else 0
            field = get_2d_field(ds, args.varname, tidx=tidx)

            tvals = get_time_values(ds)
            if tvals is not None and tidx < len(tvals):
                title = f"{args.varname}: year {tvals[tidx]} (index {tidx})"
            else:
                title = f"{args.varname}: index {tidx}"

        if args.varname in ["atmosphere", "biomass_total", "soil_total", "P1", "P10", "P100"]:
            field = field / 1e9

            # 判断是不是“差分图”
            if args.diff_year is not None or args.diff_index is not None:
                unit = "Gt C / yr"
            else:
                unit = "Gt C"

            title += f" [{unit}]"
        
        bbox = parse_bbox_arg(args.bbox)
        if bbox is not None:
            lat, lon, field = subset_by_bbox(lat, lon, field, bbox)
            title += " [Gt C/yr]"

        fig, ax = plot_field(lat, lon, field, title=title, use_cartopy=(not args.no_cartopy))

        if args.save:
            fig.savefig(args.save, dpi=args.dpi, bbox_inches="tight")
            print(f"[OK] Saved to: {args.save}")
        else:
            plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())