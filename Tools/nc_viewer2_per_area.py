#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import sys
from typing import List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import matplotlib.colors as mcolors
from netCDF4 import Dataset

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
        raise ValueError("No time variable found.")

    matches = np.where(tvals == year)[0]

    if len(matches) == 0:
        raise ValueError(
            f"Year {year} not found. Available range: {tvals.min()} .. {tvals.max()}"
        )

    return int(matches[0])


def get_2d_field(ds: Dataset, varname: str, tidx: Optional[int] = None) -> np.ndarray:
    if varname not in ds.variables:
        raise ValueError(f"Variable '{varname}' not found.")

    v = ds.variables[varname]
    dims = list(v.dimensions)

    lat_name, lon_name = find_lat_lon(ds)
    time_name = find_time(ds)

    if len(dims) == 2:
        arr = np.asarray(v[:], dtype=float)
        axes = [dims.index(lat_name), dims.index(lon_name)]
        return np.transpose(arr, axes=axes)

    if len(dims) == 3:
        if tidx is None:
            tidx = 0

        arr3 = np.asarray(v[:], dtype=float)
        axes = [dims.index(time_name), dims.index(lat_name), dims.index(lon_name)]
        arr3 = np.transpose(arr3, axes=axes)

        return arr3[tidx, :, :]

    raise ValueError(f"Unsupported dimensions: {dims}")


def parse_bbox_arg(bbox):
    if bbox is None:
        return None

    if isinstance(bbox, str):
        key = bbox.lower()
        if key not in BBOX_PRESETS:
            raise ValueError(f"Unknown bbox preset: {bbox}")
        return BBOX_PRESETS[key]

    if len(bbox) == 4:
        return tuple(map(float, bbox))

    raise ValueError("bbox must be None, preset name, or 4 numbers.")


def subset_by_bbox(lat, lon, field, bbox):
    lonmin, lonmax, latmin, latmax = bbox

    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)

    lon_work = lon.copy()

    if lon_work.min() >= 0 and lonmin < 0:
        lon_work = ((lon_work + 180) % 360) - 180

    lat_mask = (lat >= latmin) & (lat <= latmax)
    lon_mask = (lon_work >= lonmin) & (lon_work <= lonmax)

    if not lat_mask.any() or not lon_mask.any():
        raise ValueError("No points found in selected bbox.")

    return lat[lat_mask], lon[lon_mask], field[np.ix_(lat_mask, lon_mask)]


def mask_zero(field, eps=0.0):
    arr = np.asarray(field, dtype=float).copy()
    arr[np.abs(arr) <= eps] = np.nan
    return arr


def _coord_edges_from_centers(x, lower, upper):
    x = np.asarray(x, dtype=float)

    if x.size == 1:
        dx = upper - lower
        edges = np.array([x[0] - dx / 2.0, x[0] + dx / 2.0])
    else:
        mid = (x[:-1] + x[1:]) / 2.0
        first = x[0] - (x[1] - x[0]) / 2.0
        last = x[-1] + (x[-1] - x[-2]) / 2.0
        edges = np.concatenate([[first], mid, [last]])

    return np.clip(edges, lower, upper)


def grid_cell_area_m2(lat, lon, radius=6371000.0):
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)

    lat_edges = _coord_edges_from_centers(lat, -90.0, 90.0)
    lon_edges = _coord_edges_from_centers(
        lon,
        lon.min() - 180.0,
        lon.max() + 180.0
    )

    dlon = np.abs(np.diff(np.deg2rad(lon_edges)))
    sin_lat = np.sin(np.deg2rad(lat_edges))
    d_sin_lat = np.abs(np.diff(sin_lat))

    area = radius ** 2 * d_sin_lat[:, None] * dlon[None, :]

    return area


def is_carbon_like(varname: str) -> bool:
    lname = varname.lower()
    return not any(k in lname for k in NON_CARBON_KEYWORDS)


def maybe_convert_units(field, varname, is_diff, no_scale=False):
    if no_scale or not is_carbon_like(varname):
        return field, "native units"

    return field / 1e6, ("Tg C / yr" if is_diff else "Tg C")


def convert_to_area_flux(field, area_m2, varname, is_diff):
    if not is_carbon_like(varname):
        raise ValueError("per_area should only be used for carbon-like variables.")

    out = field * 1e6 / area_m2

    unit = "g C m$^{-2}$ yr$^{-1}$" if is_diff else "g C m$^{-2}$"

    return out, unit


def choose_norm_and_cmap(field, cmap_arg=None, symmetric=False):
    data = np.asarray(field, dtype=float)
    data = data[np.isfinite(data)]

    if data.size == 0:
        return Normalize(vmin=-1, vmax=1), (cmap_arg or "viridis")

    vmin = float(np.nanmin(data))
    vmax = float(np.nanmax(data))

    if np.isclose(vmin, vmax):
        eps = abs(vmin) * 0.01 if vmin != 0 else 1e-12
        vmin -= eps
        vmax += eps

    if symmetric:
        m = max(abs(vmin), abs(vmax))
        return TwoSlopeNorm(vmin=-m, vcenter=0.0, vmax=m), (cmap_arg or "RdYlGn_r")

    if vmin < 0 < vmax:
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax), (cmap_arg or "RdYlGn_r")

    if vmax <= 0:
        return Normalize(vmin=vmin, vmax=vmax), (cmap_arg or "Greens_r")

    return Normalize(vmin=vmin, vmax=vmax), (cmap_arg or "Reds")

def get_gcb_colormap(levels=None):
    """
    GCB固定分级色标：
    每个网格值落在哪个区间，就显示对应颜色。
    0 只是分界线，不单独占颜色。
    """

    if levels is None:
        fix_par = 0.2
        levels = fix_par*np.array([
            -400, -200, -100, -50, -25, -12.5,
            0,
            12.5, 25, 50, 100, 200, 400
        ], dtype=float)
    else:
        levels = np.asarray(levels, dtype=float)

    colors = [
        "#053061",  # -400 ~ -200
        "#2166AC",  # -200 ~ -100
        "#4393C3",  # -100 ~ -50
        "#92C5DE",  # -50 ~ -25
        "#D1E5F0",  # -25 ~ -12.5
        "#EAF3F8",  # -12.5 ~ 0

        "#FCEBE2",  # 0 ~ 12.5
        "#FDD0B6",  # 12.5 ~ 25
        "#F4A582",  # 25 ~ 50
        "#D6604D",  # 50 ~ 100
        "#B2182B",  # 100 ~ 200
        "#67001F",  # 200 ~ 400
    ]

    cmap = mcolors.ListedColormap(colors)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])

    norm = mcolors.BoundaryNorm(levels, cmap.N)

    return cmap, norm, levels


def plot_field(
    lat, lon, field, title, unit,
    use_cartopy=False,
    cmap=None,
    symmetric=False
):
    gcb_cmap, norm, levels = get_gcb_colormap()

    field = np.where(np.isfinite(field), field, np.nan)

    extent = [
        float(lon.min()),
        float(lon.max()),
        float(lat.min()),
        float(lat.max())
    ]

    lat_desc = bool(lat[0] > lat[-1]) if len(lat) > 1 else False

    fig = plt.figure(figsize=(12, 5.8))
    ax = plt.gca()

    origin = "upper" if lat_desc else "lower"

    im = ax.imshow(
        field,
        extent=extent,
        origin=origin,
        cmap=gcb_cmap,
        norm=norm,
        aspect="auto"
    )

    cb = plt.colorbar(
        im,
        ax=ax,
        orientation="horizontal",
        pad=0.12,
        fraction=0.055,
        shrink=0.95,
        ticks=levels,
        spacing="uniform",
        extend="both"
    )

    cb.set_label(unit)

    # 防止 tick 挤在一起
    cb.ax.tick_params(labelsize=9)
    cb.ax.set_xticklabels([f"{x:g}" for x in levels])

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)

    plt.tight_layout()

    return fig, ax
# =========================
# 1. 显式设置文件路径
# =========================
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


ncfile = r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\summary_025deg.global_time_gt1100.nc"
# ncfile = r"D:\data\atmosphere_1960_2020.nc"

# =========================
# 2. 选择变量
# =========================

varname = "Net_Emissions"

# =========================
# 3. 时间设置
# =========================

# 方式 A：看单一年份
year = None

# 方式 B：看两个年份差值，单位会自动变成年均
# 如果使用 diff_year，就把 year 设为 None
diff_year = None
# diff_year = (2020, 1960)

# 方式 C：看时间索引差值
diff_index = (70, 60)
# diff_index = (1170, 1110)

# =========================
# 4. 区域设置
# =========================

bbox = None
# bbox = "china"
# bbox = "loess"
# bbox = (73, 135, 18, 54)   # lonmin, lonmax, latmin, latmax

# =========================
# 5. 单位和显示设置
# =========================

per_area = True      # True: 转成 g C m-2 yr-1
no_scale = False      # False: 碳变量自动 /1e6 显示为 Tg C 或 Tg C/yr
mask_zero_flag = False # True: 0 值显示为空白/白色
zero_eps = 0.0

use_cartopy = False   # Notebook 中建议先 False
cmap = None           # None 自动选择
symmetric = False     # True 强制 colorbar 以 0 为中心对称

save_fig = False
save_path = r"/mnt/beegfs/product/lulc0120/outputs/output.png"
dpi = 300




ds = Dataset(ncfile, "r")

print("Available variables:")
for name in list_vars(ds):
    v = ds.variables[name]
    units = getattr(v, "units", "")
    print(f"{name:32s} dims={v.dimensions} shape={v.shape} units={units}")

tvals = get_time_values(ds)
if tvals is not None:
    print()
    print(f"Time range: {tvals.min()} .. {tvals.max()}  n={len(tvals)}")

lat_name, lon_name = find_lat_lon(ds)

lat = np.asarray(ds.variables[lat_name][:], dtype=float).squeeze()
lon = np.asarray(ds.variables[lon_name][:], dtype=float).squeeze()

is_diff = False

if diff_year is not None:
    year_n, year_m = diff_year

    idx_n = year_to_index(ds, year_n)
    idx_m = year_to_index(ds, year_m)

    field = (
        get_2d_field(ds, varname, idx_n)
        - get_2d_field(ds, varname, idx_m)
    ) / abs(year_n - year_m)

    title = f"{varname}: ({year_n} - {year_m}) / {abs(year_n - year_m)} yr"
    is_diff = True

elif diff_index is not None:
    n, m = diff_index

    tvals = get_time_values(ds)

    field = (
        get_2d_field(ds, varname, n)
        - get_2d_field(ds, varname, m)
    ) / abs(n - m)

    if tvals is not None:
        title = f"{varname}: ({int(tvals[n])} - {int(tvals[m])}) / {abs(int(tvals[n]) - int(tvals[m]))} yr"
    else:
        title = f"{varname}: index ({n} - {m}) / {abs(n - m)} step"

    is_diff = True

elif year is not None:
    idx = year_to_index(ds, year)
    field = get_2d_field(ds, varname, idx)

    title = f"{varname}: year {year} index {idx}"

else:
    field = get_2d_field(ds, varname, tidx=0)
    title = f"{varname}: 2D field"

print(title)
print("field shape:", field.shape)
print("raw min/max:", np.nanmin(field), np.nanmax(field))

area_m2 = grid_cell_area_m2(lat, lon)

bbox_value = parse_bbox_arg(bbox)

if bbox_value is not None:
    lat, lon, field = subset_by_bbox(lat, lon, field, bbox_value)
    _, _, area_m2 = subset_by_bbox(
        np.asarray(ds.variables[lat_name][:], dtype=float).squeeze(),
        np.asarray(ds.variables[lon_name][:], dtype=float).squeeze(),
        area_m2,
        bbox_value
    )

    title += f" [bbox={bbox}]"

if per_area:
    field, unit = convert_to_area_flux(field, area_m2, varname, is_diff)
    title += " [per area]"
else:
    field, unit = maybe_convert_units(field, varname, is_diff, no_scale=no_scale)

if mask_zero_flag:
    field = mask_zero(field, zero_eps)
    title += " [zero masked]"

title += f" [{unit}]"

print("plot min/max:", np.nanmin(field), np.nanmax(field))
print("unit:", unit)

fig, ax = plot_field(
    lat=lat,
    lon=lon,
    field=field,
    title=title,
    unit=unit,
    use_cartopy=use_cartopy
)

plt.show()

if save_fig:
    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    print(f"Saved to: {save_path}")