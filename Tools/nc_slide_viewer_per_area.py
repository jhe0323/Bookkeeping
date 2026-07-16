#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
from netCDF4 import Dataset
from typing import List, Optional, Sequence, Tuple


# ============================================================
# 1. 用户设置
# ============================================================

ncfile = r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\summary_025deg.global_time_ge1000_7.9.nc"

varname = "Gross_Sinks"
start_year = 1961
end_year = 2020
time_offset = 850
bbox = None
# bbox = "china"
# bbox = "loess"
# bbox = (73, 135, 18, 54)

per_area = True          # True: Mg C yr-1 -> g C m-2 yr-1
mask_zero_flag = False
zero_eps = 1.0

use_cartopy = False      # 没安装 cartopy 就用 False
legend_percentile = None  #(1, 99)
# legend_percentile = (1, 99)
# legend_percentile = (0.5, 99.5)
save_fig = False
save_path = r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\Net_Emissions_mean.png"
dpi = 300


# ============================================================
# 2. 基础参数
# ============================================================

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


# ============================================================
# 3. 工具函数
# ============================================================

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


def get_2d_field(ds: Dataset, varname: str, tidx: int = 0) -> np.ndarray:
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
        arr3 = np.asarray(v[:], dtype=float)
        axes = [dims.index(time_name), dims.index(lat_name), dims.index(lon_name)]
        arr3 = np.transpose(arr3, axes=axes)
        return arr3[tidx, :, :]

    raise ValueError(f"Unsupported dimensions: {dims}")

def get_mean_2d_field(
    ds: Dataset,
    varname: str,
    start_year: Optional[int] = None,
    end_year: Optional[int] = None,
    time_offset: int = 0,
) -> Tuple[np.ndarray, Optional[int], Optional[int], int]:
    """
    读取指定年份范围并计算多年平均。

    Parameters
    ----------
    ds
        NetCDF数据集。
    varname
        变量名。
    start_year, end_year
        公历年份范围，包含起止年份。
        如果均为None，则对全部时间层平均。
    time_offset
        time变量到公历年份的偏移量。
        例如LUH2中time=0对应850年，则设为850。

    Returns
    -------
    mean_field
        多年平均二维场。
    actual_start
        实际使用的起始年份。
    actual_end
        实际使用的结束年份。
    n_years
        使用的时间层数量。
    """
    if varname not in ds.variables:
        raise ValueError(f"Variable '{varname}' not found.")

    v = ds.variables[varname]
    dims = list(v.dimensions)

    lat_name, lon_name = find_lat_lon(ds)
    time_name = find_time(ds)

    # 二维变量无需多年平均
    if len(dims) == 2:
        arr = np.asarray(v[:], dtype=float)
        axes = [dims.index(lat_name), dims.index(lon_name)]
        arr = np.transpose(arr, axes=axes)
        return arr, None, None, 1

    if len(dims) != 3 or time_name not in dims:
        raise ValueError(
            f"Variable '{varname}' must have dimensions "
            f"(time, lat, lon) or (lat, lon). Current dimensions: {dims}"
        )

    time_values = np.asarray(ds.variables[time_name][:], dtype=int)
    calendar_years = time_values + time_offset

    if start_year is None:
        start_year = int(calendar_years.min())

    if end_year is None:
        end_year = int(calendar_years.max())

    if start_year > end_year:
        raise ValueError("start_year cannot be greater than end_year.")

    time_mask = (calendar_years >= start_year) & (calendar_years <= end_year)
    time_indices = np.where(time_mask)[0]

    if time_indices.size == 0:
        raise ValueError(
            f"No time steps found for {start_year}-{end_year}. "
            f"Available calendar years: "
            f"{calendar_years.min()}-{calendar_years.max()}"
        )

    # 只读取需要的时间段
    time_axis = dims.index(time_name)
    lat_axis = dims.index(lat_name)
    lon_axis = dims.index(lon_name)

    slices = [slice(None)] * 3
    slices[time_axis] = slice(time_indices[0], time_indices[-1] + 1)

    arr3 = np.asarray(v[tuple(slices)], dtype=float)

    # 转为 time, lat, lon 顺序
    arr3 = np.transpose(
        arr3,
        axes=[time_axis, lat_axis, lon_axis]
    )

    # 处理NetCDF掩膜值和无效值
    arr3 = np.where(np.isfinite(arr3), arr3, np.nan)

    mean_field = np.nanmean(arr3, axis=0)

    actual_start = int(calendar_years[time_indices[0]])
    actual_end = int(calendar_years[time_indices[-1]])

    return mean_field, actual_start, actual_end, len(time_indices)
    
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

    return radius ** 2 * d_sin_lat[:, None] * dlon[None, :]


def is_carbon_like(varname: str) -> bool:
    lname = varname.lower()
    return not any(k in lname for k in NON_CARBON_KEYWORDS)


def convert_to_area_flux(field, area_m2, varname):
    if not is_carbon_like(varname):
        raise ValueError("per_area should only be used for carbon-like variables.")

    # 模型输出单位：Mg C yr-1
    # 1 Mg C = 1e6 g C
    out = field * 1e6 / area_m2

    return out, "g C m$^{-2}$ yr$^{-1}$"


def maybe_keep_mg_units(field, varname):
    if not is_carbon_like(varname):
        return field, "native units"

    return field, "Mg C yr$^{-1}$ per cell"


def mask_zero(field, eps=0.0):
    arr = np.asarray(field, dtype=float).copy()
    arr[np.abs(arr) <= eps] = np.nan
    return arr

def print_percentiles(field, name="field"):
    data = np.asarray(field, dtype=float)
    data = data[np.isfinite(data)]

    n_total = data.size

    print(f"\n{name} distribution:")
    print("total valid grids:", n_total)
    print("min:", np.nanmin(data))
    print("max:", np.nanmax(data))
    print()

    percentiles = [0.05, 0.1, 0.5, 1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99, 99.5, 99.9, 99.95]

    values = []
    for p in percentiles:
        v = np.nanpercentile(data, p)
        values.append(v)

    # 累计统计
    #print("Cumulative distribution:")
    #for p, v in zip(percentiles, values):
    #    count = np.sum(data <= v)
    #    print(f"p{p:>5}: {v:>12.6g} | grids ≤ value: {count:>8} ({count/n_total:.2%})")

    #print("\nBin counts between percentiles:")
    for i in range(len(percentiles) - 1):
        v0 = values[i]
        v1 = values[i + 1]

        count = np.sum((data > v0) & (data <= v1))
        print(f"{percentiles[i]:>5}-{percentiles[i+1]:>5}% : {count:>8} ({count/n_total:.2%})")

# ============================================================
# 4. GCB 风格固定色标
# ============================================================

def get_gcb_colormap(levels=None):
    """
    GCB 风格固定分级色标：
    负值蓝色，正值红色，0 是分界线。
    """

    if levels is None:
        levels = np.array([
            -400, -200, -100, -50, -25, -12.5,
            0,
            12.5, 25, 50, 100, 200, 400
        ], dtype=float)
    else:
        levels = np.asarray(levels, dtype=float)

    colors = [
        "#053061",
        "#2166AC",
        "#4393C3",
        "#92C5DE",
        "#D1E5F0",
        "#EAF3F8",
        "#FCEBE2",
        "#FDD0B6",
        "#F4A582",
        "#D6604D",
        "#B2182B",
        "#67001F",
    ]

    cmap = mcolors.ListedColormap(colors)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    cmap.set_bad((1, 1, 1, 0))

    norm = mcolors.BoundaryNorm(levels, cmap.N)

    return cmap, norm, levels

def get_sink_emission_colormap(field, legend_percentile=None):
    data = np.asarray(field, dtype=float)
    data = data[np.isfinite(data)]

    if legend_percentile is None:
        vmin = float(np.nanmin(data))
        vmax = float(np.nanmax(data))
    else:
        pmin, pmax = legend_percentile
        vmin = float(np.nanpercentile(data, pmin))
        vmax = float(np.nanpercentile(data, pmax))

    if vmin >= 0:
        vmin = 0.0
    if vmax <= 0:
        vmax = 0.0

    if np.isclose(vmin, vmax):
        eps = abs(vmax) * 0.01 if vmax != 0 else 1.0
        vmin -= eps
        vmax += eps

    zero_pos = (0.0 - vmin) / (vmax - vmin)

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "sink_emission",
        [
            (0.0, "#006837"),
            (zero_pos, "#f7f7f7"),
            (1.0, "#b10026"),
        ],
        N=256
    )
    cmap.set_bad((1, 1, 1, 0))

    norm = Normalize(vmin=vmin, vmax=vmax)
    ticks = np.linspace(vmin, vmax, 7)

    return cmap, norm, ticks

def plot_field(lat, lon, field, title, unit, use_cartopy=False, legend_percentile=None):
    field = np.asarray(field, dtype=float)
    field = np.where(np.isfinite(field), field, np.nan)

    if use_cartopy:
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature

            cmap, norm, ticks = get_sink_emission_colormap(field, legend_percentile=legend_percentile)

            fig = plt.figure(figsize=(13, 7))
            ax = plt.axes(projection=ccrs.Robinson())
            ax.set_global()

            ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
            ax.add_feature(cfeature.LAND, facecolor="white", zorder=0)
            ax.coastlines(resolution="110m", linewidth=0.8, color="0.65", zorder=3)

            ax.gridlines(
                crs=ccrs.PlateCarree(),
                linewidth=0.4,
                color="0.85",
                linestyle="-",
                alpha=0.8
            )

            mesh = ax.pcolormesh(
                lon,
                lat,
                field,
                transform=ccrs.PlateCarree(),
                cmap=cmap,
                norm=norm,
                shading="auto",
                zorder=2
            )

            cb = plt.colorbar(
                mesh,
                ax=ax,
                orientation="horizontal",
                pad=0.08,
                fraction=0.055,
                shrink=0.62,
                extend="both",
                ticks=ticks
            )

            cb.ax.tick_params(labelsize=11)
            cb.ax.set_xticklabels([f"{x:.0f}" for x in ticks])

            cb.set_label(
                "Carbon Net Exchange (Mg C yr$^{-1}$ per 0.25° cell)\n"
                "[<0: Sink (Green), >0: Emission (Red)]",
                fontsize=12
            )

            ax.set_title(
                "Global Carbon Sink and Emission Mean (2020–2020)",
                fontsize=16,
                pad=14
            )

            plt.tight_layout()
            return fig, ax

        except ImportError:
            print("[WARN] cartopy not installed. Falling back to normal matplotlib plot.")

    # use_cartopy=False：保持原来的 GCB gC/m2/yr 配色不变
    cmap, norm, levels = get_gcb_colormap()

    extent = [
        float(lon.min()),
        float(lon.max()),
        float(lat.min()),
        float(lat.max())
    ]

    lat_desc = bool(lat[0] > lat[-1]) if len(lat) > 1 else False
    origin = "upper" if lat_desc else "lower"

    fig = plt.figure(figsize=(12, 5.8))
    ax = plt.gca()

    im = ax.imshow(
        field,
        extent=extent,
        origin=origin,
        cmap=cmap,
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
    cb.ax.tick_params(labelsize=9)
    cb.ax.set_xticklabels([f"{x:g}" for x in levels])

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)

    plt.tight_layout()

    return fig, ax


# ============================================================
# 5. 主程序
# ============================================================

def main():
    with Dataset(ncfile, "r") as ds:
        print("Available variables:")
        for name in list_vars(ds):
            v = ds.variables[name]
            units = getattr(v, "units", "")
            print(f"{name:32s} dims={v.dimensions} shape={v.shape} units={units}")

        lat_name, lon_name = find_lat_lon(ds)
        lat = np.asarray(ds.variables[lat_name][:], dtype=float).squeeze()
        lon = np.asarray(ds.variables[lon_name][:], dtype=float).squeeze()

        field, actual_start, actual_end, n_years = get_mean_2d_field(
            ds=ds,
            varname=varname,
            start_year=start_year,
            end_year=end_year,
            time_offset=time_offset,
        )

        if actual_start is not None:
            title = (
                f"{varname}: mean field, "
                f"{actual_start}-{actual_end} "
                f"(n={n_years})"
            )
        else:
            title = f"{varname}: 2D field"

        print(
            f"Mean period: {actual_start}-{actual_end}, "
            f"time steps={n_years}"
        )

        print()
        print(title)
        print("raw shape:", field.shape)
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

        if use_cartopy:
            # cartopy 图使用 Mg C yr-1 per 0.25° cell
            field, unit = maybe_keep_mg_units(field, varname)
        else:
            # 普通图保持原来的 per-area GCB 风格
            if per_area:
                field, unit = convert_to_area_flux(field, area_m2, varname)
                title += " [per area]"
            else:
                field, unit = maybe_keep_mg_units(field, varname)

        if mask_zero_flag:
            field = mask_zero(field, zero_eps)
            title += " [zero masked]"

        title += f" [{unit}]"

        print("plot min/max:", np.nanmin(field), np.nanmax(field))
        print("unit:", unit)
        
        print_percentiles(field, name=varname)

        fig, ax = plot_field(
            lat=lat,
            lon=lon,
            field=field,
            title=title,
            unit=unit,
            use_cartopy=use_cartopy,
            legend_percentile=legend_percentile
        )

        if save_fig:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
            print(f"Saved to: {save_path}")
        else:
            plt.show()


if __name__ == "__main__":
    main()