#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
from netCDF4 import Dataset
from typing import List, Optional, Sequence, Tuple
import geopandas as gpd

# ============================================================
# 1. 用户设置
# ============================================================

ncfile = r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\summary_025deg.global_time_ge1000_7.9.nc"

# 可选示例：
# varname = "Net_Emissions"
# varname = "Gross_Sources" Gross_Sinks
varname = "Gross_Sources"

# 多年平均时间范围（公历年份，包含起止年份）
start_year = 1961
end_year = 2020

# 如果 NetCDF 的 time=0 对应公元850年，则设为850；
# 如果 time 本身就是公历年份，则设为0。
time_offset = 850

# 空间范围
bbox = None
# bbox = "china"
# bbox = "loess"
# bbox = (73, 135, 18, 54)

# 普通碳变量的单位处理：
# True  : Mg C yr-1 per cell -> g C m-2 yr-1
# False : 保持 Mg C yr-1 per cell
#
# 注意：Gross_Sources 和 Gross_Sinks 会自动忽略此设置，
# 固定转换为 10^10 g C yr-1 per cell，并从0开始显示正值。
per_area = True

mask_zero_flag = False
zero_eps = 1.0

# 没安装 cartopy 就设为 False
use_cartopy = False

# None 使用完整数据范围；如 (1, 99) 则用第1和99百分位控制色标。
# Gross Sources/Sinks 只使用上限百分位。
legend_percentile = 10
# legend_percentile = (1, 99)
# legend_percentile = (0.5, 99.5)

save_fig = False
save_path = r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\Gross_Sinks_mean.png"
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
    lat_name = _find_first(ds, LAT_CANDIDATES) or _find_by_attr(
        ds, "standard_name", "latitude"
    )
    lon_name = _find_first(ds, LON_CANDIDATES) or _find_by_attr(
        ds, "standard_name", "longitude"
    )

    if lat_name is None or lon_name is None:
        raise ValueError("Could not find lat/lon variables.")

    return lat_name, lon_name


def find_time(ds: Dataset) -> Optional[str]:
    return _find_first(ds, TIME_CANDIDATES) or _find_by_attr(
        ds, "standard_name", "time"
    )


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

    if (
        len(dims) == 3
        and time_name is not None
        and time_name in dims
        and lat_name in dims
        and lon_name in dims
    ):
        return True

    return False


def list_vars(ds: Dataset) -> List[str]:
    return sorted([name for name in ds.variables if is_data_var(ds, name)])


def get_time_values(ds: Dataset) -> Optional[np.ndarray]:
    time_name = find_time(ds)
    if time_name is None or time_name not in ds.variables:
        return None
    return np.asarray(ds.variables[time_name][:]).astype(int)


def get_mean_2d_field(
    ds: Dataset,
    varname: str,
    start_year: Optional[int] = None,
    end_year: Optional[int] = None,
    time_offset: int = 0,
) -> Tuple[np.ndarray, Optional[int], Optional[int], int]:
    """
    读取指定年份范围并计算多年平均。

    支持：
    1. 二维变量 (lat, lon)
    2. 三维变量 (time, lat, lon)，维度顺序可变

    返回：
    mean_field, actual_start, actual_end, n_steps
    """
    if varname not in ds.variables:
        raise ValueError(f"Variable '{varname}' not found.")

    v = ds.variables[varname]
    dims = list(v.dimensions)

    lat_name, lon_name = find_lat_lon(ds)
    time_name = find_time(ds)

    if len(dims) == 2:
        arr = np.ma.filled(v[:], np.nan).astype(float)
        axes = [dims.index(lat_name), dims.index(lon_name)]
        arr = np.transpose(arr, axes=axes)
        return arr, None, None, 1

    if len(dims) != 3 or time_name is None or time_name not in dims:
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

    time_indices = np.where(
        (calendar_years >= start_year) & (calendar_years <= end_year)
    )[0]

    if time_indices.size == 0:
        raise ValueError(
            f"No time steps found for {start_year}-{end_year}. "
            f"Available calendar years: "
            f"{calendar_years.min()}-{calendar_years.max()}"
        )

    time_axis = dims.index(time_name)
    lat_axis = dims.index(lat_name)
    lon_axis = dims.index(lon_name)

    # 当前年份通常连续，只读取所需时间段，避免加载整个文件。
    first_idx = int(time_indices[0])
    last_idx = int(time_indices[-1])

    slices = [slice(None)] * 3
    slices[time_axis] = slice(first_idx, last_idx + 1)

    arr3 = np.ma.filled(v[tuple(slices)], np.nan).astype(float)
    arr3 = np.transpose(arr3, axes=[time_axis, lat_axis, lon_axis])
    arr3 = np.where(np.isfinite(arr3), arr3, np.nan)

    mean_field = np.nanmean(arr3, axis=0)

    actual_start = int(calendar_years[first_idx])
    actual_end = int(calendar_years[last_idx])

    return mean_field, actual_start, actual_end, int(time_indices.size)


def parse_bbox_arg(bbox_value):
    if bbox_value is None:
        return None

    if isinstance(bbox_value, str):
        key = bbox_value.lower()
        if key not in BBOX_PRESETS:
            raise ValueError(f"Unknown bbox preset: {bbox_value}")
        return BBOX_PRESETS[key]

    if len(bbox_value) == 4:
        return tuple(map(float, bbox_value))

    raise ValueError("bbox must be None, preset name, or 4 numbers.")


def subset_by_bbox(lat, lon, field, bbox_value):
    lonmin, lonmax, latmin, latmax = bbox_value

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

    if lon.size == 1:
        lon_edges = np.array([lon[0] - 0.5, lon[0] + 0.5])
    else:
        lon_edges = _coord_edges_from_centers(
            lon,
            lon.min() - abs(lon[1] - lon[0]),
            lon.max() + abs(lon[-1] - lon[-2]),
        )

    dlon = np.abs(np.diff(np.deg2rad(lon_edges)))
    sin_lat = np.sin(np.deg2rad(lat_edges))
    d_sin_lat = np.abs(np.diff(sin_lat))

    return radius**2 * d_sin_lat[:, None] * dlon[None, :]


def is_carbon_like(varname: str) -> bool:
    lname = varname.lower()
    return not any(k in lname for k in NON_CARBON_KEYWORDS)


def is_gross_source_sink(varname: str) -> bool:
    return varname.lower() in {"gross_sources", "gross_sinks"}


def convert_to_area_flux(field, area_m2, varname):
    if not is_carbon_like(varname):
        raise ValueError("per_area should only be used for carbon-like variables.")

    # 模型输出单位：Mg C yr-1 per cell
    # 1 Mg C = 1e6 g C
    out = field * 1e6 / area_m2
    return out, "g C m$^{-2}$ yr$^{-1}$"


def convert_gross_flux_units(field, varname):
    """
    Gross Sources/Sinks:
    Mg C yr-1 per cell -> 10^10 g C yr-1 per cell

    1 Mg = 10^6 g，因此除以10^4。
    Gross Sinks和Gross Sources均显示为正的强度值。
    """
    if not is_gross_source_sink(varname):
        raise ValueError("This conversion is only for Gross_Sources/Gross_Sinks.")

    arr = np.abs(np.asarray(field, dtype=float))
    arr = arr / 1e4

    unit = r"Gross Sources ($\times10^{10}$ g C yr$^{-1}$)"
    return arr, unit


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

    if data.size == 0:
        print(f"\n{name} distribution: no finite values")
        return

    n_total = data.size

    print(f"\n{name} distribution:")
    print("total valid grids:", n_total)
    print("min:", np.nanmin(data))
    print("max:", np.nanmax(data))
    print()

    percentiles = [
        0.05, 0.1, 0.5, 1, 2, 5, 10, 25, 50,
        75, 90, 95, 98, 99, 99.5, 99.9, 99.95
    ]

    values = [np.nanpercentile(data, p) for p in percentiles]

    for i in range(len(percentiles) - 1):
        v0 = values[i]
        v1 = values[i + 1]
        count = np.sum((data > v0) & (data <= v1))
        print(
            f"{percentiles[i]:>5}-{percentiles[i+1]:>5}% : "
            f"{count:>8} ({count/n_total:.2%})"
        )


# ============================================================
# 4. 色标
# ============================================================

def get_gcb_colormap(levels=None):
    """
    普通净通量变量使用的GCB风格固定分级色标：
    负值蓝色，正值红色，0为分界线。
    """
    if levels is None:
        levels = np.array(
            [-400, -200, -100, -50, -25, -12.5,
             0,
             12.5, 25, 50, 100, 200, 400],
            dtype=float,
        )
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

    if data.size == 0:
        vmin, vmax = -1.0, 1.0
    elif legend_percentile is None:
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
        N=256,
    )
    cmap.set_bad((1, 1, 1, 0))

    norm = Normalize(vmin=vmin, vmax=vmax)
    ticks = np.linspace(vmin, vmax, 7)

    return cmap, norm, ticks


def get_positive_colormap(
    field,
    varname,
    legend_percentile=None,
):
    """
    Gross Sources / Gross Sinks 专用离散色标

    单位:
    10^10 g C yr-1 per cell

    根据变量名自动选择颜色:
    Gross_Sources -> 橙红色
    Gross_Sinks   -> 绿色
    """

    # 固定分级
    levels = np.array([
        0,
        0.5,
        1,
        1.5,
        2,
        2.5,
        3,
        6,
        10
    ], dtype=float)


    name = varname.lower()


    # ==========================
    # Gross Sources
    # 白 -> 黄 -> 橙 -> 红 -> 深红
    # ==========================
    if "source" in name:

        colors = [
            "#fff7ec",  # 0-0.5
            "#fee8c8",  # 0.5-1
            "#fdd49e",  # 1-1.5
            "#fdbb84",  # 1.5-2
            "#fc8d59",  # 2-2.5
            "#ef6548",  # 2.5-3
            "#d94801",  # 3-6
            "#7f2704",  # 6-10
        ]


    # ==========================
    # Gross Sinks
    # 白 -> 浅蓝 -> 青绿 -> 深绿
    # ==========================
    elif "sink" in name:

        colors = [
            "#f7fbff",  # 0-0.5
            "#deebf7",  # 0.5-1
            "#c7e9f1",  # 1-1.5
            "#9bd3ce",  # 1.5-2
            "#7bccc4",  # 2-2.5
            "#41b6a4",  # 2.5-3
            "#238b45",  # 3-6
            "#00441b",  # 6-10
        ]


    else:

        raise ValueError(
            f"Unknown gross variable: {varname}"
        )


    cmap = mcolors.ListedColormap(colors)

    # 小于0不用显示
    cmap.set_under(colors[0])

    # 大于10使用深色
    cmap.set_over(colors[-1])

    cmap.set_bad((1,1,1,0))


    norm = mcolors.BoundaryNorm(
        boundaries=levels,
        ncolors=cmap.N,
        clip=False
    )


    # colorbar刻度
    ticks = levels[1:]


    return cmap, norm, ticks


def plot_field(
    lat,
    lon,
    field,
    title,
    unit,
    varname,
    use_cartopy=False,
    legend_percentile=None,
):
    field = np.asarray(field, dtype=float)
    field = np.where(np.isfinite(field), field, np.nan)

    gross_mode = is_gross_source_sink(varname)

    if gross_mode:
        cmap, norm, ticks = get_positive_colormap(
            field,
            varname=varname,
            legend_percentile=legend_percentile,
        )
        extend_type = "max"
    else:
        cmap, norm, ticks = get_sink_emission_colormap(
            field,
            legend_percentile=legend_percentile,
        )
        extend_type = "both"

    if use_cartopy:
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature

            fig = plt.figure(figsize=(13, 7))
            ax = plt.axes(projection=ccrs.Robinson())

            if bbox is None or (
                isinstance(bbox, str) and bbox.lower() == "global"
            ):
                ax.set_global()

            ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
            ax.add_feature(cfeature.LAND, facecolor="white", zorder=0)
            ax.coastlines(
                resolution="110m",
                linewidth=0.8,
                color="0.65",
                zorder=3,
            )

            ax.gridlines(
                crs=ccrs.PlateCarree(),
                linewidth=0.4,
                color="0.85",
                linestyle="-",
                alpha=0.8,
            )

            mesh = ax.pcolormesh(
                lon,
                lat,
                field,
                transform=ccrs.PlateCarree(),
                cmap=cmap,
                norm=norm,
                shading="auto",
                zorder=2,
            )

            cb = plt.colorbar(
                mesh,
                ax=ax,
                orientation="horizontal",
                pad=0.08,
                fraction=0.055,
                shrink=0.68,
                extend=extend_type,
                ticks=ticks,
            )

            cb.ax.tick_params(labelsize=10)
            cb.ax.set_xticklabels(
                [f"{x:g}" for x in ticks]
            )
            cb.set_label(unit, fontsize=12)

            ax.set_title(title, fontsize=15, pad=14)
            plt.tight_layout()
            return fig, ax

        except ImportError:
            print(
                "[WARN] cartopy not installed. "
                "Falling back to normal matplotlib plot."
            )

    extent = [
        float(lon.min()),
        float(lon.max()),
        float(lat.min()),
        float(lat.max()),
    ]

    lat_desc = bool(lat[0] > lat[-1]) if len(lat) > 1 else False
    origin = "upper" if lat_desc else "lower"

    fig = plt.figure(figsize=(12, 5.8))
    ax = plt.gca()

    if gross_mode:
        im = ax.imshow(
            field,
            extent=extent,
            origin=origin,
            cmap=cmap,
            norm=norm,
            aspect="auto",
        )

        cb = plt.colorbar(
            im,
            ax=ax,
            orientation="horizontal",
            pad=0.12,
            fraction=0.055,
            shrink=0.95,
            ticks=ticks,
            extend="max" if is_gross_source_sink(varname) else "both",
        )
        cb.ax.set_xticklabels([f"{x:.3g}" for x in ticks])

    else:
        # 普通变量仍保留固定GCB分级色标
        cmap_gcb, norm_gcb, levels = get_gcb_colormap()

        im = ax.imshow(
            field,
            extent=extent,
            origin=origin,
            cmap=cmap_gcb,
            norm=norm_gcb,
            aspect="auto",
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
            extend="both",
        )
        cb.ax.set_xticklabels([f"{x:g}" for x in levels])

    cb.set_label(unit)
    cb.ax.tick_params(labelsize=9)

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
            print(
                f"{name:32s} "
                f"dims={v.dimensions} shape={v.shape} units={units}"
            )

        lat_name, lon_name = find_lat_lon(ds)
        lat_full = np.asarray(
            ds.variables[lat_name][:], dtype=float
        ).squeeze()
        lon_full = np.asarray(
            ds.variables[lon_name][:], dtype=float
        ).squeeze()

        field, actual_start, actual_end, n_years = get_mean_2d_field(
            ds=ds,
            varname=varname,
            start_year=start_year,
            end_year=end_year,
            time_offset=time_offset,
        )

        if actual_start is not None:
            period_text = f"{actual_start}-{actual_end}"
            title = f"{varname}: mean, {period_text} (n={n_years})"
        else:
            period_text = "2D field"
            title = f"{varname}: {period_text}"

        print(
            f"\nMean period: {actual_start}-{actual_end}, "
            f"time steps={n_years}"
        )
        print("raw shape:", field.shape)
        print("raw min/max:", np.nanmin(field), np.nanmax(field))

        area_m2 = grid_cell_area_m2(lat_full, lon_full)
        lat = lat_full.copy()
        lon = lon_full.copy()

        bbox_value = parse_bbox_arg(bbox)
        if bbox_value is not None:
            lat, lon, field = subset_by_bbox(
                lat_full, lon_full, field, bbox_value
            )
            _, _, area_m2 = subset_by_bbox(
                lat_full, lon_full, area_m2, bbox_value
            )
            title += f" [bbox={bbox}]"

        # Gross Sources/Sinks自动使用专用转换和正值色标
        if is_gross_source_sink(varname):
            field, unit = convert_gross_flux_units(field, varname)

            if varname.lower() == "gross_sinks":
                title = title.replace(
                    "Gross_Sinks",
                    "Gross Sinks magnitude",
                )
            else:
                title = title.replace(
                    "Gross_Sources",
                    "Gross Sources magnitude",
                )

        elif use_cartopy:
            field, unit = maybe_keep_mg_units(field, varname)

        else:
            if per_area:
                field, unit = convert_to_area_flux(
                    field,
                    area_m2,
                    varname,
                )
                title += " [per area]"
            else:
                field, unit = maybe_keep_mg_units(
                    field,
                    varname,
                )

        if mask_zero_flag:
            field = mask_zero(field, zero_eps)
            title += " [zero masked]"

        title += f" [{unit}]"

        finite = field[np.isfinite(field)]
        if finite.size == 0:
            raise ValueError(
                "No finite values remain after processing/masking."
            )

        print("plot min/max:", np.nanmin(field), np.nanmax(field))
        print("unit:", unit)

        print_percentiles(field, name=varname)

        fig, _ = plot_field(
            lat=lat,
            lon=lon,
            field=field,
            title=title,
            unit=unit,
            varname=varname,
            use_cartopy=use_cartopy,
            legend_percentile=legend_percentile,
        )

        if save_fig:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
            print(f"Saved to: {save_path}")
        else:
            plt.show()


if __name__ == "__main__":
    main()
