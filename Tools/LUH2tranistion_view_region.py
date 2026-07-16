#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NC_bbox_stats.py
从 NetCDF 中提取 bbox 范围内的统计量：sum / mean / area-weighted mean
支持：
- 变量 (time, lat, lon) 或 (lat, lon)
- 时间选择：--time-slice 或 --years（若 time 可解码为 datetime）
- 输出：单值打印 或 保存 CSV（按时间输出序列）

依赖：xarray, numpy, pandas
# 1) 中国范围：对每个 time 求 bbox 内总和（输出时间序列）
python NC_bbox_stats.py summary.nc --var atmo --bbox china --stat sum --csv atmo_china_sum.csv

# 2) loess：求 bbox 内平均（非加权）
python NC_bbox_stats.py summary.nc --var soil_total --bbox loess --stat mean

# 3) xinjiang：面积加权平均（cos(lat) 近似）
python NC_bbox_stats.py summary.nc --var biomass_total --bbox xinjiang --stat mean --area-weighted

# 4) 自定义 bbox + 仅取 time 索引 31~40（END 不包含）
python NC_bbox_stats.py summary.nc --var atmo --bbox 73 135 18 54 --time-slice 31 41 --stat sum

"""

import argparse
import numpy as np
import xarray as xr
import pandas as pd


BBOX_PRESETS = {
    "global": (-180.0, 180.0, -90.0, 90.0),
    "china": (73.0, 135.0, 18.0, 54.0),
    "loess": (100.0, 114.0, 33.0, 41.5),
    "xinjiang": (73.0, 97.5, 34.0, 50.5),
}


def detect_lat_lon_names(ds: xr.Dataset, lat_name=None, lon_name=None):
    if lat_name and lon_name:
        return lat_name, lon_name

    lat_cand = ["lat", "latitude", "LAT", "nav_lat", "y", "Y"]
    lon_cand = ["lon", "longitude", "LON", "nav_lon", "x", "X"]

    for n in lat_cand:
        if n in ds.coords or n in ds.variables:
            lat_name = lat_name or n
            break
    for n in lon_cand:
        if n in ds.coords or n in ds.variables:
            lon_name = lon_name or n
            break

    if lat_name is None or lon_name is None:
        raise ValueError("找不到 lat/lon 变量名。请用 --lat-name/--lon-name 显式指定。")
    return lat_name, lon_name


def parse_bbox(bbox_args):
    # bbox_args: list[str]
    if len(bbox_args) == 1 and bbox_args[0].lower() in BBOX_PRESETS:
        return BBOX_PRESETS[bbox_args[0].lower()]
    if len(bbox_args) == 4:
        lonmin, lonmax, latmin, latmax = map(float, bbox_args)
        return lonmin, lonmax, latmin, latmax
    raise ValueError("--bbox 需要：preset 名称（china/loess/xinjiang/global）或 4 个数 lonmin lonmax latmin latmax")


def subset_bbox(da: xr.DataArray, lat_name: str, lon_name: str, bbox):
    lonmin, lonmax, latmin, latmax = bbox

    # 处理 lat 可能升序/降序
    lat = da[lat_name]
    lon = da[lon_name]

    if lat.values[0] < lat.values[-1]:
        da = da.sel({lat_name: slice(latmin, latmax)})
    else:
        da = da.sel({lat_name: slice(latmax, latmin)})

    # lon 一般升序；若是 0-360 这里不做自动转换（需要你明确 bbox 体系）
    da = da.sel({lon_name: slice(lonmin, lonmax)})
    return da


def cell_area_weights(lat: xr.DataArray):
    """
    简单面积权重：w = cos(lat)（适用于规则经纬网的“面积加权平均”）
    注意：这是权重平均，不是“总量”；总量需要乘 cell area（m2）才严谨。
    """
    w = np.cos(np.deg2rad(lat))
    # 防止极点数值问题
    w = xr.where(np.isfinite(w), w, 0.0)
    return w


def stats_over_space(da2d: xr.DataArray, lat_name: str, lon_name: str, stat: str, area_weighted: bool):
    if stat not in ("sum", "mean"):
        raise ValueError("--stat 只能是 sum 或 mean")

    if area_weighted and stat == "mean":
        w = cell_area_weights(da2d[lat_name])
        # xarray weighted: 需要 weights 能广播到 (lat, lon)
        w2 = w.broadcast_like(da2d)
        return da2d.weighted(w2).mean(dim=(lat_name, lon_name), skipna=True)

    # 非加权
    if stat == "sum":
        return da2d.sum(dim=(lat_name, lon_name), skipna=True)
    else:
        return da2d.mean(dim=(lat_name, lon_name), skipna=True)


def main():
    ap = argparse.ArgumentParser(description="Compute bbox sum/mean from NetCDF variable.")
    ap.add_argument("nc_path", help="NetCDF 文件路径")
    ap.add_argument("--var", required=True, help="变量名")
    ap.add_argument("--bbox", nargs="+", required=True,
                    help="bbox: preset 名称(china/loess/xinjiang/global) 或 4 个数 lonmin lonmax latmin latmax")
    ap.add_argument("--lat-name", default=None, help="lat 变量名（可选）")
    ap.add_argument("--lon-name", default=None, help="lon 变量名（可选）")

    # 时间选择（两种方式二选一或都不选）
    ap.add_argument("--time-slice", nargs=2, type=int, metavar=("START", "END"),
                    help="按 time 索引切片（python 切片语义）：START END（END 不包含）")
    ap.add_argument("--years", nargs=2, type=int, metavar=("Y1", "Y2"),
                    help="按年份选择（需要 time 能解码成 datetime 或带 year 坐标）")

    ap.add_argument("--stat", choices=["sum", "mean"], default="sum", help="空间统计：sum 或 mean")
    ap.add_argument("--area-weighted", action="store_true",
                    help="仅对 mean 有效：使用 cos(lat) 做面积加权平均（规则经纬网近似）")
    ap.add_argument("--csv", default=None, help="输出 CSV（如果变量有 time 维，将输出时间序列）")

    args = ap.parse_args()

    ds = xr.open_dataset(args.nc_path)

    if args.var not in ds:
        raise ValueError(f"{args.var} 不在文件中。可用变量：{list(ds.data_vars)}")

    lat_name, lon_name = detect_lat_lon_names(ds, args.lat_name, args.lon_name)

    da = ds[args.var]

    # --- 时间选择 ---
    if "time" in da.dims:
        if args.years is not None:
            y1, y2 = args.years
            # 优先使用可解码 datetime
            if np.issubdtype(da["time"].dtype, np.datetime64):
                years = da["time"].dt.year
                da = da.sel(time=(years >= y1) & (years <= y2))
            elif "year" in da.coords:
                da = da.sel(year=slice(y1, y2))
            else:
                raise ValueError("time 不是 datetime 且没有 year 坐标，无法用 --years。请改用 --time-slice。")

        elif args.time_slice is not None:
            s, e = args.time_slice
            da = da.isel(time=slice(s, e))

    else:
        # 没有 time 维：忽略时间参数
        pass

    # --- bbox 裁剪 ---
    bbox = parse_bbox(args.bbox)
    da = subset_bbox(da, lat_name, lon_name, bbox)

    # --- 统计 ---
    if "time" in da.dims:
        # 对每个时间步做空间统计，输出序列
        out = []
        for i in range(da.sizes["time"]):
            da2d = da.isel(time=i)
            val = stats_over_space(da2d, lat_name, lon_name, args.stat, args.area_weighted)
            out.append(float(val.values))

        # time 标签
        time_vals = da["time"].values
        df = pd.DataFrame({"time": time_vals, f"{args.var}_{args.stat}": out})

        if args.csv:
            df.to_csv(args.csv, index=False)
            print(f"[OK] Saved CSV: {args.csv}")
        else:
            print(df.to_string(index=False))

    else:
        # 2D：输出单值
        val = stats_over_space(da, lat_name, lon_name, args.stat, args.area_weighted)
        v = float(val.values)
        if args.csv:
            pd.DataFrame({f"{args.var}_{args.stat}": [v]}).to_csv(args.csv, index=False)
            print(f"[OK] Saved CSV: {args.csv}")
        else:
            print(f"{args.var} {args.stat} over bbox {bbox} = {v}")

    ds.close()


if __name__ == "__main__":
    main()
