#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NCOutput_to_csv3_updated.py

将当前 bookkeeping summary NetCDF 输出汇总为 CSV。

主要更新：
1. 自动识别所有 3D(time,lat,lon) 和 2D(lat,lon) 数据变量，不再依赖固定变量清单。
2. 兼容旧诊断变量和新的 LUCE/BLUE/GCB 分类变量。
3. 对年度通量变量直接汇总；对库变量可额外计算逐年差分。
4. 自动生成核心守恒诊断：C_sys_total, C_bookkeep_total。
5. 如果存在旧变量，则额外生成一组兼容性分类列，便于和最新输出口径对照。

用法：
python NCOutput_to_csv3_updated.py summary_1deg.global.nc summary_1deg.global.csv
python NCOutput_to_csv3_updated.py summary_1deg.global.nc summary.csv --no-diff
"""

from __future__ import annotations

import argparse
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import netCDF4 as nc
import numpy as np
import pandas as pd


COORD_NAMES = {
    "time", "Time", "TIME", "t",
    "lat", "latitude", "LAT", "nav_lat", "y", "Y",
    "lon", "longitude", "LON", "nav_lon", "x", "X",
}

TIME_CANDIDATES = ["time", "Time", "TIME", "t"]
LAT_CANDIDATES = ["lat", "latitude", "LAT", "nav_lat", "y", "Y"]
LON_CANDIDATES = ["lon", "longitude", "LON", "nav_lon", "x", "X"]

# 状态库变量：通常是累计库量，逐年差分才对应年变化。
STOCK_VARS = {"biomass_total", "soil_total", "P1", "P10", "P100", "atmosphere"}

# 明确不按 1e9 转换的面积/掩膜类变量。
NON_CARBON_KEYWORDS = (
    "area", "frac", "fraction", "mask", "done", "lat", "lon", "land", "cell",
)

# 推荐优先显示的变量顺序；不存在的自动跳过，其余变量追加在后面。
PREFERRED_ORDER = [
    "biomass_total", "soil_total", "P1", "P10", "P100", "atmosphere",
    "Gross_Sources", "Gross_Sinks", "Net_Emissions",
    "Flux_FD", "Flux_NFC", "Flux_FR", "Flux_NFR", "Flux_CAL", "Flux_WHp",
    "Flux_Clearing", "Flux_Abandonment", "Flux_Harvest_Net",
    "Flux_Harvest_SoilSlash", "Flux_Harvest_Regrowth", "Flux_Other",
    "area_clearing", "area_abandonment", "area_other", "area_harvest", "area_deforestation",
    "harvest_requested_biomass", "harvest_met_biomass", "harvest_unmet_biomass",
    "unmet_harvest",
]


def _find_first(ds: nc.Dataset, names: Sequence[str]) -> Optional[str]:
    for name in names:
        if name in ds.variables:
            return name
    return None


def _find_by_standard_name(ds: nc.Dataset, standard_name: str) -> Optional[str]:
    for name, var in ds.variables.items():
        if getattr(var, "standard_name", "").lower() == standard_name.lower():
            return name
    return None


def find_time_name(ds: nc.Dataset) -> Optional[str]:
    return _find_first(ds, TIME_CANDIDATES) or _find_by_standard_name(ds, "time")


def find_lat_lon_names(ds: nc.Dataset) -> Tuple[Optional[str], Optional[str]]:
    lat = _find_first(ds, LAT_CANDIDATES) or _find_by_standard_name(ds, "latitude")
    lon = _find_first(ds, LON_CANDIDATES) or _find_by_standard_name(ds, "longitude")
    return lat, lon


def get_years(ds: nc.Dataset) -> np.ndarray:
    time_name = find_time_name(ds)
    if time_name is None:
        # 如果没有 time，就按单个时间步处理。
        return np.array([0], dtype=int)
    return np.asarray(ds.variables[time_name][:]).astype(int)


def is_data_var(ds: nc.Dataset, name: str) -> bool:
    if name in COORD_NAMES:
        return False

    var = ds.variables[name]
    dims = tuple(var.dimensions)
    if len(dims) == 0:
        return False

    time_name = find_time_name(ds)
    lat_name, lon_name = find_lat_lon_names(ds)

    if len(dims) == 3 and time_name in dims and lat_name in dims and lon_name in dims:
        return True
    if len(dims) == 2 and lat_name in dims and lon_name in dims:
        return True
    return False


def sorted_data_vars(ds: nc.Dataset) -> List[str]:
    existing = [name for name in ds.variables if is_data_var(ds, name)]
    preferred = [name for name in PREFERRED_ORDER if name in existing]
    rest = sorted([name for name in existing if name not in preferred])
    return preferred + rest


def read_as_time_lat_lon(ds: nc.Dataset, name: str) -> np.ndarray:
    """返回 shape=(time, lat, lon) 的数组；2D 变量会扩展为单时间步。"""
    var = ds.variables[name]
    arr = np.asarray(var[:], dtype=float)
    dims = list(var.dimensions)

    time_name = find_time_name(ds)
    lat_name, lon_name = find_lat_lon_names(ds)

    if arr.ndim == 2:
        # 2D 输出文件，如 atmosphere_1960_2020.nc
        return arr[np.newaxis, :, :]

    if arr.ndim != 3:
        raise ValueError(f"Unsupported variable dims for {name}: {dims}")

    # 允许维度顺序不是严格 time,lat,lon，自动转置。
    target_dims = [time_name, lat_name, lon_name]
    axes = [dims.index(d) for d in target_dims]
    return np.transpose(arr, axes=axes)


def spatial_sum(arr: np.ndarray) -> np.ndarray:
    return np.nansum(arr, axis=(1, 2))


def is_non_carbon_var(name: str) -> bool:
    lname = name.lower()
    return any(k in lname for k in NON_CARBON_KEYWORDS)


def add_if_available(row: Dict[str, float], out_name: str, source_names: Iterable[str]) -> None:
    values = []
    for s in source_names:
        col = f"{s}_total"
        if col in row:
            values.append(row[col])
    if values:
        row[out_name] = float(np.nansum(values))


def add_compatibility_categories(row: Dict[str, float]) -> None:
    """
    如果新代码已经直接输出这些变量，本函数不会覆盖。
    如果还在使用旧诊断变量，则给出可对照的近似分类列。
    """
    # BLUE event-level compatibility
    if "Flux_Clearing_total" not in row:
        add_if_available(row, "Flux_Clearing_total", ["emit_clearing"])
    if "Flux_Abandonment_total" not in row:
        add_if_available(row, "Flux_Abandonment_total", ["emit_abandonment"])
    if "Flux_Other_total" not in row:
        add_if_available(row, "Flux_Other_total", ["emit_other"])
    if "Flux_WHp_total" not in row:
        add_if_available(row, "Flux_WHp_total", ["emit_products", "emit_products_clearing", "emit_products_harvest"])
    if "Flux_Harvest_Net_total" not in row:
        add_if_available(row, "Flux_Harvest_Net_total", ["emit_harvest", "emit_products_harvest"])
    if "Flux_Harvest_SoilSlash_total" not in row:
        add_if_available(row, "Flux_Harvest_SoilSlash_total", ["harvest_to_soil", "harvest_R", "harvest_delta_SS_h"])
    if "Flux_Harvest_Regrowth_total" not in row:
        # harvest_delta_B_h 往往可为负，作为采伐后生物量恢复/亏损的诊断项。
        add_if_available(row, "Flux_Harvest_Regrowth_total", ["harvest_delta_B_h"])

    # deforestation compatibility
    if "Flux_FD_total" not in row:
        add_if_available(row, "Flux_FD_total", [
            "deforestation_to_A", "deforestation_to_products", "deforestation_to_soil",
            "emit_products_clearing",
        ])

    # Macro flux compatibility：优先用 event/flux 类年度通量，不用 biomass/soil 等状态库。
    if "Gross_Sources_total" not in row or "Gross_Sinks_total" not in row or "Net_Emissions_total" not in row:
        candidates = []
        for name in [
            "Flux_FD", "Flux_NFC", "Flux_FR", "Flux_NFR", "Flux_CAL", "Flux_WHp",
            "Flux_Clearing", "Flux_Abandonment", "Flux_Harvest_Net",
            "Flux_Harvest_SoilSlash", "Flux_Harvest_Regrowth", "Flux_Other",
            "emit_clearing", "emit_abandonment", "emit_other", "emit_harvest", "emit_products",
        ]:
            col = f"{name}_total"
            if col in row:
                candidates.append(row[col])
        if candidates:
            arr = np.asarray(candidates, dtype=float)
            row.setdefault("Gross_Sources_total", float(np.nansum(arr[arr > 0])))
            row.setdefault("Gross_Sinks_total", float(np.nansum(arr[arr < 0])))
            row.setdefault("Net_Emissions_total", float(np.nansum(arr)))


def summarize_summary_nc(
    nc_path: str,
    out_csv: str,
    add_diff: bool = True,
    add_gtc_columns: bool = True,
) -> pd.DataFrame:
    with nc.Dataset(nc_path) as ds:
        years = get_years(ds)
        vars_to_read = sorted_data_vars(ds)
        if not vars_to_read:
            raise ValueError("No 2D/3D data variables found in this NetCDF file.")

        print("检测到的数据变量：")
        for name in vars_to_read:
            v = ds.variables[name]
            units = getattr(v, "units", "")
            print(f"  {name:32s} dims={v.dimensions} shape={v.shape} units={units}")

        data = {name: read_as_time_lat_lon(ds, name) for name in vars_to_read}
        nt = max(arr.shape[0] for arr in data.values())
        if len(years) != nt:
            if len(years) == 1 and nt > 1:
                years = np.arange(nt, dtype=int)
            elif nt == 1:
                years = years[:1]
            else:
                raise ValueError(f"time length mismatch: len(time)={len(years)}, data nt={nt}")

        totals_by_var = {name: spatial_sum(arr) for name, arr in data.items()}

        rows: List[Dict[str, float]] = []
        for t in range(nt):
            row: Dict[str, float] = {"Year": int(years[t])}
            for name in vars_to_read:
                vals = totals_by_var[name]
                idx = min(t, len(vals) - 1)
                row[f"{name}_total"] = float(vals[idx])

            # 核心库量与守恒检查。
            b = row.get("biomass_total_total", 0.0)
            s = row.get("soil_total_total", 0.0)
            p1 = row.get("P1_total", 0.0)
            p10 = row.get("P10_total", 0.0)
            p100 = row.get("P100_total", 0.0)
            a = row.get("atmosphere_total", 0.0)
            row["C_sys_total"] = b + s + p1 + p10 + p100
            row["C_bookkeep_total"] = row["C_sys_total"] + a

            add_compatibility_categories(row)
            rows.append(row)

        df = pd.DataFrame(rows)

    # 差分列：主要用于库变量和守恒检查；年度通量变量本身不要用 diff 解读。
    if add_diff:
        for col in list(df.columns):
            if col == "Year":
                continue
            df[f"{col}_diff"] = df[col].diff()

        if "atmosphere_total" in df.columns:
            df["ELUC_from_atmosphere_diff"] = df["atmosphere_total"].diff()
        if "C_sys_total" in df.columns:
            df["C_sys_change"] = df["C_sys_total"].diff()
        if "C_bookkeep_total" in df.columns:
            df["C_bookkeep_change"] = df["C_bookkeep_total"].diff()

    # 额外添加 Gt C 或 Gt C/yr 版本，面积变量不转换。
    if add_gtc_columns:
        for col in list(df.columns):
            if col == "Year" or not np.issubdtype(df[col].dtype, np.number):
                continue
            base = col
            for suffix in ["_diff", "_change"]:
                if base.endswith(suffix):
                    base = base[: -len(suffix)]
            var_name = base.replace("_total", "")
            if not is_non_carbon_var(var_name):
                df[f"{col}_GtC"] = df[col] / 1e9

    df.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"\n✅ 汇总完成：{out_csv}")
    print(f"行数={len(df)}, 列数={len(df.columns)}")
    print(df.head())
    return df


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize bookkeeping summary NetCDF to CSV.")
    parser.add_argument("ncfile", help="Input summary NetCDF file")
    parser.add_argument("out_csv", help="Output CSV path")
    parser.add_argument("--no-diff", action="store_true", help="Do not add year-to-year difference columns")
    parser.add_argument("--no-gtc", action="store_true", help="Do not add /1e9 GtC columns")
    args = parser.parse_args()

    summarize_summary_nc(
        args.ncfile,
        args.out_csv,
        add_diff=(not args.no_diff),
        add_gtc_columns=(not args.no_gtc),
    )
    return 0


if __name__ == "__main__":
    summarize_summary_nc(
        "summary_1deg.global_1950_4.28.nc",
        "summary_1deg.global_1950_4.28.csv",
        add_diff=True,
    )
