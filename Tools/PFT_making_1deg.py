#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_dominant_pft_to_match_LUH2.py

功能：
- 读取 PFT 源文件（PFTfrac 或 maxvegetfrac 的份额栈，或已有 dominant_pft）
- 若为 0.25°/0.5° 规则格且与 1°整齐对齐，则按面积权重聚合到 1°
- 读取 LUH2 的 states_1deg.nc 作为对齐模板，精确复制其 lat/lon 变量（1D或2D）
- 输出与 LUH2 完全一致的 1° dominant_pft(lat,lon)（int16），便于后续索引一致

依赖：
  xarray, numpy, netCDF4 (或 h5netcdf)
"""

import os
import numpy as np
import xarray as xr


# ========= 用户参数 =========
LUH2_STATES_1DEG = r"./states_1deg.nc"          # 对齐模板
PFT_SRC          = r"./ESA_CCI_LC_PFT_halfdeg_2000.nc"           # 源：含 PFTfrac / maxvegetfrac(veget,lat,lon) / dominant_pft
OUT_NC           = r"./dominant_pft_1deg_aligned.nc"

# 源数据中候选变量名（多通道份额）
CAND_MULTI = [("PFTfrac", "PFT"), ("maxvegetfrac", "veget")]  # (var_name, pft_dim_name)
# 是否将 ORCHIDEE 1..15 映射到你的内部PFT编码（如果源 already 是你的编码，可设为 False）
APPLY_ORCHIDEE_MAPPING = True
# ORCHIDEE(1..15) → 你的内部 PFT 索引的映射（按你之前的版本）
ORCHIDEE_TO_INTERNAL = np.array([0, 14, 1, 2, 4, 3, 5, 6, 8, 7, 10, 10, 10, 10, 9, 10], dtype=np.int32)
# ==========================


def to_1d_coord(da, kind):
    """把 lat/lon 坐标（可能是 1D 或 2D）收敛为 1D 维度坐标。"""
    arr = np.array(da.values, dtype=float)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2:
        if kind == "lat":
            # 每行经向是否全相等？是则取第一列
            if np.allclose(arr, arr[:, :1], equal_nan=True):
                return arr[:, 0]
            return np.sort(np.unique(np.round(arr, 10)))
        else:
            # lon：每列纬向相等？取第一行
            if np.allclose(arr, arr[:1, :], equal_nan=True):
                return arr[0, :]
            return np.sort(np.unique(np.round(arr, 10)))
    raise ValueError(f"{kind} coord ndim={arr.ndim} 不支持")


def cell_area_m2(lat_2d, dlat_deg, dlon_deg, R=6371000.0):
    """规则格的每格面积（m²）。"""
    dlat = np.deg2rad(dlat_deg); dlon = np.deg2rad(dlon_deg)
    return (R**2) * dlat * dlon * np.cos(np.deg2rad(lat_2d))


def detect_regular_res(coords):
    """从 1D 坐标推断规则步长（度）。"""
    u = np.unique(np.sort(coords))
    if u.size < 2:
        return None
    dif = np.diff(u)
    dif = dif[dif > 0]
    if dif.size == 0:
        return None
    d = float(np.median(dif))
    return d


def align_lon(lon, target_is_0360):
    """将经度坐标调整到与目标一致的体系（0..360 或 -180..180）。"""
    lon = np.asarray(lon, dtype=float).copy()
    if target_is_0360:
        lon = lon % 360.0
    else:
        lon = ((lon + 180.0) % 360.0) - 180.0
    order = np.argsort(lon)
    return lon[order], order


def main():
    # 读 LUH2 1° 模板（不解码时间，避免 cftime 问题）
    ds_tmpl = xr.open_dataset(LUH2_STATES_1DEG, decode_times=False)
    # 获取模板的 lat/lon 名与数据（优先 lat/lon，其次 latitude/longitude）
    lat_name = "lat" if "lat" in ds_tmpl.coords or "lat" in ds_tmpl.variables else "latitude"
    lon_name = "lon" if "lon" in ds_tmpl.coords or "lon" in ds_tmpl.variables else "longitude"
    lat_tmpl = ds_tmpl[lat_name]
    lon_tmpl = ds_tmpl[lon_name]
    # 目标网格是 1D 还是 2D？
    lat1d = to_1d_coord(lat_tmpl, "lat")
    lon1d = to_1d_coord(lon_tmpl, "lon")
    target_is_0360 = (np.nanmin(lon1d) >= 0) and (np.nanmax(lon1d) > 180)
    # 推断目标分辨率
    dlat_t = detect_regular_res(lat1d); dlon_t = detect_regular_res(lon1d)
    if not (abs(dlat_t - 1.0) < 1e-6 and abs(dlon_t - 1.0) < 1e-6):
        print(f"[warn] 模板分辨率检测为 {dlat_t}°×{dlon_t}°，不是 1°，仍按模板对齐写出。")
    # 构造目标 2D 纬度场用于面积（1D→2D；2D模板我们只复用坐标不需要其格点纬度场）
    LAT2D = np.repeat(lat1d[:, None], lon1d.size, axis=1)
    area_1deg_m2 = cell_area_m2(LAT2D, dlat_t, dlon_t)

    # 读 PFT 源
    ds_src = xr.open_dataset(PFT_SRC, decode_times=False)

    # 情况 A：已有 dominant_pft(lat,lon) —— 直接对齐输出
    if "dominant_pft" in ds_src.variables:
        dom = ds_src["dominant_pft"].load()
        # squeeze 成 2D
        while dom.ndim > 2:
            # 尝试去掉最前面的非lat/lon维
            dims = list(dom.dims)
            drop_axis = 0
            dom = dom.isel({dims[drop_axis]: 0})
        # 调整经度到与模板一致
        src_lon_name = "lon" if "lon" in dom.coords or "lon" in ds_src.variables else "longitude"
        src_lat_name = "lat" if "lat" in dom.coords or "lat" in ds_src.variables else "latitude"
        lon_src = ds_src[src_lon_name].values
        lon_adj, order = align_lon(lon_src, target_is_0360)
        # 重排经度
        dom = dom.transpose(src_lat_name, src_lon_name)
        dom = dom.isel({src_lon_name: xr.DataArray(order, dims=(src_lon_name,))})
        # 重新标注目标坐标名和值（严格对齐 LUH2 模板）
        dom = dom.assign_coords({src_lat_name: lat1d, src_lon_name: lon1d})
        # ORCHIDEE→内部映射（如需）
        arr = dom.values
        if APPLY_ORCHIDEE_MAPPING and np.nanmin(arr) >= 1 and np.nanmax(arr) <= 15:
            arr_int = np.array(np.rint(arr), dtype=np.int32)
            mapped = ORCHIDEE_TO_INTERNAL[arr_int]
        else:
            mapped = np.array(np.rint(arr), dtype=np.int32)

    else:
        # 情况 B：多通道份额，先确定变量名与 PFT 维
        varname = None; pft_dim = None
        for vname, dname in CAND_MULTI:
            if vname in ds_src.variables:
                dims = list(ds_src[vname].dims)
                if dname in dims:
                    varname, pft_dim = vname, dname
                    break
        if varname is None:
            raise KeyError("未找到可用的 PFT 份额变量（期望 PFTfrac(PFT,lat,lon) 或 maxvegetfrac(veget,lat,lon)，或 dominant_pft）。")

        da = ds_src[varname].load()  # (PFT,lat,lon)
        # 统一经纬度名
        src_lon_name = "lon" if "lon" in da.coords or "lon" in ds_src.variables else "longitude"
        src_lat_name = "lat" if "lat" in da.coords or "lat" in ds_src.variables else "latitude"
        # 经度对齐（0..360 或 -180..180）
        lon_src = ds_src[src_lon_name].values
        lon_adj, order = align_lon(lon_src, target_is_0360)
        da = da.transpose(pft_dim, src_lat_name, src_lon_name)
        da = da.isel({src_lon_name: xr.DataArray(order, dims=(src_lon_name,))})

        # 检测源分辨率，判断是否 0.25°/0.5° 与 1°整齐对齐
        lat_src_1d = to_1d_coord(da[src_lat_name], "lat")
        lon_src_1d = to_1d_coord(da[src_lon_name], "lon")
        dlat_s = detect_regular_res(lat_src_1d)
        dlon_s = detect_regular_res(lon_src_1d)

        if dlat_s is None or dlon_s is None:
            raise ValueError("源坐标不是规则格，当前脚本不支持（可改为更通用的分箱版）。")

        fac_lat = int(round(dlat_t / dlat_s)) if dlat_t and dlat_s else None
        fac_lon = int(round(dlon_t / dlon_s)) if dlon_t and dlon_s else None

        if not (fac_lat in (2,4) and fac_lon in (2,4) and abs(dlat_t - fac_lat*dlat_s) < 1e-6 and abs(dlon_t - fac_lon*dlon_s) < 1e-6):
            raise ValueError(f"暂仅支持 0.5°(×2) 或 0.25°(×4) → 1° 的整齐聚合。检测到源 {dlat_s}×{dlon_s}，目标 {dlat_t}×{dlon_t}，无法整齐聚合。")

        # 面积权重聚合（对每个 PFT 通道求 sum(fraction * area)）
        LAT2D_src, LON2D_src = xr.broadcast(da[src_lat_name], da[src_lon_name])
        area_src_m2 = cell_area_m2(LAT2D_src, dlat_s, dlon_s)
        weighted = da * area_src_m2  # (PFT,lat,lon)
        summed = weighted.coarsen({src_lat_name: fac_lat, src_lon_name: fac_lon}, boundary="trim").sum()
        # 在 PFT 维上 argmax
        # 先把 NaN → -inf，避免 argmax 受 NaN 干扰
        summed_filled = summed.where(np.isfinite(summed), other=-np.inf)
        dom_idx = summed_filled.argmax(dim=pft_dim)  # 0-based
        orch = (dom_idx + 1).astype(np.int32)        # 1..N（按 ORCHIDEE 习惯）

        if APPLY_ORCHIDEE_MAPPING:
            mapped = ORCHIDEE_TO_INTERNAL[orch.values]
        else:
            mapped = orch.values

        # 此时的 (lat,lon) 是聚合后的中心；用 LUH2 模板坐标替换
        # 确保纬度方向与模板一致（升/降序）。
        mapped = np.array(mapped, dtype=np.int32)

    # ===== 写出，与 LUH2 模板完全一致的坐标与维度 =====
    # 若模板 lat/lon 是 2D 变量，则保持 1D 坐标 + 2D 数据也是 OK（CF 兼容）；更严格可复制 2D 变量。
    # 这里采用：数据 2D(dims=lat_name,lon_name)，坐标使用模板的 lat/lon DataArray（保留属性）
    da_out = xr.DataArray(
        data=mapped.astype("int16"),
        dims=(lat_name, lon_name),
        coords={lat_name: lat_tmpl, lon_name: lon_tmpl},
        name="dominant_pft",
        attrs={
            "long_name": "Dominant PFT aligned to LUH2 1-deg grid",
            "units": "-",
            "description": "Area-weighted aggregation to 1° (if needed), then argmax over PFT; coordinates copied from LUH2 states_1deg.nc",
            #"_FillValue": np.int16(-32767),
        },
    )

    # 如果你希望将无效值写成填充值，可在写盘前先确保没有 NaN（int16 不支持 NaN）
    fill = np.int16(-32767)
    if np.issubdtype(da_out.dtype, np.integer):
        # 如果 mapped 里可能有 <0 或非常大的无效编码，也可以统一改成 fill
        pass
    else:
        # 一般不会进来，因为上面已经 astype("int16")
        da_out = da_out.fillna(fill).astype("int16")
    
    ds_out = xr.Dataset({"dominant_pft": da_out})
    
    # 只在 encoding 里指定填充值
    encoding = {"dominant_pft": {"_FillValue": fill}}
    
    ds_out.to_netcdf(OUT_NC, encoding=encoding)
    print("✅ Saved:", OUT_NC)
    print("shape:", ds_out["dominant_pft"].shape)
    print("lat head/tail:", np.asarray(lat1d)[:3], np.asarray(lat1d)[-3:])
    print("lon head/tail:", np.asarray(lon1d)[:3], np.asarray(lon1d)[-3:])


if __name__ == "__main__":
    main()
