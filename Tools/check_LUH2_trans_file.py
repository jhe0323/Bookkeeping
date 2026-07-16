# diagnose_luh2_double_count.py
# -*- coding: utf-8 -*-

from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional
from pathlib import Path
import netCDF4 as nc
import numpy as np
import pandas as pd
import yaml


# =========================================================
# 1. 配置读取
# =========================================================

def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_luh2_to_lulc_mapping(cfg: dict) -> Dict[str, str]:
    """
    返回类似:
    primf -> v
    primn -> v
    secdf -> s
    ...
    """
    return dict(cfg["LUH2toLULC"])


def get_num_to_luh2_mapping(cfg: dict) -> Dict[int, str]:
    """
    返回类似:
    1 -> primf
    2 -> primn
    ...
    """
    return {int(k): v for k, v in cfg["numtoLUH2Type"].items()}


# =========================================================
# 2. 坐标/索引工具
# =========================================================

@dataclass
class GridMeta:
    lat: np.ndarray
    lon: np.ndarray
    lat_desc: bool
    lon_0_360: bool
    lat_name: str
    lon_name: str


def guess_lat_lon_names(ds: nc.Dataset) -> Tuple[str, str]:
    candidates_lat = ("lat", "latitude", "y")
    candidates_lon = ("lon", "longitude", "x")
    lat_name = next((n for n in candidates_lat if n in ds.variables), None)
    lon_name = next((n for n in candidates_lon if n in ds.variables), None)
    if lat_name is None or lon_name is None:
        raise KeyError(f"Cannot find lat/lon in dataset: {list(ds.variables.keys())[:50]}")
    return lat_name, lon_name


def load_grid_meta(ds: nc.Dataset) -> GridMeta:
    lat_name, lon_name = guess_lat_lon_names(ds)
    lat = np.asarray(ds.variables[lat_name][:], dtype=float).ravel()
    lon = np.asarray(ds.variables[lon_name][:], dtype=float).ravel()
    lat_desc = bool(lat.size >= 2 and np.all(np.diff(lat) < 0))
    lon_0_360 = bool(np.nanmin(lon) >= 0 and np.nanmax(lon) > 180)
    return GridMeta(
        lat=lat,
        lon=lon,
        lat_desc=lat_desc,
        lon_0_360=lon_0_360,
        lat_name=lat_name,
        lon_name=lon_name,
    )


def wrap_lon_180(lon: float) -> float:
    return ((lon + 180.0) % 360.0) - 180.0


def normalize_lon_for_dataset(lon: float, lon_0_360: bool) -> float:
    if lon_0_360:
        return (lon + 360.0) % 360.0
    return wrap_lon_180(lon)


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def ds_lat_idx(meta: GridMeta, i_internal: int) -> int:
    n = len(meta.lat)
    return (n - 1 - i_internal) if meta.lat_desc else i_internal


def build_internal_grid_from_states(state_meta: GridMeta) -> Tuple[np.ndarray, np.ndarray, float, float]:
    lat_clean = np.asarray(state_meta.lat, dtype=float).ravel()
    lon_clean = np.asarray(state_meta.lon, dtype=float).ravel()
    lat_clean = lat_clean[np.isfinite(lat_clean)]
    lon_clean = lon_clean[np.isfinite(lon_clean)]

    lat_u = np.unique(np.sort(lat_clean))
    lon_u = np.unique(np.sort(lon_clean))

    if lat_u.size >= 2:
        dlat = float(np.median(np.diff(lat_u)[np.diff(lat_u) > 0]))
    else:
        dlat = 0.0
    if lon_u.size >= 2:
        dlon = float(np.median(np.diff(lon_u)[np.diff(lon_u) > 0]))
    else:
        dlon = 0.0

    lat_res = abs(dlat)
    lon_res = abs(dlon)

    nlat = lat_u.size
    nlon = int(round(360.0 / lon_res))
    lat_internal = np.linspace(-90.0 + 0.5 * lat_res, 90.0 - 0.5 * lat_res, nlat)
    lon_internal = np.linspace(-180.0 + 0.5 * lon_res, 180.0 - 0.5 * lon_res, nlon)

    return lat_internal, lon_internal, lat_res, lon_res


def ds_lon_idx_from_internal_j(
    j_internal: int,
    lon_res: float,
    trans_meta: GridMeta,
) -> int:
    lon_val_internal = -180.0 + 0.5 * lon_res + j_internal * lon_res
    lon_for_ds = normalize_lon_for_dataset(lon_val_internal, trans_meta.lon_0_360)
    return nearest_index(np.asarray(trans_meta.lon, dtype=float), lon_for_ds)


def ds_lon_idx_from_lon_value(
    lon_value_180: float,
    trans_meta: GridMeta,
) -> int:
    lon_for_ds = normalize_lon_for_dataset(lon_value_180, trans_meta.lon_0_360)
    return nearest_index(np.asarray(trans_meta.lon, dtype=float), lon_for_ds)


# =========================================================
# 3. 面积网格
# =========================================================

def build_area_grid_ha(lat_internal: np.ndarray, lon_internal: np.ndarray, lat_res: float, lon_res: float) -> np.ndarray:
    """
    返回每个网格面积，单位 ha
    """
    R = 6371000.0
    dphi = np.deg2rad(lat_res)
    dlmb = np.deg2rad(lon_res)

    LAT = np.repeat(lat_internal[:, None], len(lon_internal), axis=1)
    area_m2 = (R ** 2) * dphi * dlmb * np.cos(np.deg2rad(LAT))
    return area_m2 / 1e4


# =========================================================
# 4. 状态/转移读取
# =========================================================

COARSE_KEYS = ["v", "s", "p", "c", "U"]


def safe_read_scalar(var, t: int, i: int, j: int) -> float:
    try:
        v = var[t, i, j]
        if hasattr(v, "filled"):
            v = v.filled(np.nan)
        out = float(v)
    except Exception:
        out = np.nan

    if not np.isfinite(out):
        return 0.0
    return out


def get_state_frac_for_cell(
    state_ds: nc.Dataset,
    state_meta: GridMeta,
    cfg: dict,
    t: int,
    i_internal: int,
    j_internal: int,
    lon_res: float,
) -> Dict[str, float]:
    """
    返回指定时空格点的粗类面积分数: v/s/p/c/U
    """
    mapping_num = get_num_to_luh2_mapping(cfg)
    mapping_lulc = get_luh2_to_lulc_mapping(cfg)

    ii = ds_lat_idx(state_meta, i_internal)
    jj = ds_lon_idx_from_internal_j(j_internal, lon_res, state_meta)

    out = {k: 0.0 for k in COARSE_KEYS}

    for idx in range(1, 13):
        var_name = mapping_num[idx]
        coarse = mapping_lulc[var_name]
        frac = safe_read_scalar(state_ds.variables[var_name], t, ii, jj)
        frac = min(max(frac, 0.0), 1.0)
        out[coarse] += frac

    return out


def get_transition_outgoing_frac_for_cell(
    trans_ds: nc.Dataset,
    trans_meta: GridMeta,
    cfg: dict,
    t: int,
    i_internal: int,
    j_internal: int,
    lon_res: float,
) -> Dict[str, float]:
    """
    统计每个 source 粗类的“普通 transition 总转出 fraction”
    """
    mapping_num = get_num_to_luh2_mapping(cfg)
    mapping_lulc = get_luh2_to_lulc_mapping(cfg)

    ii = ds_lat_idx(trans_meta, i_internal)
    jj = ds_lon_idx_from_internal_j(j_internal, lon_res, trans_meta)

    outgoing = {k: 0.0 for k in COARSE_KEYS}

    for src_idx in range(1, 13):
        for dst_idx in range(1, 13):
            src_name = mapping_num[src_idx]
            dst_name = mapping_num[dst_idx]
            coarse_src = mapping_lulc[src_name]
            coarse_dst = mapping_lulc[dst_name]

            var_name = f"{src_name}_to_{dst_name}_frac"
            var = trans_ds.variables.get(var_name)
            if var is None:
                continue

            frac = safe_read_scalar(var, t, ii, jj)
            frac = min(max(frac, 0.0), 1.0)

            if frac > 0.0 and coarse_src in outgoing:
                outgoing[coarse_src] += frac

    return outgoing


def get_transition_pair_frac_for_cell(
    trans_ds: nc.Dataset,
    trans_meta: GridMeta,
    cfg: dict,
    t: int,
    i_internal: int,
    j_internal: int,
    lon_res: float,
) -> Dict[str, float]:
    """
    返回聚合后的 pair transition:
    v_to_s, v_to_p, ...
    """
    mapping_num = get_num_to_luh2_mapping(cfg)
    mapping_lulc = get_luh2_to_lulc_mapping(cfg)

    ii = ds_lat_idx(trans_meta, i_internal)
    jj = ds_lon_idx_from_internal_j(j_internal, lon_res, trans_meta)

    pair_frac: Dict[str, float] = {}

    for src_idx in range(1, 13):
        for dst_idx in range(1, 13):
            src_name = mapping_num[src_idx]
            dst_name = mapping_num[dst_idx]
            coarse_src = mapping_lulc[src_name]
            coarse_dst = mapping_lulc[dst_name]

            var_name = f"{src_name}_to_{dst_name}_frac"
            var = trans_ds.variables.get(var_name)
            if var is None:
                continue

            frac = safe_read_scalar(var, t, ii, jj)
            frac = min(max(frac, 0.0), 1.0)
            if frac <= 0.0:
                continue

            key = f"{coarse_src}_to_{coarse_dst}"
            pair_frac[key] = pair_frac.get(key, 0.0) + frac

    return pair_frac


def get_wood_harvest_for_cell(
    trans_ds: nc.Dataset,
    trans_meta: GridMeta,
    t: int,
    i_internal: int,
    j_internal: int,
    lon_res: float,
    cell_area_ha: float,
) -> Dict[str, Dict[str, float]]:
    """
    返回:
    {
        'v': {'area_frac':..., 'area_ha':..., 'biomass_tC':...},
        's': {'area_frac':..., 'area_ha':..., 'biomass_tC':...}
    }
    """
    ii = ds_lat_idx(trans_meta, i_internal)
    jj = ds_lon_idx_from_internal_j(j_internal, lon_res, trans_meta)

    family_to_src = {
        "primn": "v",
        "primf": "v",
        "secdf": "s",
        "secdn": "s",
        "secmf": "s",
        "secyf": "s",
        "secnf": "s",
    }

    out: Dict[str, Dict[str, float]] = {}

    for fam, src in family_to_src.items():
        v_harv = trans_ds.variables.get(f"{fam}_harv")
        v_bioh = trans_ds.variables.get(f"{fam}_bioh")

        frac = 0.0
        if v_harv is not None:
            frac = safe_read_scalar(v_harv, t, ii, jj)
            frac = min(max(frac, 0.0), 1.0)

        biomass_tC = 0.0
        if v_bioh is not None:
            bioh_kgC = safe_read_scalar(v_bioh, t, ii, jj)
            biomass_tC = max(0.0, bioh_kgC * 1e-3)

        if frac <= 0.0 and biomass_tC <= 0.0:
            continue

        rec = out.get(src, {"area_frac": 0.0, "area_ha": 0.0, "biomass_tC": 0.0})
        rec["area_frac"] += frac
        rec["area_ha"] += frac * cell_area_ha
        rec["biomass_tC"] += biomass_tC
        out[src] = rec

    return out


# =========================================================
# 5. 诊断逻辑
# =========================================================

def infer_next_state_from_state_plus_transitions(
    state_now: Dict[str, float],
    pair_frac: Dict[str, float],
) -> Dict[str, float]:
    """
    仅用普通 transition 推算 t+1 的粗类 fraction
    harvest 不改变 land cover，因此不进这里
    """
    out = state_now.copy()

    for key, frac in pair_frac.items():
        src, _, dst = key.partition("_to_")
        if src not in out or dst not in out:
            continue
        out[src] -= frac
        out[dst] += frac

    return out


def diagnose_luh2(
    state_path: str,
    trans_path: str,
    config_path: str,
    out_csv: str = "luh2_diagnosis.csv",
    lat_min: float = -90.0,
    lat_max: float = 90.0,
    lon_min: float = -180.0,
    lon_max: float = 180.0,
    start_t: Optional[int] = None,
    end_t: Optional[int] = None,
    excess_tol_frac: float = 1e-8,
    state_mismatch_tol_frac: float = 1e-6,
) -> pd.DataFrame:
    """
    主诊断函数
    """
    cfg = load_config(config_path)

    state_ds = nc.Dataset(state_path, "r")
    trans_ds = nc.Dataset(trans_path, "r")

    try:
        state_meta = load_grid_meta(state_ds)
        trans_meta = load_grid_meta(trans_ds)

        lat_internal, lon_internal, lat_res, lon_res = build_internal_grid_from_states(state_meta)
        area_grid_ha = build_area_grid_ha(lat_internal, lon_internal, lat_res, lon_res)

        time_name_state = "time" if "time" in state_ds.variables else list(state_ds.dimensions.keys())[0]
        time_name_trans = "time" if "time" in trans_ds.variables else list(trans_ds.dimensions.keys())[0]

        ntime_state = len(state_ds.dimensions[time_name_state])
        ntime_trans = len(trans_ds.dimensions[time_name_trans])

        ntime = min(ntime_state, ntime_trans)

        if start_t is None:
            start_t = 0
        if end_t is None:
            end_t = ntime - 2  # 因为还要比较 t+1
        end_t = min(end_t, ntime - 2)

        lat_mask = (lat_internal >= lat_min) & (lat_internal <= lat_max)
        lon_mask = (lon_internal >= lon_min) & (lon_internal <= lon_max)

        i_list = np.where(lat_mask)[0]
        j_list = np.where(lon_mask)[0]

        rows: List[dict] = []

        total_steps = (end_t - start_t + 1) * len(i_list) * len(j_list)
        step = 0
        for t in range(start_t, end_t + 1):
            print(f"\n=== Year index {t} ===")
            for i in i_list:
                for j in j_list:
                    step += 1
                    if step % 100 == 0:
                        pct = step / total_steps * 100
                        print(f"Progress: {pct:.2f}% ({step}/{total_steps})")
                    cell_area_ha = float(area_grid_ha[i, j])
                    if cell_area_ha <= 0.0:
                        continue

                    state_now = get_state_frac_for_cell(
                        state_ds, state_meta, cfg, t, i, j, lon_res
                    )
                    state_next_true = get_state_frac_for_cell(
                        state_ds, state_meta, cfg, t + 1, i, j, lon_res
                    )

                    pair_frac = get_transition_pair_frac_for_cell(
                        trans_ds, trans_meta, cfg, t, i, j, lon_res
                    )
                    outgoing = get_transition_outgoing_frac_for_cell(
                        trans_ds, trans_meta, cfg, t, i, j, lon_res
                    )
                    wood = get_wood_harvest_for_cell(
                        trans_ds, trans_meta, t, i, j, lon_res, cell_area_ha
                    )

                    state_next_pred = infer_next_state_from_state_plus_transitions(state_now, pair_frac)

                    lat_val = float(lat_internal[i])
                    lon_val = float(lon_internal[j])

                    row = {
                        "time_index": t,
                        "lat": lat_val,
                        "lon": lon_val,
                        "cell_area_ha": cell_area_ha,
                    }

                    # 当前面积分数
                    for k in COARSE_KEYS:
                        row[f"state_{k}"] = state_now.get(k, 0.0)

                    # 下一年真实 state
                    for k in COARSE_KEYS:
                        row[f"state_next_true_{k}"] = state_next_true.get(k, 0.0)

                    # 仅由普通 transition 推算的下一年
                    for k in COARSE_KEYS:
                        row[f"state_next_pred_{k}"] = state_next_pred.get(k, 0.0)

                    # 普通 transition 总转出
                    for k in COARSE_KEYS:
                        row[f"trans_out_frac_{k}"] = outgoing.get(k, 0.0)
                        row[f"trans_out_area_ha_{k}"] = outgoing.get(k, 0.0) * cell_area_ha

                    # wood harvest
                    for k in ["v", "s"]:
                        row[f"wood_area_frac_{k}"] = wood.get(k, {}).get("area_frac", 0.0)
                        row[f"wood_area_ha_{k}"] = wood.get(k, {}).get("area_ha", 0.0)
                        row[f"wood_biomass_tC_{k}"] = wood.get(k, {}).get("biomass_tC", 0.0)

                    # 核心检查 1：普通 transition 是否已超过 source state
                    for k in COARSE_KEYS:
                        excess = outgoing.get(k, 0.0) - state_now.get(k, 0.0)
                        row[f"flag_trans_exceed_state_{k}"] = int(excess > excess_tol_frac)
                        row[f"trans_minus_state_frac_{k}"] = excess

                    # 核心检查 2：普通 transition + wood harvest 面积 是否超过 source state
                    for k in ["v", "s"]:
                        combined = outgoing.get(k, 0.0) + wood.get(k, {}).get("area_frac", 0.0)
                        excess = combined - state_now.get(k, 0.0)
                        row[f"combined_out_plus_wood_frac_{k}"] = combined
                        row[f"flag_combined_exceed_state_{k}"] = int(excess > excess_tol_frac)
                        row[f"combined_minus_state_frac_{k}"] = excess

                    # 核心检查 3：普通 transitions 推到下一年，是否和真实 next state 对得上
                    for k in COARSE_KEYS:
                        mismatch = state_next_pred.get(k, 0.0) - state_next_true.get(k, 0.0)
                        row[f"next_state_mismatch_{k}"] = mismatch
                        row[f"flag_next_state_mismatch_{k}"] = int(abs(mismatch) > state_mismatch_tol_frac)

                    # 木材收获强度诊断: biomass / area
                    for k in ["v", "s"]:
                        area_h = wood.get(k, {}).get("area_ha", 0.0)
                        bio_h = wood.get(k, {}).get("biomass_tC", 0.0)
                        row[f"wood_biomass_density_tCha_{k}"] = (bio_h / area_h) if area_h > 0 else np.nan

                    rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(out_csv, index=False, encoding="utf-8-sig")
        return df

    finally:
        state_ds.close()
        trans_ds.close()


# =========================================================
# 6. 汇总打印
# =========================================================

def print_summary(df: pd.DataFrame) -> None:
    print("\n================ DIAGNOSIS SUMMARY ================\n")
    if df.empty:
        print("No records.")
        return

    print(f"Total checked records: {len(df)}")

    for k in COARSE_KEYS:
        col = f"flag_trans_exceed_state_{k}"
        if col in df.columns:
            n = int(df[col].sum())
            if n > 0:
                print(f"[WARN] trans_out > state for {k}: {n} records")

    for k in ["v", "s"]:
        col = f"flag_combined_exceed_state_{k}"
        if col in df.columns:
            n = int(df[col].sum())
            if n > 0:
                print(f"[WARN] transition + wood harvest > state for {k}: {n} records")

    for k in COARSE_KEYS:
        col = f"flag_next_state_mismatch_{k}"
        if col in df.columns:
            n = int(df[col].sum())
            if n > 0:
                print(f"[WARN] predicted next state != true next state for {k}: {n} records")

    print("\nTop suspicious rows for v/s combined exceed:")
    sus_cols = [
        "time_index", "lat", "lon",
        "state_v", "trans_out_frac_v", "wood_area_frac_v", "combined_out_plus_wood_frac_v",
        "state_s", "trans_out_frac_s", "wood_area_frac_s", "combined_out_plus_wood_frac_s",
        "combined_minus_state_frac_v", "combined_minus_state_frac_s",
    ]
    sus_cols = [c for c in sus_cols if c in df.columns]

    tmp = df.copy()
    tmp["max_combined_excess"] = np.nanmax(
        np.vstack([
            tmp["combined_minus_state_frac_v"].to_numpy() if "combined_minus_state_frac_v" in tmp else np.full(len(tmp), np.nan),
            tmp["combined_minus_state_frac_s"].to_numpy() if "combined_minus_state_frac_s" in tmp else np.full(len(tmp), np.nan),
        ]),
        axis=0
    )
    tmp = tmp.sort_values("max_combined_excess", ascending=False)
    print(tmp[sus_cols + ["max_combined_excess"]].head(20).to_string(index=False))


# =========================================================
# 7. 直接运行
# =========================================================

if __name__ == "__main__":
    # ===== 这里改成你的路径 =====
    BASE = Path(r"D:\Work\Research_doc\LULCC\BookKeeping\In_ncfile")
    state_path = BASE / "states_1deg.nc"
    trans_path = BASE / "transitions_1deg.nc"
    config_path = BASE / "config.yml"
    out_csv = BASE / "luh2_diagnosis.csv"

    # 例子：全国
    df = diagnose_luh2(
        state_path=str(state_path),
        trans_path=str(trans_path),
        config_path=str(config_path),
        out_csv=str(out_csv),
        #lat_min=-90.0,
        #lat_max=90.0,
        #lon_min=-180.0,
        #lon_max=180.0,
        lat_min = 60.0,
        lat_max = 65.0, 
        lon_min = 70.0,
        lon_max = 75.0,
        start_t=0,
        end_t=20,   # 先只查前20年，避免太慢
    )
    print_summary(df)
    print(f"\nSaved to: {out_csv}")