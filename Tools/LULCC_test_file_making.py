#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import numpy as np
import netCDF4 as nc


# =========================================================
# User settings
# =========================================================
OUT_DIR = Path("./synthetic_case_current_model")

# 模拟长度
N_YEARS = 100
START_YEAR = 0

# 目标测试格点：选在中国区域内，便于你后面 bbox 或单格点测试
TARGET_LAT = 35.5
TARGET_LON = 105.5

# 初始 land-cover fraction
# 这里按照你当前 config.yml 的映射：
#   secdf -> s
#   pastr -> c
S0 = 0.70   # secondary fraction
C0 = 0.30   # cropland/pasture fraction

# 两次双向转换事件
EVENTS = [
    {"year": 20, "secdf_to_pastr": 0.10, "pastr_to_secdf": 0.04},
    {"year": 50, "secdf_to_pastr": 0.03, "pastr_to_secdf": 0.08},
]

# 是否生成所有 12 个 state 变量
WRITE_ALL_STATE_VARS = True


# =========================================================
# Helpers
# =========================================================
STATE_VAR_NAMES = [
    "primf", "primn", "secdf", "secdn", "urban",
    "c3ann", "c4ann", "c3per", "c4per", "c3nfx",
    "pastr", "range"
]

# 只需要这两个 transition 变量就能驱动你要的 s <-> c 测试
TRANS_VAR_NAMES = [
    "secdf_to_pastr_frac",
    "pastr_to_secdf_frac",
]

PFT_VAR_NAME = "maxvegetfrac"


def build_global_1deg_grid():
    """
    生成与你当前代码兼容的全球1°中心点网格：
    lat: -89.5 ... 89.5  (180)
    lon: -179.5 ... 179.5 (360)
    """
    lat = np.arange(-89.5, 90.0, 1.0, dtype=np.float32)   # 180
    lon = np.arange(-179.5, 180.0, 1.0, dtype=np.float32) # 360
    return lat, lon


def find_nearest_index(coord, values):
    values = np.asarray(values)
    return int(np.argmin(np.abs(values - coord)))


def build_landcover_time_series(n_years, s0, c0, events):
    """
    构造逐年状态：
    state[t] 表示第 t 年开始时的 land-cover fraction
    transition[t] 在第 t 年发生，用于推进到 state[t+1]
    """
    if not np.isclose(s0 + c0, 1.0):
        raise ValueError(f"S0 + C0 must equal 1.0, got {s0 + c0}")

    secdf = np.zeros(n_years, dtype=np.float64)
    pastr = np.zeros(n_years, dtype=np.float64)
    secdf_to_pastr = np.zeros(n_years, dtype=np.float64)
    pastr_to_secdf = np.zeros(n_years, dtype=np.float64)

    event_map = {e["year"]: e for e in events}

    secdf[0] = s0
    pastr[0] = c0

    for t in range(n_years):
        if t in event_map:
            e = event_map[t]
            s2c = float(e.get("secdf_to_pastr", 0.0))
            c2s = float(e.get("pastr_to_secdf", 0.0))

            if s2c > secdf[t] + 1e-12:
                raise ValueError(
                    f"Year {t}: secdf_to_pastr={s2c} exceeds secdf={secdf[t]}"
                )
            if c2s > pastr[t] + 1e-12:
                raise ValueError(
                    f"Year {t}: pastr_to_secdf={c2s} exceeds pastr={pastr[t]}"
                )

            secdf_to_pastr[t] = s2c
            pastr_to_secdf[t] = c2s

        if t < n_years - 1:
            secdf[t + 1] = secdf[t] - secdf_to_pastr[t] + pastr_to_secdf[t]
            pastr[t + 1] = pastr[t] + secdf_to_pastr[t] - pastr_to_secdf[t]

            if secdf[t + 1] < -1e-12 or pastr[t + 1] < -1e-12:
                raise ValueError(
                    f"Negative fraction at year {t+1}: secdf={secdf[t+1]}, pastr={pastr[t+1]}"
                )
            if not np.isclose(secdf[t + 1] + pastr[t + 1], 1.0, atol=1e-10):
                raise ValueError(
                    f"Fractions do not sum to 1 at year {t+1}: "
                    f"{secdf[t + 1] + pastr[t + 1]}"
                )

    return secdf, pastr, secdf_to_pastr, pastr_to_secdf


def create_nc_dims(ds, n_time, lat, lon):
    ds.createDimension("time", n_time)
    ds.createDimension("lat", len(lat))
    ds.createDimension("lon", len(lon))


def write_coords(ds, years, lat, lon):
    vtime = ds.createVariable("time", "i4", ("time",))
    vlat = ds.createVariable("lat", "f4", ("lat",))
    vlon = ds.createVariable("lon", "f4", ("lon",))

    vtime[:] = years
    vlat[:] = lat
    vlon[:] = lon

    vtime.long_name = "simulation year index"
    vlat.units = "degrees_north"
    vlon.units = "degrees_east"


def make_states_nc(out_path, years, lat, lon, target_i, target_j, secdf_ts, pastr_ts):
    """
    生成 states 文件。
    你的 _get_initial_fractions() 会从这里读取 12 类 LUH2 state。
    这里只给 secdf/pastr 写非零值，其它变量全零。
    """
    with nc.Dataset(out_path, "w", format="NETCDF4") as ds:
        create_nc_dims(ds, len(years), lat, lon)
        write_coords(ds, years, lat, lon)

        ds.title = "Synthetic states for current LULCC model"
        ds.note = (
            "Global 1deg sparse synthetic dataset. "
            "Only one target cell has nonzero secdf/pastr fractions."
        )

        if WRITE_ALL_STATE_VARS:
            var_names = STATE_VAR_NAMES
        else:
            var_names = ["secdf", "pastr"]

        vars_dict = {}
        for name in var_names:
            v = ds.createVariable(
                name,
                "f4",
                ("time", "lat", "lon"),
                zlib=True,
                complevel=4,
                fill_value=0.0
            )
            v.units = "1"
            v.long_name = f"synthetic state fraction for {name}"
            vars_dict[name] = v

        # 只在目标格点写 secdf / pastr
        vars_dict["secdf"][:, target_i, target_j] = secdf_ts.astype(np.float32)
        vars_dict["pastr"][:, target_i, target_j] = pastr_ts.astype(np.float32)


def make_transitions_nc(
    out_path,
    years,
    lat,
    lon,
    target_i,
    target_j,
    secdf_to_pastr_ts,
    pastr_to_secdf_ts
):
    """
    生成 transitions 文件。
    你的 _parse_transitions() 只会找 xxx_to_xxx_frac 变量。
    """
    with nc.Dataset(out_path, "w", format="NETCDF4") as ds:
        create_nc_dims(ds, len(years), lat, lon)
        write_coords(ds, years, lat, lon)

        ds.title = "Synthetic transitions for current LULCC model"
        ds.note = (
            "Global 1deg sparse synthetic dataset. "
            "Only one target cell has nonzero secdf<->pastr transitions."
        )

        v_s2c = ds.createVariable(
            "secdf_to_pastr_frac",
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            fill_value=0.0
        )
        v_c2s = ds.createVariable(
            "pastr_to_secdf_frac",
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            fill_value=0.0
        )

        v_s2c.units = "1"
        v_c2s.units = "1"
        v_s2c.long_name = "synthetic transition fraction from secdf to pastr"
        v_c2s.long_name = "synthetic transition fraction from pastr to secdf"

        v_s2c[:, target_i, target_j] = secdf_to_pastr_ts.astype(np.float32)
        v_c2s[:, target_i, target_j] = pastr_to_secdf_ts.astype(np.float32)


def make_pft_nc(out_path, lat, lon, target_i, target_j):
    """
    生成 PFT 文件：
    - 变量名必须是 maxvegetfrac
    - 维度做成 (time_counter, veget, lat, lon)
    - 15个PFT
    - 仅目标格点 PFT1=1，其余全0
    """
    with nc.Dataset(out_path, "w", format="NETCDF4") as ds:
        ds.createDimension("time_counter", 1)
        ds.createDimension("veget", 15)
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))

        vtime = ds.createVariable("time_counter", "i4", ("time_counter",))
        vlat = ds.createVariable("lat", "f4", ("lat",))
        vlon = ds.createVariable("lon", "f4", ("lon",))
        vpft = ds.createVariable(
            PFT_VAR_NAME,
            "f4",
            ("time_counter", "veget", "lat", "lon"),
            zlib=True,
            complevel=4,
            fill_value=0.0
        )

        vtime[:] = np.array([0], dtype=np.int32)
        vlat[:] = lat
        vlon[:] = lon

        vlat.units = "degrees_north"
        vlon.units = "degrees_east"
        vpft.units = "1"
        vpft.long_name = "maximum vegetation fraction"

        ds.title = "Synthetic PFT map for current LULCC model"
        ds.note = (
            "Global 1deg sparse synthetic PFT map. "
            "At the target cell only veget-1 has fraction 1."
        )

        # 只在目标格点设置 PFT1=1
        vpft[0, 0, target_i, target_j] = 1.0


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    lat, lon = build_global_1deg_grid()
    years = np.arange(START_YEAR, START_YEAR + N_YEARS, dtype=np.int32)

    target_i = find_nearest_index(TARGET_LAT, lat)
    target_j = find_nearest_index(TARGET_LON, lon)

    secdf_ts, pastr_ts, secdf_to_pastr_ts, pastr_to_secdf_ts = build_landcover_time_series(
        n_years=N_YEARS,
        s0=S0,
        c0=C0,
        events=EVENTS
    )

    states_path = OUT_DIR / "synthetic_states.nc"
    trans_path = OUT_DIR / "synthetic_transitions.nc"
    pft_path = OUT_DIR / "synthetic_pft.nc"

    make_states_nc(
        out_path=states_path,
        years=years,
        lat=lat,
        lon=lon,
        target_i=target_i,
        target_j=target_j,
        secdf_ts=secdf_ts,
        pastr_ts=pastr_ts
    )

    make_transitions_nc(
        out_path=trans_path,
        years=years,
        lat=lat,
        lon=lon,
        target_i=target_i,
        target_j=target_j,
        secdf_to_pastr_ts=secdf_to_pastr_ts,
        pastr_to_secdf_ts=pastr_to_secdf_ts
    )

    make_pft_nc(
        out_path=pft_path,
        lat=lat,
        lon=lon,
        target_i=target_i,
        target_j=target_j
    )

    print("Done.")
    print(f"states      : {states_path}")
    print(f"transitions : {trans_path}")
    print(f"pft         : {pft_path}")
    print()
    print("Target test cell:")
    print(f"  lat index = {target_i}, lat = {lat[target_i]}")
    print(f"  lon index = {target_j}, lon = {lon[target_j]}")
    print()
    print("Quick check:")
    print(f"  year 0  : secdf={secdf_ts[0]:.3f}, pastr={pastr_ts[0]:.3f}")
    print(f"  year 20 event: secdf->pastr={secdf_to_pastr_ts[20]:.3f}, "
          f"pastr->secdf={pastr_to_secdf_ts[20]:.3f}")
    print(f"  year 21 : secdf={secdf_ts[21]:.3f}, pastr={pastr_ts[21]:.3f}")
    print(f"  year 50 event: secdf->pastr={secdf_to_pastr_ts[50]:.3f}, "
          f"pastr->secdf={pastr_to_secdf_ts[50]:.3f}")
    print(f"  year 51 : secdf={secdf_ts[51]:.3f}, pastr={pastr_ts[51]:.3f}")
    print(f"  final   : secdf={secdf_ts[-1]:.3f}, pastr={pastr_ts[-1]:.3f}")


if __name__ == "__main__":
    main()