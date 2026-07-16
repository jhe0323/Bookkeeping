import netCDF4 as nc
import numpy as np
import pandas as pd

def list_transition_vars(nc_path, prefix=None, include=None, exclude=None):
    """
    列出 transitions.nc 中符合条件的变量名。
    - prefix: 仅保留以 prefix 开头的变量
    - include: 仅保留名称包含 include 子串的变量（可选）
    - exclude: 排除名称包含 exclude 子串的变量（可选）
    """
    ds = nc.Dataset(nc_path, "r")
    var_names = list(ds.variables.keys())
    ds.close()

    # 排除坐标变量
    var_names = [v for v in var_names if v not in ("time", "lat", "lon")]

    if prefix is not None:
        var_names = [v for v in var_names if v.startswith(prefix)]
    if include is not None:
        var_names = [v for v in var_names if include in v]
    if exclude is not None:
        var_names = [v for v in var_names if exclude not in v]

    return sorted(var_names)


def global_sum_transition(
    nc_path,
    var_name,
    start_time_idx=0,
    end_time_idx=None,
    time_var="time"
):
    """
    单变量：逐年(逐time索引)对全球(lat,lon)做 nansum，返回 DataFrame
    """
    ds = nc.Dataset(nc_path, "r")
    var = ds.variables[var_name]  # (time, lat, lon)
    time = ds.variables[time_var][:]

    if end_time_idx is None:
        end_time_idx = var.shape[0]

    results = []
    for ti in range(start_time_idx, end_time_idx):
        total = np.nansum(var[ti, :, :])
        results.append({
            "time_index": ti,
            "time_value": float(time[ti]),
            var_name: float(total)
        })

    ds.close()
    return pd.DataFrame(results)


def global_sum_transitions_by_prefix(
    nc_path,
    prefix,
    start_time_idx=0,
    end_time_idx=None,
    time_var="time",
    include=None,
    exclude=None,
    chunk_years=1
):
    """
    prefix 多变量：输出一个表，行=年份(时间索引)，列=变量名，值=该变量的全球 fraction 和。
    - chunk_years: 每次读取的时间步数（防止内存爆炸）。默认1最稳。
    """
    var_list = list_transition_vars(nc_path, prefix=prefix, include=include, exclude=exclude)
    if len(var_list) == 0:
        raise ValueError(f"No variables found with prefix='{prefix}' (include={include}, exclude={exclude}).")

    ds = nc.Dataset(nc_path, "r")
    time = ds.variables[time_var][:]
    ntime = ds.variables[var_list[0]].shape[0]

    if end_time_idx is None:
        end_time_idx = ntime

    # 输出表的时间轴
    tids = list(range(start_time_idx, end_time_idx))
    out = pd.DataFrame({
        "time_index": tids,
        "time_value": [float(time[ti]) for ti in tids]
    })

    # 逐变量计算（更稳健，内存占用小）
    for v in var_list:
        var = ds.variables[v]
        vals = []

        # 分块读 time 维
        for t0 in range(start_time_idx, end_time_idx, chunk_years):
            t1 = min(end_time_idx, t0 + chunk_years)
            block = var[t0:t1, :, :]  # shape=(chunk, lat, lon)
            # 对每个时间步求和
            block_sum = np.nansum(block, axis=(1, 2))
            vals.extend([float(x) for x in block_sum])

        out[v] = vals

    ds.close()
    return out


def run_transition_summary(
    nc_path,
    mode="single",
    var_name=None,
    prefix=None,
    start_time_idx=0,
    end_time_idx=None,
    out_csv="transition_global_sum.csv",
    include=None,
    exclude=None
):
    """
    统一入口：
    - mode="single": 需要 var_name
    - mode="prefix": 需要 prefix
    """
    if mode == "single":
        if not var_name:
            raise ValueError("mode='single' requires var_name.")
        df = global_sum_transition(
            nc_path, var_name=var_name,
            start_time_idx=start_time_idx, end_time_idx=end_time_idx
        )
    elif mode == "prefix":
        if not prefix:
            raise ValueError("mode='prefix' requires prefix.")
        df = global_sum_transitions_by_prefix(
            nc_path, prefix=prefix,
            start_time_idx=start_time_idx, end_time_idx=end_time_idx,
            include=include, exclude=exclude,
            chunk_years=1
        )
    else:
        raise ValueError("mode must be 'single' or 'prefix'.")

    df.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}  (rows={len(df)}, cols={df.shape[1]})")
    return df
    
# =========================
# 示例用法
# =========================
if __name__ == "__main__":
    nc_file = "transitions_1deg.nc"
    start_time_idx = 1000
    run_transition_summary(
        nc_path=nc_file,
        mode="prefix",
        prefix="primn_",
        start_time_idx=start_time_idx,
        out_csv="global_primn_prefix_fraction.csv"
    )


    # # 1) 单变量：primn_to_secdf
    # run_transition_summary(
        # nc_path=nc_file,
        # mode="single",
        # var_name="primn_to_secdf",
        # start_time_idx=start_time_idx,
        # out_csv="global_primn_to_secdf_fraction.csv"
    # )

    # # 2) prefix 多变量：比如所有 primn_*（包括 primn_to_secdf / primn_to_...）
    # run_transition_summary(
        # nc_path=nc_file,
        # mode="prefix",
        # prefix="primn_",
        # start_time_idx=start_time_idx,
        # out_csv="global_primn_prefix_fraction.csv"
    # )

    # # 3) 更精细筛选：所有 prim* 但只要包含 "_to_secdf"
    # run_transition_summary(
        # nc_path=nc_file,
        # mode="prefix",
        # prefix="prim",
        # include="_to_secdf",
        # start_time_idx=start_time_idx,
        # out_csv="global_prim_to_secdf_fraction.csv"
    # )