import netCDF4 as nc
import numpy as np
import pandas as pd
import time


def summarize_summary_nc(nc_path: str, out_csv: str, progress_every: int = 10):
    """
    汇总指定变量。
    对每一年所有网格求和，并显示计算进度。

    Parameters
    ----------
    nc_path : str
        输入 nc 文件路径
    out_csv : str
        输出 csv 文件路径
    progress_every : int
        每隔多少个年份打印一次进度
    """

    vars_to_read = [
        "Gross_Sources",
        "Gross_Sinks",
        "Net_Emissions",
        "Flux_Harvest_Net",
        "Flux_WHp",
        "Flux_Harvest_SoilSlash",
        "Flux_Harvest_Regrowth",
        "Flux_Clearing",
        "Flux_Abandonment",
        "Flux_Other",
        "harvest_requested_biomass",
        "harvest_met_biomass",
        "harvest_unmet_raw_biomass",
        "harvest_forced_biomass",
        "harvest_to_products",
        "harvest_to_SR_from_biomass",
        "harvest_to_SR_from_soil",
        "harvest_delta_B_h",
        "harvest_delta_SS_h",
    ]

    start_time = time.time()

    with nc.Dataset(nc_path) as ds:
        years = ds.variables["time"][:]
        existing_vars = list(ds.variables.keys())

        missing_vars = [v for v in vars_to_read if v not in existing_vars]
        vars_to_read = [v for v in vars_to_read if v in existing_vars]

        if missing_vars:
            print("\n以下变量在 nc 中不存在，已跳过：")
            for v in missing_vars:
                print("  ", v)

        if not vars_to_read:
            raise ValueError("没有检测到需要汇总的变量。")

        print("\n检测到并汇总的变量：")
        for v in vars_to_read:
            print("  ", v, ds.variables[v].dimensions, ds.variables[v].shape)

        n_years = len(years)
        n_vars = len(vars_to_read)
        total_tasks = n_years * n_vars

        print("\n开始汇总：")
        print(f"  年份数量: {n_years}")
        print(f"  变量数量: {n_vars}")
        print(f"  总任务数: {total_tasks}")

        totals = []

        for t in range(n_years):
            year_start = time.time()
            row = {"Year": int(years[t])}

            for v in vars_to_read:
                var = ds.variables[v]
                dims = var.dimensions

                if "time" not in dims:
                    continue

                # 构造只读取当前年份的切片，避免每次 var[:] 读取整个变量
                index = []
                for d in dims:
                    if d == "time":
                        index.append(t)
                    else:
                        index.append(slice(None))

                arr_t = var[tuple(index)]
                arr_t = np.ma.filled(arr_t, np.nan)

                row[f"{v}_total"] = float(np.nansum(arr_t))

            totals.append(row)

            finished_years = t + 1
            finished_tasks = finished_years * n_vars
            percent = finished_tasks / total_tasks * 100.0

            elapsed = time.time() - start_time
            avg_per_year = elapsed / finished_years
            remain_years = n_years - finished_years
            eta = avg_per_year * remain_years

            if (
                finished_years == 1
                or finished_years % progress_every == 0
                or finished_years == n_years
            ):
                print(
                    f"进度: {finished_years}/{n_years} 年 "
                    f"({percent:.1f}%) | "
                    f"当前年份: {int(years[t])} | "
                    f"本年耗时: {time.time() - year_start:.2f}s | "
                    f"已用: {elapsed / 60:.1f} min | "
                    f"预计剩余: {eta / 60:.1f} min"
                )

    df = pd.DataFrame(totals)

    df.to_csv(out_csv, index=False, encoding="utf-8-sig")

    total_elapsed = time.time() - start_time

    print(f"\n✅ 汇总完成，结果已保存到 {out_csv}")
    print(f"总耗时: {total_elapsed / 60:.2f} min")
    print(df.head())

    return df


if __name__ == "__main__":
    summarize_summary_nc(
        "summary_025deg.global_time_ge1000.nc",
        "summary_025deg.global_time_ge1000_GrossNet.csv",
        progress_every=50,
    )