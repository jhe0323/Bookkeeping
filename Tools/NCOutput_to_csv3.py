import netCDF4 as nc
import numpy as np
import pandas as pd


def summarize_summary_nc(nc_path: str, out_csv: str, add_diff: bool = True):
    """
    读取 summary.nc 文件，将每一年所有网格的输出结果按变量汇总求和。
    自动适配当前 summary.nc 中存在的新增 diag 变量。

    Parameters
    ----------
    nc_path : str
        summary.nc 路径
    out_csv : str
        输出 csv 路径
    add_diff : bool
        是否增加逐年差值列（t - t-1）

    Returns
    -------
    pd.DataFrame
    """
    ds = nc.Dataset(nc_path)

    try:
        years = ds.variables["time"][:]

        # 必选核心变量
        core_vars = [
            "biomass_total",
            "soil_total",
            "P1",
            "P10",
            "P100",
            "atmosphere",
            "unmet_harvest",
        ]

        # 当前代码里可能存在的过程诊断变量
        diag_vars = [
            "area_clearing",
            "area_abandonment",
            "area_other",
            "area_harvest",
            "clearing_to_A",
            "clearing_to_products",
            "clearing_to_soil",
            "harvest_to_products",
            "harvest_to_soil",
            "harvest_biomass_removed",
            "emit_clearing",
            "emit_abandonment",
            "emit_other",
            "emit_harvest",
            "emit_products",
            "emit_abandonment_biomass",
            "emit_abandonment_soil",
            "area_deforestation",
            "deforestation_to_A",
            "deforestation_to_products",
            "deforestation_to_soil",
            "emit_products_clearing",
            "emit_products_harvest",
            "harvest_requested_biomass",
            "harvest_met_biomass", "harvest_unmet_biomass", "harvest_beta_delta_sum",
            "harvest_delta_B_h", "harvest_R", "harvest_delta_SS_h","harvest_unmet_to_products", "harvest_unmet_to_soil",
        ]

        existing_vars = list(ds.variables.keys())
        vars_to_read = [v for v in core_vars + diag_vars if v in existing_vars]

        print("检测到的变量：")
        for v in vars_to_read:
            print("  ", v)

        # 读取所有存在的变量
        data = {v: ds.variables[v][:] for v in vars_to_read}

        # 检查一个代表变量形状
        ref_name = "biomass_total" if "biomass_total" in data else vars_to_read[0]
        ref = data[ref_name]
        print(
            f"\n数据维度: time={ref.shape[0]}, lat={ref.shape[1]}, lon={ref.shape[2]}"
        )

        totals = []
        for t in range(len(years)):
            row = {"Year": int(years[t])}

            # 对所有变量做全局求和
            for v in vars_to_read:
                row[f"{v}_total"] = float(np.nansum(data[v][t, :, :]))

            # 派生核心总量
            biomass_total = row.get("biomass_total_total", 0.0)
            soil_total = row.get("soil_total_total", 0.0)
            p1_total = row.get("P1_total", 0.0)
            p10_total = row.get("P10_total", 0.0)
            p100_total = row.get("P100_total", 0.0)
            atmosphere_total = row.get("atmosphere_total", 0.0)

            c_sys_total = biomass_total + soil_total + p1_total + p10_total + p100_total
            c_bookkeep_total = c_sys_total + atmosphere_total

            row["C_sys_total"] = c_sys_total
            row["C_bookkeep_total"] = c_bookkeep_total

            totals.append(row)

        df = pd.DataFrame(totals)

        # 增加逐年差值列
        if add_diff:
            diff_cols = []
            for col in df.columns:
                if col == "Year":
                    continue
                diff_col = f"{col}_diff"
                df[diff_col] = df[col].diff()
                diff_cols.append(diff_col)

            # 更直观地给 atmosphere 单独起一个名字
            if "atmosphere_total" in df.columns:
                df["ELUC_from_atmosphere_diff"] = df["atmosphere_total"].diff()

            if "C_sys_total" in df.columns:
                df["C_sys_change"] = df["C_sys_total"].diff()

            if "C_bookkeep_total" in df.columns:
                df["C_bookkeep_change"] = df["C_bookkeep_total"].diff()

        df.to_csv(out_csv, index=False, encoding="utf-8-sig")

        print(f"\n✅ 汇总完成，结果已保存到 {out_csv}")
        print(df.head())

        return df

    finally:
        ds.close()


if __name__ == "__main__":
    summarize_summary_nc(
        "summary_1deg.global_1950_4.25.nc",
        "summary_1deg.global_1950_4.25.csv",
        add_diff=True,
    )