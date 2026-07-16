import netCDF4 as nc
import numpy as np
import pandas as pd

def summarize_summary_nc(nc_path: str, out_csv: str):
    """
    读取 summary.nc 文件，将每一年所有网格的输出结果按类型汇总求和。
    生成每年总量表（单位与原文件一致）。
    """
    # 打开 NetCDF 文件
    ds = nc.Dataset(nc_path)

    # 读取变量
    years = ds.variables["time"][:]
    biomass = ds.variables["biomass_total"][:]
    soil    = ds.variables["soil_total"][:]
    P1      = ds.variables["P1"][:]
    P10     = ds.variables["P10"][:]
    P100    = ds.variables["P100"][:]
    atmosphere    = ds.variables["atmosphere"][:]

    # 检查形状
    print(f"数据维度: time={biomass.shape[0]}, lat={biomass.shape[1]}, lon={biomass.shape[2]}")

    # 汇总每年所有格点
    totals = []
    for t in range(len(years)):
        total_biomass = np.nansum(biomass[t, :, :])
        total_soil    = np.nansum(soil[t, :, :])
        total_P1      = np.nansum(P1[t, :, :])
        total_P10     = np.nansum(P10[t, :, :])
        total_P100    = np.nansum(P100[t, :, :])
        total_atmosphere    = np.nansum(atmosphere[t, :, :])

        totals.append({
            "Year": int(years[t]),
            "Biomass_total": total_biomass,
            "Soil_total": total_soil,
            "P1_total": total_P1,
            "P10_total": total_P10,
            "P100_total": total_P100,
            "Atmo_total": total_atmosphere
        })

    df = pd.DataFrame(totals)
    df.to_csv(out_csv, index=False)
    print(f"✅ 汇总完成，结果已保存到 {out_csv}")
    print(df.head())

    ds.close()
    return df

if __name__ == "__main__":
    summarize_summary_nc("summary_1deg_2.6.nc", "summary_1deg_2.6.csv")
