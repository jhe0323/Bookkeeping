import netCDF4 as nc
import numpy as np
import csv
import sys

# --- 用户设置 ---
file_path = r"D:\Work\Research_doc\LULCC\BookKeeping\In_ncfile\states.nc"  # LUH2 数据文件路径
lat_idx = 360  # 网格索引或行号
lon_idx = 800  # 网格索引或列号
output_csv = "grid_last_secma_secmb.csv"
n_check_years = 50  # 检测前 n 年是否有有效值

# --- 打开 LUH2 数据 ---
ds = nc.Dataset(file_path)

# --- 读取变量 ---
secma = ds.variables['secma'][:]  # shape: (time, lat, lon)
secmb = ds.variables['secmb'][:]
time = ds.variables['time'][:]    # 时间序列

# --- 取指定网格的时间序列 ---
secma_series = secma[:, lat_idx, lon_idx]
secmb_series = secmb[:, lat_idx, lon_idx]

# --- 检测前 n 年是否有有效值 ---
# 这里认为有效值 > 0
if not (np.any(secma_series[:n_check_years] > 0) and np.any(secmb_series[:n_check_years] > 0)):
    print(f"网格 ({lat_idx}, {lon_idx}) 前 {n_check_years} 年没有有效 secma/secmb 值，程序退出")
    ds.close()
    sys.exit(1)

# --- 保存整个时间序列到 CSV ---
with open(output_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['year', 'secma', 'secmb'])  # 表头
    for t, s_ma, s_mb in zip(time, secma_series, secmb_series):
        writer.writerow([t, s_ma, s_mb])

print(f"网格 ({lat_idx}, {lon_idx}) 的时间序列已保存到 {output_csv}")

ds.close()