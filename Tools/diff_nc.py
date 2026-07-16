import xarray as xr

# 读取全局或分区结果
ds = xr.open_dataset("summary_1deg_1107.global.nc")  # 或单个 rank 文件

# 沿时间维度做差分（后一时刻 - 前一时刻）
diff_ds = ds.diff(dim="time", label="upper")

# 保存结果
diff_ds.to_netcdf("summary_1deg_1107.global.diff.nc")
print("✅ 已保存差值文件 summary_1deg.global.diff.nc")
