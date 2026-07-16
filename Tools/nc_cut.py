import netCDF4 as nc
import numpy as np


def average_nc(input_file, output_file, start_year, end_year):
    ds_in = nc.Dataset(input_file, "r")

    # ---- 找时间变量 ----
    time_var = None
    for name in ds_in.variables:
        if name.lower() in ["time", "year"]:
            time_var = name
            break

    if time_var is None:
        raise ValueError("No time variable found.")

    time_data = ds_in.variables[time_var][:]

    # ---- 找对应年份索引 ----
    idx = np.where((time_data >= start_year) & (time_data <= end_year))[0]

    if len(idx) == 0:
        raise ValueError("No data in selected year range.")

    print(f"Selected years: {time_data[idx][0]} - {time_data[idx][-1]}")
    print(f"Number of time steps: {len(idx)}")

    # ---- 创建输出文件 ----
    ds_out = nc.Dataset(output_file, "w")

    # ---- 复制维度 ----
    for dim_name, dim in ds_in.dimensions.items():
        if dim_name == time_var:
            ds_out.createDimension(dim_name, 1)
        else:
            ds_out.createDimension(dim_name, len(dim))

    # ---- 复制变量 ----
    for var_name, var in ds_in.variables.items():
        dims = var.dimensions
        dtype = var.datatype

        fill_value = None
        if "_FillValue" in var.ncattrs():
            fill_value = var.getncattr("_FillValue")

        # 创建变量时就传入
        if fill_value is not None:
            out_var = ds_out.createVariable(var_name, dtype, dims, fill_value=fill_value)
        else:
            out_var = ds_out.createVariable(var_name, dtype, dims)

        # ---- 再复制其他属性（跳过 _FillValue）----
        for attr in var.ncattrs():
            if attr == "_FillValue":
                continue
            out_var.setncattr(attr, var.getncattr(attr))

        data = var[:]

        # ---- 如果有时间维，做平均 ----
        if time_var in dims:
            axis = dims.index(time_var)
            mean_data = np.nanmean(data[idx, ...], axis=axis)

            # 保持维度一致（time=1）
            mean_data = np.expand_dims(mean_data, axis=axis)

            out_var[:] = mean_data
        else:
            out_var[:] = data

    # ---- 更新时间变量 ----
    ds_out.variables[time_var][:] = np.array([int((start_year + end_year) / 2)])

    ds_in.close()
    ds_out.close()

    print(f"Saved to: {output_file}")


# ===== 使用示例 =====
if __name__ == "__main__":
    input_nc = r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\summary_025deg.global_time_gt1100.nc"
    output_nc = r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\summary_025deg.global_time_gt1100_cut_1160-1170.nc"

    start_year = 1160
    end_year = 1170

    average_nc(input_nc, output_nc, start_year, end_year)