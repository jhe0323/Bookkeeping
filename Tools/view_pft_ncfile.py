from netCDF4 import Dataset
import numpy as np

nc_path = r"D:/Work/Research_doc/LULCC/BookKeeping/In_ncfile/PFTmap_IBIS_1deg.nc"

with Dataset(nc_path, "r") as ds:
    # 这里的 frac 形状应该是 (1, 1, 360, 720)
    frac = ds.variables["vegtype"][:] 
    
    # 转换为普通 numpy 数组并处理掩码
    if np.ma.isMaskedArray(frac):
        data = frac.filled(np.nan)
    else:
        data = frac

    # 获取所有不重复的数值（类别索引）
    # 排除掉填充值和 NaN
    valid_data = data[~np.isnan(data)]
    vtypes = np.unique(valid_data)
    
    print("--- 文件数值统计 ---")
    print(f"数据形状 (time, level, lat, lon): {data.shape}")
    print(f"出现的唯一值（类别索引）: \n{vtypes}")
    print(f"最小值: {np.nanmin(data)}")
    print(f"最大值: {np.nanmax(data)}")