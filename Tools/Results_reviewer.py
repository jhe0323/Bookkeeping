import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import numpy as np

def plot_lulcc_spatial_mean(file_path, variable_name, start_year=None, end_year=None):
    """
    查看指定参数在特定时间段内的空间分布均值。
    逻辑：负值越小越绿，正值越大越红，0值为白色。
    """
    # 1. 加载数据
    ds = xr.open_dataset(file_path)
    
    if variable_name not in ds.variables:
        print(f"错误: 变量 '{variable_name}' 不在文件中。")
        return

    # 2. 时间切片与均值计算
    if start_year is None: start_year = int(ds.time.min())
    if end_year is None: end_year = int(ds.time.max())
    
    data_slice = ds[variable_name].sel(time=slice(start_year, end_year))
    spatial_mean = data_slice.mean(dim='time')

    # 3. 极值计算
    d_min = float(spatial_mean.min())
    d_max = float(spatial_mean.max())

    # 4. 核心逻辑：动态裁剪色带并设置 Norm
    # 原始全色带：绿色(负) -> 白色(0) -> 红色(正)
    full_colors = ["green", "white", "red"]
    
    if d_min < 0 < d_max:
        # 跨越 0 值：使用完整色带，强制 0 为白色
        custom_cmap = mcolors.LinearSegmentedColormap.from_list("GnWhRd", full_colors)
        norm = mcolors.TwoSlopeNorm(vmin=d_min, vcenter=0, vmax=d_max)
    elif d_max <= 0:
        # 纯负值（如 Gross_Sinks）：只需 绿色 -> 白色 部分
        custom_cmap = mcolors.LinearSegmentedColormap.from_list("GnWh", ["green", "white"])
        # 即使没有正值，也将 vmax 设为 0，确保最接近 0 的值是白色而不是绿色
        norm = mcolors.Normalize(vmin=d_min, vmax=0)
    else:
        # 纯正值：只需 白色 -> 红色 部分
        custom_cmap = mcolors.LinearSegmentedColormap.from_list("WhRd", ["white", "red"])
        # 即使没有负值，也将 vmin 设为 0，确保最接近 0 的值是白色
        norm = mcolors.Normalize(vmin=0, vmax=d_max)

    # 5. 绘图
    plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    
    units = ds[variable_name].attrs.get('units', 'n/a')
    
    img = spatial_mean.plot(
        ax=ax, 
        transform=ccrs.PlateCarree(),
        cmap=custom_cmap, 
        norm=norm,
        add_colorbar=True,
        cbar_kwargs={
            'label': f'Mean {variable_name} ({units})',
            'extend': 'both'
        }
    )

    plt.title(f"Spatial Mean: {variable_name}\n({start_year}-{end_year}) | 0 is White")
    plt.show()

# 调用测试
plot_lulcc_spatial_mean("summary_1deg.global_1950_4.28.nc", "Net_Emissions", start_year=1110, end_year=1120)