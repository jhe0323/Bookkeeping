import xarray as xr
import glob

files = sorted(glob.glob("summary_1deg.rank*.nc"))
ds_list = [xr.open_dataset(f) for f in files]

# 按每个文件的 lon 最小值排序，而不是按文件名排序
ds_list = sorted(ds_list, key=lambda d: float(d["lon"].min()))

ds = xr.concat(ds_list, dim="lon")

ds.to_netcdf("summary_1deg.global.nc")
print("merged:", len(files), "files")