import numpy as np
import xarray as xr
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds
from rasterio.crs import CRS

# --- 标准 WGS84 WKT，避免 Windows/conda 的 proj.db 问题 ---
WGS84_WKT = (
'GEOGCRS["WGS 84",'
' DATUM["World Geodetic System 1984",'
'  ELLIPSOID["WGS 84",6378137,298.257223563,'
'   LENGTHUNIT["metre",1]]],'
' PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],'
' CS[ellipsoidal,2],'
'  AXIS["latitude",north],'
'  AXIS["longitude",east],'
'  ANGLEUNIT["degree",0.0174532925199433],'
' ID["EPSG",4326]]'
)
crs_wgs84 = CRS.from_wkt(WGS84_WKT)

def edges_from_centers(centers: np.ndarray) -> np.ndarray:
    """由中心点坐标推算像元边界：内部边界=相邻中心中点，首尾外推半个格距。"""
    centers = np.asarray(centers, dtype=float)
    d = np.diff(centers)
    if centers.size < 2 or not np.all(np.isfinite(d)):
        raise ValueError("centers 至少需要两个且为有限数。")
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0]  = centers[0]  - 0.5 * d[0]
    edges[-1] = centers[-1] + 0.5 * d[-1]
    return edges

def nc_to_tiff_resample_025(
    nc_path: str,
    var_name: str = "atmo",
    time_slice = slice(0, None),     # 例如 slice(31, 41)
    years: tuple | None = None,      # 例如 (1961, 2020)；若 time 可映射到年份
    out_tif: str = "var_mean_025.tif",
    dst_res_deg: float = 0.25,       # 目标角分辨率
    resampling: str = "bilinear",    # nearest / bilinear / cubic
    nodata_val: float = -9999.0,
    compress: str = "LZW",
    out_dtype: str = "float32",
):
    # 1) 读取 & 选时段
    ds = xr.open_dataset(nc_path)
    if var_name not in ds:
        raise ValueError(f"{var_name} 不在文件中。可用变量：{list(ds.data_vars)}")
    da = ds[var_name]  # (time, lat, lon)

    if years is not None:
        if 'year' in da.coords:
            da_sel = da.sel(year=slice(years[0], years[1]))
        else:
            # 根据你的实际起始年修改
            start_year = 1950
            years_arr = (start_year + ds['time'].values).astype(int)
            da = da.assign_coords(year=("time", years_arr))
            da_sel = da.sel(year=slice(years[0], years[1]))
    else:
        da_sel = da.isel(time=time_slice)

    da_mean = da_sel.mean(dim="time", skipna=True)

    # 2) 源网格信息（1°）
    lat = ds['lat'].values
    lon = ds['lon'].values
    lat_edges = edges_from_centers(lat)
    lon_edges = edges_from_centers(lon)
    src_left, src_right = lon_edges.min(), lon_edges.max()
    src_bottom, src_top = lat_edges.min(), lat_edges.max()
    src_h, src_w = len(lat), len(lon)

    src_transform = from_bounds(src_left, src_bottom, src_right, src_top, src_w, src_h)

    # 源数组（注意：rasterio 期望行从北到南，如果 lat 升序，需要翻转）
    src_arr = da_mean.values.astype(out_dtype)
    src_arr = np.flipud(src_arr) if lat[0] < lat[-1] else src_arr
    # NaN -> nodata
    src_mask = ~np.isfinite(src_arr)
    src_arr = src_arr.copy()
    src_arr[src_mask] = nodata_val

    # 3) 目标 0.25° 网格（统一覆盖与源一致的经纬范围）
    dst_w = int(round((src_right - src_left) / dst_res_deg))
    dst_h = int(round((src_top   - src_bottom) / dst_res_deg))
    # 防止四舍五入导致边界外溢，重新用 from_bounds 构造目标 transform
    dst_transform = from_bounds(src_left, src_bottom, src_right, src_top, dst_w, dst_h)

    # 4) 选择重采样方法
    resample_map = {
        "nearest": Resampling.nearest,
        "bilinear": Resampling.bilinear,
        "cubic": Resampling.cubic,
    }
    rs = resample_map.get(resampling.lower(), Resampling.bilinear)

    # 5) 重采样到 0.25°
    dst_arr = np.full((dst_h, dst_w), nodata_val, dtype=src_arr.dtype)
    reproject(
        source=src_arr,
        destination=dst_arr,
        src_transform=src_transform,
        src_crs=crs_wgs84,
        src_nodata=nodata_val,
        dst_transform=dst_transform,
        dst_crs=crs_wgs84,
        dst_nodata=nodata_val,
        resampling=rs
    )

    # 6) 写 GeoTIFF
    profile = {
        "driver": "GTiff",
        "height": dst_h,
        "width":  dst_w,
        "count": 1,
        "dtype": out_dtype,
        "crs": crs_wgs84,
        "transform": dst_transform,
        "nodata": nodata_val,
        "compress": compress
    }
    with rasterio.open(out_tif, "w", **profile) as dst:
        dst.write(dst_arr, 1)

    ds.close()
    return out_tif

# ==== 调用示例 ====
nc_to_tiff_resample_025(
    nc_path="summary_1deg_1104.global.diff.nc",
    var_name="atmo",
    time_slice=slice(31, 41),        # 任选
    out_tif="atmo_mean_025.tif",
    dst_res_deg=0.25,                # 目标 0.25°
    resampling="bilinear",           # 插值平滑；若想严格保持块状可用 "nearest"
    nodata_val=-9999.0,
    out_dtype="float32"
)
