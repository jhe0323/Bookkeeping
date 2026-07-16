#!/usr/bin/env python3
"""
plot_cci_pft_halfdeg.py

读取 ESA_CCI_LC_PFT_halfdeg_2000.nc 中的 PFTfrac(PFT,lat,lon)，掩掉 missing_value，
打印每层统计，并绘制指定 PFT、全部 PFT 小图和主导 PFT 地图。

依赖:
  pip install xarray netcdf4 numpy matplotlib
可选（更好的地图）:
  pip install cartopy
"""

import os
import sys
import math
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except Exception:
    HAVE_CARTOPY = False

# ============= 用户参数（按需修改） =============
nc_path = r"D:/Work/ORCHIDEE/Work/LULCC/ESA_CCI_LC_PFT_halfdeg_2000.nc"
output_prefix = "ESA_CCI_PFT_2000"
pft_value_to_plot = 9    # 想单独查看哪个 PFT（按 PFT 变量的实际值，例如 1..15）
clip_high_percentile = 99.5
clip_low_percentile = 1.0
# =================================================

def ensure_lon_range(lon):
    lon = np.array(lon)
    # 将 lon 转为 -180..180 并返回排序索引
    if lon.min() >= 0 and lon.max() > 180:
        lon2 = (lon + 180) % 360 - 180
    else:
        lon2 = lon.copy()
    order = np.argsort(lon2)
    return lon2[order], order

def plot_map(lon, lat, arr2d, title, outfn, vmin=None, vmax=None, cmap='viridis', discrete=False, pft_labels=None):
    fig = plt.figure(figsize=(12,6))
    if HAVE_CARTOPY:
        ax = plt.axes(projection=ccrs.Robinson())
        ax.set_global()
        pcm = ax.pcolormesh(lon, lat, arr2d, transform=ccrs.PlateCarree(), shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.add_feature(cfeature.COASTLINE.with_scale('110m'), linewidth=0.4)
        ax.add_feature(cfeature.BORDERS.with_scale('110m'), linewidth=0.2)
    else:
        ax = fig.add_subplot(1,1,1)
        extent = [lon.min(), lon.max(), lat.min(), lat.max()]
        pcm = ax.imshow(arr2d, origin='lower', extent=extent, aspect='auto', vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
    if discrete and pft_labels is not None:
        cb = plt.colorbar(pcm, ax=ax, orientation='horizontal', pad=0.05, ticks=np.arange(len(pft_labels)) + 0.5)
        cb.set_ticklabels([str(int(x)) for x in pft_labels])
        cb.set_label('PFT value')
    else:
        cbar = plt.colorbar(pcm, ax=ax, orientation='horizontal', pad=0.05)
        cbar.set_label('PFT fraction')
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outfn, dpi=300)
    print("Saved:", outfn)
    plt.show()

def main():
    if not os.path.exists(nc_path):
        print("NetCDF 文件不存在：", nc_path)
        sys.exit(1)

    ds = xr.open_dataset(nc_path)
    print(ds)  # quick print

    # 读取坐标与 PFT 列表
    lat = ds['lat'].values
    lon = ds['lon'].values
    PFT_vals = ds['PFT'].values  # e.g. [1,2,...,15]

    lon2, lon_order = ensure_lon_range(lon)

    # 读取 PFTfrac，并掩掉 missing_value
    da = ds['PFTfrac']  # dims (PFT, lat, lon)
    mv = None
    if 'missing_value' in da.attrs:
        mv = float(da.attrs['missing_value'])
        print("Detected missing_value:", mv)
    else:
        print("No missing_value attribute found; will also mask extremely large values >1e30.")
        mv = None

    arr = da.values.astype('float64')  # shape (PFT, lat, lon)
    # mask missing_value and extremely large values
    if mv is not None:
        arr[arr == mv] = np.nan
    arr[np.abs(arr) > 1e30] = np.nan

    n_pft = arr.shape[0]
    nlat = arr.shape[1]; nlon = arr.shape[2]
    total_cells = nlat * nlon

    # 每层统计
    print("\nPer-PFT stats (after masking): PFT, finite_cells, min, max, p50, p99")
    for i in range(n_pft):
        layer = arr[i,:,:]
        finite_mask = np.isfinite(layer)
        finite_count = finite_mask.sum()
        if finite_count == 0:
            print(f"PFT={int(PFT_vals[i])}: all NaN")
            continue
        mn = np.nanmin(layer); mx = np.nanmax(layer)
        p50 = np.nanpercentile(layer,50)
        p99 = np.nanpercentile(layer,99)
        print(f"PFT={int(PFT_vals[i])}: finite={finite_count} ({finite_count/total_cells*100:.2f}%), min={mn:.3e}, max={mx:.3e}, p50={p50:.3e}, p99={p99:.3e}")

    # 计算主导 PFT（argmax），先把 NaN 视为 -inf
    arr_nonan = np.where(np.isfinite(arr), arr, -np.inf)
    dom_idx = np.argmax(arr_nonan, axis=0)  # 0-based indices of PFT array
    # 标记所有层都为 NaN 的格点
    all_nan_mask = ~np.any(np.isfinite(arr), axis=0)
    valid_cells = np.sum(~all_nan_mask)
    print(f"\nValid gridcells with >=1 finite PFT fraction: {valid_cells} / {total_cells} ({valid_cells/total_cells*100:.2f}%)")

    # 统计每个 PFT 被选为主导的格点数量
    print("\nDominant counts per PFT (value, count, pct_of_valid):")
    for i in range(n_pft):
        cnt = int(np.sum((dom_idx == i) & (~all_nan_mask)))
        pct = cnt / valid_cells * 100 if valid_cells>0 else 0
        print(f"PFT={int(PFT_vals[i])}: count={cnt}, pct_valid={pct:.3f}%")

    # ---- 绘制单个指定 PFT ----
    # 找到 pft_value_to_plot 在 PFT_vals 的索引
    try:
        pv = float(pft_value_to_plot)
        idxs = np.where(PFT_vals == pv)[0]
        if idxs.size == 0:
            raise ValueError(f"PFT value {pft_value_to_plot} not found in PFT variable")
        idx = int(idxs[0])
    except Exception as e:
        print("Error finding requested PFT:", e)
        idx = 0
        print("Will plot first PFT by default.")

    layer = arr[idx,:,:][:, lon_order]  # reorder lon to lon2
    finite_vals = layer[np.isfinite(layer)]
    if finite_vals.size > 0:
        vmin = np.percentile(finite_vals, clip_low_percentile)
        vmax = np.percentile(finite_vals, clip_high_percentile)
    else:
        vmin = None; vmax = None

    plot_map(lon2, lat, layer, f"PFT fraction for PFT={int(PFT_vals[idx])}", f"{output_prefix}_PFT{int(PFT_vals[idx])}.png", vmin=vmin, vmax=vmax)

    # ---- 绘制全部 PFT 的小图面板 ----
    nonzero_idxs = [i for i in range(n_pft) if np.isfinite(arr[i,:,:]).sum() > 0]
    if len(nonzero_idxs) > 0:
        cols = 5
        rows = math.ceil(len(nonzero_idxs) / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3*rows), squeeze=False)
        for k, i in enumerate(nonzero_idxs):
            r = k // cols; c = k % cols
            ax = axes[r][c]
            dat = arr[i,:,:][:, lon_order]
            im = ax.imshow(dat, origin='lower', extent=[lon2.min(), lon2.max(), lat.min(), lat.max()], aspect='auto')
            ax.set_title(f"PFT {int(PFT_vals[i])}")
            ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        # hide unused axes
        for j in range(len(nonzero_idxs), rows*cols):
            r = j // cols; c = j % cols
            axes[r][c].axis('off')
        plt.tight_layout()
        fn_all = f"{output_prefix}_allPFT_smallpanels.png"
        fig.savefig(fn_all, dpi=200)
        print("Saved:", fn_all)
        plt.show()

    # ---- 绘制主导 PFT 地图（显示 PFT value） ----
    dom_vals = np.full_like(dom_idx, np.nan, dtype='float32')
    for i in range(n_pft):
        dom_vals[dom_idx==i] = PFT_vals[i]
    dom_vals[all_nan_mask] = np.nan
    dom_vals = dom_vals[:, lon_order]

    # discrete cmap
    cmap = plt.cm.get_cmap('tab20', n_pft)
    fig = plt.figure(figsize=(12,6))
    if HAVE_CARTOPY:
        ax = plt.axes(projection=ccrs.Robinson())
        ax.set_global()
        pcm = ax.pcolormesh(lon2, lat, dom_vals, transform=ccrs.PlateCarree(), shading='auto', cmap=cmap)
        ax.add_feature(cfeature.COASTLINE.with_scale('110m'), linewidth=0.4)
    else:
        ax = fig.add_subplot(1,1,1)
        pcm = ax.imshow(dom_vals, origin='lower', extent=[lon2.min(), lon2.max(), lat.min(), lat.max()], aspect='auto', cmap=cmap)
        ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')

    cb = plt.colorbar(pcm, ax=ax, orientation='horizontal', pad=0.05, ticks=np.arange(n_pft)+0.5)
    cb.set_ticklabels([str(int(x)) for x in PFT_vals])
    cb.set_label('Dominant PFT value')
    ax.set_title("Dominant PFT (argmax of PFTfrac)")
    plt.tight_layout()
    fn_dom = f"{output_prefix}_dominantPFT.png"
    fig.savefig(fn_dom, dpi=300)
    print("Saved:", fn_dom)
    plt.show()

    ds.close()

if __name__ == "__main__":
    main()
