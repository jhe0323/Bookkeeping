import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# ======================
# 用户输入
# ======================
pft_nc = "D:\Work\Research_doc\LULCC\BookKeeping\In_ncfile\IBIS_PFT_dominant_1deg.nc"   # 修改路径
var_name = "maxvegetfrac"    # 你的变量名

# ======================
# 读取数据
# ======================
ds = nc.Dataset(pft_nc)

pft = ds.variables[var_name][0]   # (pft, lat, lon) or (time,pft,lat,lon)
lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]

# 转成 (lat, lon, pft)
if pft.shape[0] == 15:
    pft = np.transpose(pft, (1, 2, 0))
elif pft.shape[-1] == 15:
    pass
else:
    raise ValueError("PFT 维度不识别")

nlat, nlon, npft = pft.shape
print("Shape:", pft.shape)

# ======================
# 1. PFT sum 检查
# ======================
pft_sum = np.sum(pft, axis=2)

print("\n=== PFT SUM CHECK ===")
print("min:", np.nanmin(pft_sum))
print("max:", np.nanmax(pft_sum))
print("mean:", np.nanmean(pft_sum))

# ======================
# 2. dominant PFT
# ======================
dominant = np.argmax(pft, axis=2) + 1

plt.figure(figsize=(8,4))
plt.imshow(dominant, origin="lower")
plt.colorbar(label="Dominant PFT")
plt.title("Dominant PFT")
plt.savefig("dominant_pft.png", dpi=150)

# ======================
# 3. 区域统计
# ======================
def region_mean(lat, lon, data, lat_min, lat_max, lon_min, lon_max):
    mask = (lat[:, None] >= lat_min) & (lat[:, None] <= lat_max) & \
           (lon[None, :] >= lon_min) & (lon[None, :] <= lon_max)
    return data[mask]

regions = {
    "Amazon": (-15, 5, -75, -50),
    "Congo": (-8, 5, 10, 30),
    "India": (8, 30, 68, 90),
    "Europe": (45, 55, 0, 20),
}

print("\n=== REGION PFT MEAN ===")
for name, (la1, la2, lo1, lo2) in regions.items():
    vals = region_mean(lat, lon, pft, la1, la2, lo1, lo2)
    mean_pft = vals.mean(axis=0)
    print(name, np.round(mean_pft, 3))

# ======================
# 4. 简化 biomass density
# ======================
# ⚠️ 这里先用“假设参数”，你可以换成真实 params

# 假设森林高、草地低、crop更低
# ======================
# 4. 修复 biomass density (Orchidee PFT 逻辑)
# ======================
# 为不同 PFT 分配符合生态学逻辑的生物量密度 (tC/ha)
# PFT 1-9: 森林 (高密度)
# PFT 10-13: 草地/农田 (低密度)
# PFT 14-15: 裸地/冰雪 (零密度)
rho = np.array([
    200, 200, 200, 200, 200, 200, 150, 150, 150, # PFT 1-9 (Forests)
    20,  20,  10,  10,                          # PFT 10-13 (Grasses/Crops)
    0,   0                                      # PFT 14-15 (Bare soil/Ice)
])

# 确保 rho 的长度与 npft 一致
if len(rho) < npft:
    rho = np.pad(rho, (0, npft - len(rho)), 'constant')
elif len(rho) > npft:
    rho = rho[:npft]

bio_map = np.tensordot(pft, rho, axes=([2],[0]))

plt.figure(figsize=(8,4))
plt.imshow(bio_map, origin="lower")
plt.colorbar(label="Relative Biomass")
plt.title("Potential Biomass Density")
plt.savefig("biomass_map.png", dpi=150)

# ======================
# 5. 区域 biomass 检查
# ======================
print("\n=== REGION BIOMASS ===")
for name, (la1, la2, lo1, lo2) in regions.items():
    vals = region_mean(lat, lon, bio_map, la1, la2, lo1, lo2)
    print(name, np.nanmean(vals))

# ======================
# 6. 输出 summary
# ======================
print("\n=== QUICK DIAGNOSIS ===")

if np.nanmean(pft_sum) < 0.95 or np.nanmean(pft_sum) > 1.05:
    print("⚠️ PFT sum not normalized")
else:
    print("✓ PFT sum OK")

print("Check if biomass high in Amazon/Congo and lower in India/Europe.")
print("If not → PFT map likely inconsistent.")