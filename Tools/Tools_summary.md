# LULCC Bookkeeping Model Tools

Updated: 2026-05  
Number of tools summarized: 27

This document summarizes the auxiliary tools used in the LULCC bookkeeping model workflow. These scripts support data preprocessing, spatial conversion, diagnostic validation, regional and global statistical aggregation, NetCDF processing, and visualization.

---

## 1. Visualization and Plotting Tools

### 1. `nc_viewer2.py`

`nc_viewer2.py` is the most complete command-line NetCDF visualization tool in the current toolbox.

It supports 2D and 3D NetCDF variables, year-specific extraction, year-to-year or period-difference calculation, and bounding-box-based regional clipping. For carbon-related variables, it can automatically convert values to Gt C units. It also provides flexible colormap options and optional Cartopy-based map projection support.

### 2. `nc_viewer2_per_area.py`

`nc_viewer2_per_area.py` is designed for visualizing area-normalized carbon fluxes, such as g C m⁻² yr⁻¹.

It builds on `nc_viewer2.py` by adding dynamic grid-cell area calculation through `grid_cell_area_m2`. It also includes a fixed colorbar scheme similar to the Global Carbon Budget (GCB) mapping style, making it suitable for comparing spatial carbon flux intensity.

### 3. `nc_slide_viewer_per_area.py`

`nc_slide_viewer_per_area.py` is a non-command-line plotting script derived from the per-area viewer.

It is intended for generating presentation-quality figures for slides and reports. Key parameters are hard-coded at the beginning of the script for easier repeated use. It also prints data distribution information, including extremes and percentiles, which is useful for identifying outliers and setting colorbar limits.

### 4. `Results_reviewer.py`

`Results_reviewer.py` is a lightweight xarray-based tool for quickly reviewing spatial or temporal averages.

It can slice a target time range, calculate the mean field, and visualize the result using dynamic colormap logic. The color scheme automatically adapts to positive and negative values, typically using a green-white-red structure while keeping zero strictly white.

---

## 2. Global and Regional Statistics / Data Extraction

### 5. `NCOutput_to_csv4.py`

`NCOutput_to_csv4.py` is the core tool for converting global NetCDF model outputs into time-series CSV tables.

It supports dynamic variable recognition and does not rely on hard-coded variable lists. It can calculate key carbon budget closure metrics such as `C_sys_total` and `C_bookkeep_total`, generate year-to-year differences for stock variables, remain compatible with both old and new output variable names, and create Gt C-converted columns.

### 6. `LUH2tranistion_view_globe.py`

`LUH2tranistion_view_globe.py` summarizes global LUH2 transition data.

It accumulates transition-area variables from `transitions.nc`, such as variables with `primn_*` patterns, and exports global annual time series to CSV. It supports batch processing by variable prefix or by selected variable names.

### 7. `LUH2tranistion_view_region.py`

`LUH2tranistion_view_region.py`, internally based on `NC_bbox_stats`, extracts spatial statistics for a specific bounding box.

It can calculate the sum, arithmetic mean, or area-weighted mean of selected NetCDF variables within a given region, such as China or the Loess Plateau. It supports time slicing and outputs regional time-series CSV files.

---

## 3. Spatial Data Conversion and Projection Tools

### 8. `NC_to_tif.py`

`NC_to_tif.py` converts NetCDF data to GeoTIFF format.

It first averages the selected NetCDF variable along the time dimension and then uses `rasterio` to export the result as a standard WGS84 GeoTIFF. It can also resample the data to a target resolution, such as 0.25°, for use in GIS software.

### 9. `convert_ibis_vegtype_to_luce_pft_1deg.py`

`convert_ibis_vegtype_to_luce_pft_1deg.py` converts IBIS 0.5° vegetation types to LUCE-compatible 1° one-hot PFT format.

It uses a 2 × 2 majority-vote approach. When there is a tie, cosine-latitude-based area weighting is used to determine the dominant PFT class.

### 10. `convert_ibis_vegtype_to_luce_pft_0.25deg.py`

`convert_ibis_vegtype_to_luce_pft_0.25deg.py` converts IBIS 0.5° vegetation data to LUCE-compatible 0.25° PFT format.

The script uses nearest-neighbor downscaling, effectively splitting each 0.5° grid cell into 2 × 2 0.25° cells. It then converts the dominant vegetation type into a one-hot PFT fraction matrix compatible with the LUCE-style model input.

### 11. `convert_orchidee_pft_to_luce_1deg.py`

`convert_orchidee_pft_to_luce_1deg.py` maps ORCHIDEE 0.25° PFT fractions to LUCE 1° PFT fractions.

It contains ecological mapping rules from the 15 ORCHIDEE PFT classes to the 15 LUCE-compatible PFT classes. Upscaling is performed using cosine-latitude area-weighted averaging.

### 12. `convert_orchidee_pft_to_luce_1deg_diminant.py`

`convert_orchidee_pft_to_luce_1deg_diminant.py` is an extended version of the ORCHIDEE-to-LUCE conversion tool.

In addition to generating weighted-average PFT fractions, it provides a `--dominant` option to extract the single dominant vegetation class and generate a dominant-PFT mask.

### 13. `PFT_making_1deg.py`

`PFT_making_1deg.py` aligns processed PFT fraction maps to the LUH2 `states_1deg.nc` grid.

It strictly copies the latitude and longitude coordinates and array structure from the LUH2 template file, ensuring that downstream model indexing is fully consistent.

### 14. `PFT_making_0.25.py`

`PFT_making_0.25.py` is a specialized script for analyzing and visualizing ESA CCI land-cover-derived PFT layers.

It reads high-resolution PFT fraction layers, handles missing-value masks, reports the spatial distribution of each PFT layer, and plots the dominant PFT map.

---

## 4. Diagnostics, Validation, and Unit Testing

### 15. `check_LUH2_trans_file.py`

`check_LUH2_trans_file.py` is a core diagnostic tool for checking consistency between LUH2 `states.nc` and `transitions.nc`.

It verifies land-area conservation at the grid-cell level and checks whether outgoing transition areas exceed the initial land state. It is especially useful for identifying possible double counting or logical overflow in wood harvest transitions.

### 16. `pft_diagnosis.py`

`pft_diagnosis.py` evaluates the ecological realism of PFT input data.

It checks whether PFT fractions are normalized, generates a potential biomass map using prior biomass-density parameters, and compares key regions such as the Amazon and Europe against ecological expectations to assess whether the PFT distribution is reasonable.

### 17. `ibis_pft_0.25deg_check.py`

`ibis_pft_0.25deg_check.py` performs quality control on LUCE-compatible 0.25° IBIS PFT files.

It checks the dimensions of `maxvegetfrac`, detects negative or out-of-range values, verifies whether land-grid PFT fractions sum to one, and reports mixed-grid statistics.

### 18. `view_pft_ncfile.py`

`view_pft_ncfile.py` is a minimal data-probing script for PFT NetCDF files.

It quickly prints array shape, valid PFT class indices after removing NaN values, and basic minimum/maximum values. It is useful for fast checks during preprocessing.

### 19. `LULCC_test_file_making.py`

`LULCC_test_file_making.py` generates synthetic unit-test data for the bookkeeping model.

It creates a clean 1° artificial dataset and injects user-defined land states and transition events at selected grid cells. This provides a controlled sandbox where model behavior can be manually calculated and compared with code output.

---

## 5. Basic File Operations and NetCDF Processing

### 20. `merge.py`

`merge.py` merges parallel model output files.

It combines multiple `summary_1deg.rank*.nc` files generated by cluster or parallel runs along the longitude dimension (`dim="lon"`) to produce a complete global NetCDF output.

### 21. `diff_nc.py`

`diff_nc.py` calculates NetCDF differences along the time dimension.

It applies xarray `ds.diff()` to generate a `*.diff.nc` file, which directly represents annual changes or time-step differences.

### 22. `nc_cut.py`

`nc_cut.py` extracts and compresses a target time range from a NetCDF file.

It slices the selected year interval and averages along the time axis, producing a static mean field that can be used as a background map or period-mean diagnostic.

---

## 6. Legacy and Archived Scripts

The following scripts are older versions retained for code provenance. They have generally been replaced by newer tools such as `nc_viewer2.py`, `NCOutput_to_csv4.py`, or later model utilities.

### 23. `nc_viewer1.py`

`nc_viewer1.py` is an earlier version of `nc_viewer2.py`.

Some carbon-pool variable names and Gt C conversion rules are hard-coded.

### 24. `nc_viewer.py`

`nc_viewer.py` is the most basic visualization script.

It only provides simple mapping functionality and has been superseded by more flexible plotting tools.

### 25. `NCOutput_to_csv3.py`

`NCOutput_to_csv3.py` is the predecessor of `NCOutput_to_csv4.py`.

It relies on hard-coded lists to extract core variables and diagnostic variables.

### 26. `NCOutput_to_csv2.py`

`NCOutput_to_csv2.py` is an early CSV summary script.

It only extracts a fixed set of six basic carbon-pool variables.

### 27. `NCOutput_to_csv.py`

`NCOutput_to_csv.py` is the earliest CSV conversion tool.

It focuses mainly on atmospheric carbon-pool output, such as `Atmo_total`, rather than full system-level carbon accounting.

---

# LULCC 簿记模型辅助工具总结

更新日期：2026-05  
工具数量：27 个

本文档总结 LULCC 簿记模型中使用的辅助工具脚本。这些脚本主要用于数据预处理、空间格式转换、诊断验证、区域和全球统计汇总、NetCDF 文件处理以及模型结果可视化。

---

## 一、核心可视化工具

### 1. `nc_viewer2.py`

`nc_viewer2.py` 是当前工具箱中功能最完整的命令行 NetCDF 可视化核心工具。

该脚本支持 2D 和 3D 变量读取，支持按年份提取数据，也支持计算年份差值或时间段差值，并可以通过 bounding box 对目标区域进行裁剪。对于碳相关变量，脚本可以自动转换为 Gt C 单位。此外，它还内置了灵活的色带设置和可选的 Cartopy 地图投影支持。

### 2. `nc_viewer2_per_area.py`

`nc_viewer2_per_area.py` 主要用于单位面积碳通量的可视化，例如 g C m⁻² yr⁻¹。

该工具在 `nc_viewer2.py` 的基础上加入了动态网格面积计算函数 `grid_cell_area_m2`，可以将总量变量转换为单位面积通量。同时，它内置了接近 Global Carbon Budget (GCB) 风格的固定色标方案，适合用于空间碳通量强度的比较。

### 3. `nc_slide_viewer_per_area.py`

`nc_slide_viewer_per_area.py` 是基于 per-area 版本修改的免命令行交互绘图脚本。

它主要用于生成 PPT 汇报或论文展示级别的高质量配图。常用参数被硬编码在脚本头部，便于反复使用。同时，它会打印数据分布信息，包括极值和百分位数，方便检查异常值并设置合理的 colorbar 范围。

### 4. `Results_reviewer.py`

`Results_reviewer.py` 是一个轻量级的 xarray 结果快速查看工具。

它可以对目标时间段进行切片并求均值，然后绘制空间分布图。脚本具有动态色带判断逻辑，可以根据数据是否同时包含正值和负值自动选择配色，并严格保证 0 值显示为白色。

---

## 二、全球与区域数据统计汇总

### 5. `NCOutput_to_csv4.py`

`NCOutput_to_csv4.py` 是将簿记模型全球 NetCDF 输出汇总为时间序列 CSV 表格的核心工具。

它支持动态变量识别，不再依赖硬编码变量列表。脚本可以计算核心碳库守恒指标，例如 `C_sys_total` 和 `C_bookkeep_total`，并对库变量计算逐年差分。同时，它兼容新旧输出变量名称，并支持生成 Gt C 单位转换列。

### 6. `LUH2tranistion_view_globe.py`

`LUH2tranistion_view_globe.py` 用于对 LUH2 transitions 数据进行全球统计汇总。

该脚本可以对 `transitions.nc` 中的转移面积变量进行全球累加，例如 `primn_*` 类型变量，并输出年度时间序列 CSV 表格。它支持按变量前缀或变量名称进行批量处理。

### 7. `LUH2tranistion_view_region.py`

`LUH2tranistion_view_region.py` 的内部逻辑基于 `NC_bbox_stats`，用于提取指定区域内的空间统计结果。

在给定区域范围内，例如中国或黄土高原，脚本可以对指定 NetCDF 变量计算总和、简单平均或面积加权平均，并支持按时间切片输出区域时间序列 CSV 文件。

---

## 三、空间数据格式与投影转换

### 8. `NC_to_tif.py`

`NC_to_tif.py` 用于将 NetCDF 文件转换为 GeoTIFF 文件。

该脚本会先沿时间维度对目标变量求平均，然后利用 `rasterio` 将结果输出为标准 WGS84 GeoTIFF 格式。它也可以进行目标分辨率重采样，例如重采样到 0.25°，便于在 GIS 软件中读取和使用。

### 9. `convert_ibis_vegtype_to_luce_pft_1deg.py`

`convert_ibis_vegtype_to_luce_pft_1deg.py` 用于将 IBIS 0.5° 植被类型转换为 LUCE 兼容的 1° one-hot PFT 格式。

该脚本采用 2 × 2 网格多数投票方法。如果出现平局，则使用基于纬度余弦的面积权重来确定主导 PFT 类型。

### 10. `convert_ibis_vegtype_to_luce_pft_0.25deg.py`

`convert_ibis_vegtype_to_luce_pft_0.25deg.py` 用于将 IBIS 0.5° 植被数据转换为 LUCE 兼容的 0.25° PFT 格式。

该脚本采用最邻近降尺度方法，相当于将每个 0.5° 网格拆分为 2 × 2 个 0.25° 网格。随后，它将单一主导植被类型转换为 LUCE 模型兼容的 one-hot PFT 分数矩阵。

### 11. `convert_orchidee_pft_to_luce_1deg.py`

`convert_orchidee_pft_to_luce_1deg.py` 用于将 ORCHIDEE 0.25° PFT 比例图映射为 LUCE 1° PFT 比例图。

脚本内置 ORCHIDEE 15 类 PFT 到 LUCE 15 类 PFT 的生态学映射规则，并通过纬度余弦面积加权平均进行升尺度。

### 12. `convert_orchidee_pft_to_luce_1deg_diminant.py`

`convert_orchidee_pft_to_luce_1deg_diminant.py` 是 ORCHIDEE 到 LUCE 转换工具的进阶版本。

除了生成面积加权平均的 PFT 比例图外，它还提供 `--dominant` 开关，用于提取最高比例的单一植被类型，并生成主导 PFT 掩膜数据。

### 13. `PFT_making_1deg.py`

`PFT_making_1deg.py` 用于将处理后的 PFT 比例图严格对齐到 LUH2 的 `states_1deg.nc` 模板网格。

该脚本会精确复制模板文件中的纬度、经度坐标和数组结构，确保后续模型运行中的数组索引完全一致。

### 14. `PFT_making_0.25.py`

`PFT_making_0.25.py` 是针对 ESA CCI 土地覆盖源数据的专用分析与可视化脚本。

它可以读取高分辨率 PFT 分数层，处理缺失值掩膜，输出各层 PFT 的空间分布状态，并绘制主导 PFT 空间图。

---

## 四、模型诊断、检验与单元测试

### 15. `check_LUH2_trans_file.py`

`check_LUH2_trans_file.py` 是检查 LUH2 `states.nc` 和 `transitions.nc` 数据一致性的核心诊断工具。

该脚本可以逐网格验证土地面积守恒，并检查转移转出面积是否超过初始状态面积。它尤其适合用于排查木材采伐过程中的重复计算和逻辑越界问题。

### 16. `pft_diagnosis.py`

`pft_diagnosis.py` 用于检查 PFT 输入数据的生态学合理性。

该脚本会检查 PFT 分数是否归一化，并利用预设的生物量密度参数生成潜在生物量地图。通过结合重点区域常识，例如亚马逊和欧洲，可以反向判断 PFT 分布是否合理。

### 17. `ibis_pft_0.25deg_check.py`

`ibis_pft_0.25deg_check.py` 是 LUCE 兼容 IBIS 0.25° PFT 文件的质量控制脚本。

它会检查 `maxvegetfrac` 的数据维度、负值和超范围值，验证陆地网格上的 PFT 分数是否求和为 1，并统计混合网格情况。

### 18. `view_pft_ncfile.py`

`view_pft_ncfile.py` 是一个极简 PFT NetCDF 数据探针脚本。

它可以快速打印 PFT 文件的数组形状、去除 NaN 后的有效类别索引以及基本极值，适合在预处理流程中进行快速抽查。

### 19. `LULCC_test_file_making.py`

`LULCC_test_file_making.py` 用于生成簿记模型的合成单元测试数据。

它会构建一个干净的 1° 虚拟数据集，并在指定经纬度网格中注入用户定义的土地状态和转移事件。这样可以提供一个可手算、无干扰的测试环境，用于验证主程序核心逻辑是否正确。

---

## 五、基础文件合并与处理

### 20. `merge.py`

`merge.py` 是并行输出文件合并工具。

它可以将集群或并行计算生成的多个 `summary_1deg.rank*.nc` 分块文件沿经度维度 `dim="lon"` 拼接为完整的全球 NetCDF 输出文件。

### 21. `diff_nc.py`

`diff_nc.py` 用于沿时间维度计算 NetCDF 差分。

它调用 xarray 的 `ds.diff()` 生成 `*.diff.nc` 文件，结果可以直接表示年度变化量或相邻时间步差值。

### 22. `nc_cut.py`

`nc_cut.py` 用于对 NetCDF 时间序列进行切片和均值压缩。

该脚本可以从原始数据中裁剪目标年份区间，并沿时间轴求平均，从而生成静态平均场，可用于背景图、阶段平均或诊断分析。

---

## 六、早期旧版与存档文件

以下脚本为早期版本，主要用于代码溯源。目前大多已被 `nc_viewer2.py`、`NCOutput_to_csv4.py` 或后续版本工具替代。

### 23. `nc_viewer1.py`

`nc_viewer1.py` 是 `nc_viewer2.py` 的前身。

其中部分碳库变量名称和 Gt C 单位转换逻辑仍采用硬编码方式。

### 24. `nc_viewer.py`

`nc_viewer.py` 是最基础版本的可视化脚本。

它只提供简单的空间绘图功能，已被后续更灵活的绘图工具替代。

### 25. `NCOutput_to_csv3.py`

`NCOutput_to_csv3.py` 是 `NCOutput_to_csv4.py` 的前身。

它依赖硬编码变量列表来提取核心变量和诊断变量。

### 26. `NCOutput_to_csv2.py`

`NCOutput_to_csv2.py` 是更早期的 CSV 汇总脚本。

它只提取固定的 6 个基础碳库变量。

### 27. `NCOutput_to_csv.py`

`NCOutput_to_csv.py` 是最早期的 CSV 转换工具。

它主要关注大气碳库输出，例如 `Atmo_total`，而不是完整的系统碳收支核算。
