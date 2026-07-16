import json
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask


def guess_coord_name(ds, candidates):
    for name in candidates:
        if name in ds.coords or name in ds.variables:
            return name
    raise ValueError(f"Cannot find coordinate name from candidates: {candidates}")


def make_region_mask_from_shp(
    template_nc,
    shp_path,
    region_col,
    out_mask_nc,
    out_region_csv=None,
):
    """
    Convert a regional shapefile to a NetCDF region mask using the lat/lon grid
    from a template model output nc file.

    Parameters
    ----------
    template_nc : str
        Path to model output nc file, e.g. summary_025deg.global.nc.
        This file only provides lat/lon grid information.
    shp_path : str
        Path to shapefile.
    region_col : str
        Column name in the shapefile that stores region names or region IDs.
    out_mask_nc : str
        Output nc file path for the region mask.
    out_region_csv : str or None
        Output csv file for region ID-name mapping.
    """

    # ============================================================
    # 1. Read template nc grid
    # ============================================================
    ds = xr.open_dataset(template_nc)

    lat_name = guess_coord_name(ds, ["lat", "latitude"])
    lon_name = guess_coord_name(ds, ["lon", "longitude"])

    lat = ds[lat_name]
    lon = ds[lon_name]

    if lat.ndim != 1 or lon.ndim != 1:
        raise ValueError("This script currently supports only 1D lat/lon coordinates.")

    lat_values = lat.values
    lon_values = lon.values

    print("Template NC grid:")
    print(f"  lat name: {lat_name}, size: {lat_values.size}, range: {lat_values.min()} to {lat_values.max()}")
    print(f"  lon name: {lon_name}, size: {lon_values.size}, range: {lon_values.min()} to {lon_values.max()}")

    # ============================================================
    # 2. Read shapefile
    # ============================================================
    gdf = gpd.read_file(shp_path)

    print("\nShapefile columns:")
    print(list(gdf.columns))

    # 清理区域名称中的空格、换行符
    gdf[region_col] = (
        gdf[region_col]
        .astype(str)
        .str.replace("\r", "", regex=False)
        .str.replace("\n", "", regex=False)
        .str.strip()
    )

    print("\nCONTINENT values after cleaning:")
    print(gdf[region_col].unique())

    if region_col not in gdf.columns:
        raise ValueError(
            f"Cannot find region_col='{region_col}' in shapefile. "
            f"Available columns are: {list(gdf.columns)}"
        )

    if gdf.crs is None:
        raise ValueError(
            "The shapefile has no CRS. Please define its CRS first, usually EPSG:4326."
        )

    # Convert shapefile to WGS84 lon/lat
    gdf = gdf.to_crs("EPSG:4326")

    # Keep only region column and geometry
    gdf = gdf[[region_col, "geometry"]].copy()
    gdf = gdf.dropna(subset=[region_col, "geometry"])
    gdf[region_col] = gdf[region_col].astype(str)

    # Dissolve polygons by region name.
    # This is important if one region contains multiple polygons.
    gdf = gdf.dissolve(by=region_col, as_index=False)

    # Remove empty geometry
    gdf = gdf[~gdf.geometry.is_empty].copy()

    # Try to fix invalid geometries if present
    if not gdf.is_valid.all():
        print("\nWarning: invalid geometries found. Trying buffer(0) fix...")
        gdf["geometry"] = gdf.geometry.buffer(0)

    # Sort regions to make region ID stable
    gdf = gdf.sort_values(region_col).reset_index(drop=True)

    region_names = gdf[region_col].tolist()
    region_ids = list(range(len(region_names)))

    print("\nRegions:")
    for rid, rname in zip(region_ids, region_names):
        print(f"  {rid}: {rname}")

    # ============================================================
    # 3. Handle longitude convention
    # ============================================================
    # Shapefile is normally -180 to 180.
    # If nc longitude is 0 to 360, convert lon values only for masking.
    # The output mask still keeps the original nc lon coordinates.
    if np.nanmax(lon_values) > 180:
        lon_for_mask = ((lon_values + 180) % 360) - 180
        print("\nDetected 0-360 longitude in NC. Using converted -180 to 180 longitude for masking.")
    else:
        lon_for_mask = lon_values.copy()

    # ============================================================
    # 4. Create region mask
    # ============================================================
    regions = regionmask.Regions(
        numbers=region_ids,
        names=region_names,
        abbrevs=region_names,
        outlines=gdf.geometry.tolist(),
    )

    mask = regions.mask(lon_for_mask, lat_values)

    # regionmask usually returns dims named lat/lon
    # Rename and assign original nc coordinates
    rename_dict = {}
    if "lat" in mask.dims and lat_name != "lat":
        rename_dict["lat"] = lat_name
    if "lon" in mask.dims and lon_name != "lon":
        rename_dict["lon"] = lon_name

    if rename_dict:
        mask = mask.rename(rename_dict)

    mask = mask.assign_coords({
        lat_name: lat_values,
        lon_name: lon_values,
    })

    # Convert NaN outside polygons to -9999
    mask_int = mask.fillna(-9999).astype("int16")
    mask_int.name = "region_mask"

    region_mapping = {int(rid): str(rname) for rid, rname in zip(region_ids, region_names)}

    mask_int.attrs["long_name"] = "region mask from shapefile"
    mask_int.attrs["description"] = "Grid cells are assigned by cell center location."
    mask_int.attrs["missing_value"] = -9999
    mask_int.attrs["outside_region_value"] = -9999
    mask_int.attrs["region_col"] = region_col
    mask_int.attrs["region_mapping_json"] = json.dumps(region_mapping, ensure_ascii=False)

    # ============================================================
    # 5. Save mask nc
    # ============================================================
    out_ds = xr.Dataset(
        data_vars={
            "region_mask": mask_int,
            "region_name": xr.DataArray(
                np.array(region_names, dtype=str),
                dims=("region",),
                coords={"region": region_ids},
            ),
        },
        coords={
            lat_name: lat_values,
            lon_name: lon_values,
            "region": region_ids,
        },
    )

    # Copy simple coordinate attributes if available
    if hasattr(lat, "attrs"):
        out_ds[lat_name].attrs.update(lat.attrs)
    if hasattr(lon, "attrs"):
        out_ds[lon_name].attrs.update(lon.attrs)

    encoding = {
        "region_mask": {
            "dtype": "int16",
            "_FillValue": -9999,
            "zlib": True,
            "complevel": 4,
        }
    }

    out_ds.to_netcdf(out_mask_nc, encoding=encoding)

    print(f"\nRegion mask saved to:")
    print(f"  {out_mask_nc}")

    # ============================================================
    # 6. Save region ID table
    # ============================================================
    if out_region_csv is None:
        out_region_csv = out_mask_nc.replace(".nc", "_region_table.csv")

    region_table = pd.DataFrame({
        "region_id": region_ids,
        "region_name": region_names,
    })

    region_table.to_csv(out_region_csv, index=False, encoding="utf-8-sig")

    print(f"\nRegion table saved to:")
    print(f"  {out_region_csv}")

    # ============================================================
    # 7. Basic check
    # ============================================================
    print("\nGrid cell count by region:")
    values, counts = np.unique(mask_int.values, return_counts=True)

    for val, count in zip(values, counts):
        if val == -9999:
            print(f"  outside / ocean / no region: {count}")
        else:
            print(f"  {val} - {region_names[int(val)]}: {count}")

    ds.close()

    return out_ds


if __name__ == "__main__":
    make_region_mask_from_shp(
        template_nc="summary_025deg.global_time_gt1100_cut_1160-1170_5.5.nc",
        shp_path="World_Continents_guo.shp",
        region_col="CONTINENT",
        out_mask_nc="global_6_regions_mask_025deg.nc",
        out_region_csv="global_6_regions_mask_025deg_region_table.csv",
    )