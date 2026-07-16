import argparse
import os
import numpy as np
import pandas as pd
import netCDF4 as nc

SELECTED_VARS = [
    "Gross_Sources",
    "Gross_Sinks",
    "Net_Emissions",
    "Closure_Error",

    "Flux_FD",
    "Flux_NFC",
    "Flux_FR",
    "Flux_NFR",
    "Flux_CAL",
    "Flux_WHp",

    "Flux_Clearing",
    "Flux_Abandonment",

    "harvest_requested_biomass",
    "harvest_met_biomass",
    "harvest_unmet_biomass",
    "harvest_unmet_raw_biomass",
    "harvest_forced_biomass",
    "harvest_luh2_area_frac",
    "harvest_effective_area_frac",
]

def guess_name(ds, candidates):
    for name in candidates:
        if name in ds.variables or name in ds.dimensions:
            return name
    raise ValueError(f"Cannot find name from candidates: {candidates}")


def read_region_names(mask_ds, region_ids):
    """
    Try to read region names from mask nc.
    If not found, use Region_<id>.
    """
    if "region_name" not in mask_ds.variables:
        return {rid: f"Region_{rid}" for rid in region_ids}

    raw = mask_ds.variables["region_name"][:]

    try:
        names = [str(x) for x in raw]
    except Exception:
        names = [f"Region_{rid}" for rid in region_ids]

    mapping = {}
    for rid in region_ids:
        if rid < len(names):
            mapping[rid] = names[rid].replace("\r", "").replace("\n", "").strip()
        else:
            mapping[rid] = f"Region_{rid}"

    return mapping


def get_time_values(ds, time_name):
    time_values = ds.variables[time_name][:]

    years = []
    for x in time_values:
        try:
            years.append(int(x))
        except Exception:
            years.append(x)

    return years


def get_vars_to_process(
    ds,
    time_name,
    lat_name,
    lon_name,
    selected_vars=None,
):
    """
    Select numeric variables containing time, lat, lon dimensions.
    """
    if selected_vars is not None and len(selected_vars) > 0:
        missing = [v for v in selected_vars if v not in ds.variables]
        if missing:
            print("\nWarning: selected variables not found and will be skipped:")
            for v in missing:
                print("  ", v)

        candidates = [v for v in selected_vars if v in ds.variables]
    else:
        candidates = list(ds.variables.keys())

    vars_to_process = []

    for v in candidates:
        var = ds.variables[v]

        if v in [time_name, lat_name, lon_name]:
            continue

        dims = var.dimensions

        if time_name not in dims:
            continue
        if lat_name not in dims:
            continue
        if lon_name not in dims:
            continue

        if not np.issubdtype(var.dtype, np.number):
            continue

        vars_to_process.append(v)

    return vars_to_process


def read_var_time_as_latlon(
    var,
    t,
    time_name,
    lat_name,
    lon_name,
):
    """
    Read one time slice from a variable and return a 2D lat-lon array.

    Supports:
      time, lat, lon
      time, lon, lat
      time, pft, lat, lon
      time, lat, lon, pft

    Any extra dimensions will be summed.
    """

    dims = var.dimensions

    index = []
    for d in dims:
        if d == time_name:
            index.append(t)
        else:
            index.append(slice(None))

    arr = var[tuple(index)]
    arr = np.ma.filled(arr, np.nan).astype("float64", copy=False)

    remaining_dims = [d for d in dims if d != time_name]

    lat_axis = remaining_dims.index(lat_name)
    lon_axis = remaining_dims.index(lon_name)

    # Move lat and lon to first two axes
    arr = np.moveaxis(arr, [lat_axis, lon_axis], [0, 1])

    # Sum all extra dimensions, e.g. pft/category
    if arr.ndim > 2:
        extra_axes = tuple(range(2, arr.ndim))
        arr = np.nansum(arr, axis=extra_axes)

    return arr


def summarize_nc_by_region_mask(
    nc_path,
    mask_path,
    out_prefix,
    selected_vars=None,
):
    """
    Summarize model nc outputs by region mask.

    Parameters
    ----------
    nc_path : str
        Model output nc file.
    mask_path : str
        Region mask nc file generated from shp.
    out_prefix : str
        Output file prefix.
    selected_vars : list or None
        If None, all numeric time-lat-lon variables will be summarized.
    """

    print("Opening files...")
    ds = nc.Dataset(nc_path)
    mask_ds = nc.Dataset(mask_path)

    try:
        time_name = guess_name(ds, ["time"])
        lat_name = guess_name(ds, ["lat", "latitude"])
        lon_name = guess_name(ds, ["lon", "longitude"])

        if "region_mask" not in mask_ds.variables:
            raise ValueError("Cannot find variable 'region_mask' in mask nc.")

        mask = mask_ds.variables["region_mask"][:]
        mask = np.ma.filled(mask, -9999).astype("int32")

        lat = ds.variables[lat_name][:]
        lon = ds.variables[lon_name][:]

        mask_lat_name = guess_name(mask_ds, ["lat", "latitude"])
        mask_lon_name = guess_name(mask_ds, ["lon", "longitude"])

        mask_lat = mask_ds.variables[mask_lat_name][:]
        mask_lon = mask_ds.variables[mask_lon_name][:]

        print("\nModel NC grid:")
        print(f"  lat: {lat_name}, size={len(lat)}, range={float(np.min(lat))} to {float(np.max(lat))}")
        print(f"  lon: {lon_name}, size={len(lon)}, range={float(np.min(lon))} to {float(np.max(lon))}")

        print("\nMask NC grid:")
        print(f"  lat: {mask_lat_name}, size={len(mask_lat)}, range={float(np.min(mask_lat))} to {float(np.max(mask_lat))}")
        print(f"  lon: {mask_lon_name}, size={len(mask_lon)}, range={float(np.min(mask_lon))} to {float(np.max(mask_lon))}")

        if mask.shape != (len(lat), len(lon)):
            raise ValueError(
                f"Mask shape {mask.shape} does not match nc grid {(len(lat), len(lon))}."
            )

        if not np.allclose(lat, mask_lat):
            raise ValueError("Latitude values in model nc and mask nc do not match.")

        if not np.allclose(lon, mask_lon):
            raise ValueError("Longitude values in model nc and mask nc do not match.")

        region_ids = sorted([int(x) for x in np.unique(mask) if int(x) >= 0])
        region_names = read_region_names(mask_ds, region_ids)

        print("\nRegions detected:")
        for rid in region_ids:
            count = int(np.sum(mask == rid))
            print(f"  {rid}: {region_names[rid]}, grid cells={count}")

        years = get_time_values(ds, time_name)

        vars_to_process = get_vars_to_process(
            ds,
            time_name=time_name,
            lat_name=lat_name,
            lon_name=lon_name,
            selected_vars=selected_vars,
        )

        if not vars_to_process:
            raise ValueError("No valid time-lat-lon variables found.")

        print("\nVariables to summarize:")
        for v in vars_to_process:
            print(f"  {v}: dims={ds.variables[v].dimensions}, shape={ds.variables[v].shape}")

        records = []

        print("\nStart summarizing...")

        for v in vars_to_process:
            var = ds.variables[v]
            print(f"\nProcessing variable: {v}")

            for t_idx, year in enumerate(years):
                arr = read_var_time_as_latlon(
                    var=var,
                    t=t_idx,
                    time_name=time_name,
                    lat_name=lat_name,
                    lon_name=lon_name,
                )

                for rid in region_ids:
                    region_bool = mask == rid
                    value = float(np.nansum(arr[region_bool]))

                    records.append({
                        "Year": year,
                        "region_id": rid,
                        "region": region_names[rid],
                        "variable": v,
                        "value": value,
                    })

                if (t_idx + 1) % 50 == 0 or (t_idx + 1) == len(years):
                    print(f"  finished {t_idx + 1}/{len(years)} years")

        long_df = pd.DataFrame(records)

        wide_df = long_df.pivot_table(
            index=["Year", "region_id", "region"],
            columns="variable",
            values="value",
        ).reset_index()

        wide_df.columns.name = None

        long_csv = f"{out_prefix}_long.csv"
        wide_csv = f"{out_prefix}_wide.csv"
        xlsx_path = f"{out_prefix}.xlsx"

        long_df.to_csv(long_csv, index=False, encoding="utf-8-sig")
        wide_df.to_csv(wide_csv, index=False, encoding="utf-8-sig")

        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            long_df.to_excel(writer, sheet_name="long", index=False)
            wide_df.to_excel(writer, sheet_name="wide", index=False)

        print("\nDone.")
        print(f"Long table: {long_csv}")
        print(f"Wide table: {wide_csv}")
        print(f"Excel file: {xlsx_path}")

        return long_df, wide_df

    finally:
        ds.close()
        mask_ds.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--nc",
        required=True,
        help="Model output nc file, e.g. summary_025deg.global.nc",
    )

    parser.add_argument(
        "--mask",
        required=True,
        help="Region mask nc file, e.g. global_6_regions_mask_025deg.nc",
    )

    parser.add_argument(
        "--out_prefix",
        default="region_summary",
        help="Output prefix",
    )

    parser.add_argument(
        "--vars",
        nargs="*",
        default=None,
        help="Optional variable names to summarize. If omitted, use default selected variables.",
    )

    args = parser.parse_args()

    if args.vars is None:
        selected_vars = SELECTED_VARS
    else:
        selected_vars = args.vars

    summarize_nc_by_region_mask(
        nc_path=args.nc,
        mask_path=args.mask,
        out_prefix=args.out_prefix,
        selected_vars=selected_vars,
    )