#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import regionmask
import pandas as pd
import numpy as np


def country_timeseries(
    nc_path,
    out_csv,
    country_name="Indonesia",
    vars_to_sum=None,
):
    """
    Extract country-level annual time series from a gridded NetCDF file
    using Natural Earth country boundaries via regionmask.

    Parameters
    ----------
    nc_path : str
        Input NetCDF path.
    out_csv : str
        Output CSV path.
    country_name : str
        Country name, e.g., "China", "United States of America", "Brazil".
    vars_to_sum : list[str]
        Variables to aggregate. If None, a default list is used.
    """

    if vars_to_sum is None:
        vars_to_sum = [
            "Gross_Sources",
            "Gross_Sinks",
            "Net_Emissions",
            "Flux_Harvest_Net",
            "Flux_WHp",
            "Flux_Harvest_SoilSlash",
            "Flux_Harvest_Regrowth",
            "Flux_Clearing",
            "Flux_Abandonment",
            "Flux_Other",
            "Flux_Products_Total",
            "area_harvest",
            "harvest_luh2_area",
            "harvest_effective_area",
        ]

    ds = xr.open_dataset(nc_path)

    # ---- identify coordinate names ----
    lat_name = "lat" if "lat" in ds.coords else "latitude"
    lon_name = "lon" if "lon" in ds.coords else "longitude"

    lat = ds[lat_name]
    lon = ds[lon_name]

    # ---- regionmask uses lon usually in -180 to 180 ----
    # If your lon is 0–360, convert it to -180–180 and sort.
    if float(lon.max()) > 180:
        ds = ds.assign_coords(
            {lon_name: (((ds[lon_name] + 180) % 360) - 180)}
        ).sortby(lon_name)
        lon = ds[lon_name]
        lat = ds[lat_name]

    # ---- Natural Earth country boundaries ----
    countries = regionmask.defined_regions.natural_earth_v5_0_0.countries_110

    # Find country index
    names = list(countries.names)
    matches = [i for i, name in enumerate(names) if country_name.lower() in name.lower()]

    if len(matches) == 0:
        print("Available example country names:")
        print(names[:30])
        raise ValueError(f"Country not found: {country_name}")

    if len(matches) > 1:
        print("Multiple matches found:")
        for i in matches:
            print(i, names[i])
        raise ValueError("Please use a more specific country name.")

    country_idx = matches[0]
    print(f"Selected country: {countries.names[country_idx]}")

    # ---- make 2D country mask ----
    mask = countries.mask(lon, lat)

    # country_mask: True inside selected country
    country_mask = mask == country_idx

    # ---- aggregate variables ----
    out = pd.DataFrame()

    # time coordinate
    if "time" in ds.coords:
        out["time"] = ds["time"].values
    else:
        out["time"] = np.arange(ds.dims["time"])

    for var in vars_to_sum:
        if var not in ds:
            print(f"[Skip] variable not found: {var}")
            continue

        da = ds[var]

        # Only process variables with lat/lon dimensions
        if lat_name not in da.dims or lon_name not in da.dims:
            print(f"[Skip] no spatial dims: {var}")
            continue

        # Sum all grid cells inside the country
        country_ts = da.where(country_mask).sum(dim=[lat_name, lon_name], skipna=True)

        out[var + "_country_total"] = country_ts.values

    out.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    country_timeseries(
        nc_path=r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\summary_025deg.global_time_ge1000_7.9.nc",
        out_csv="Dem.Rep.Congo_ELUC_timeseries.csv",
        country_name="Dem. Rep. Congo",
    )