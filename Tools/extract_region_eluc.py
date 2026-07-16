#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import regionmask
import pandas as pd
import numpy as np


def extract_ar6_region_groups(
    nc_path,
    out_csv,
    region_groups=None,
    vars_to_sum=None,
):
    """
    Aggregate gridded NetCDF variables to IPCC AR6 regional groups.

    This method does not require user-provided shapefiles.
    It uses regionmask built-in AR6 land regions.

    Notes
    -----
    This assumes variables are grid-cell total fluxes.
    If variables are per-area fluxes, multiply by grid-cell area before summing.
    """

    if region_groups is None:
        region_groups = {
            "Europe_AR6": ["NEU", "WCE", "EEU", "MED"],
            "Southeast_Asia_AR6": ["SEA"],
            "North_America_AR6": ["NWN", "NEN", "WNA", "CNA", "ENA"],
            "South_America_AR6": ["NWS", "NSA", "NES", "SAM", "SWS", "SES", "SSA"],
            "Africa_AR6": ["SAH", "WAF", "CAF", "NEAF", "SEAF", "WSAF", "ESAF", "MDG"],
        }

    if vars_to_sum is None:
        vars_to_sum = [
            "Gross_Sources",
            "Gross_Sinks",
            "Net_Emissions",
        ]

    ds = xr.open_dataset(nc_path)

    lat_name = "lat" if "lat" in ds.coords else "latitude"
    lon_name = "lon" if "lon" in ds.coords else "longitude"

    # Convert longitude from 0–360 to -180–180 if needed
    if float(ds[lon_name].max()) > 180:
        ds = ds.assign_coords(
            {lon_name: (((ds[lon_name] + 180) % 360) - 180)}
        ).sortby(lon_name)

    lat = ds[lat_name]
    lon = ds[lon_name]

    # IPCC AR6 land regions
    ar6 = regionmask.defined_regions.ar6.land

    print("Available AR6 regions:")
    for num, abbr, name in zip(ar6.numbers, ar6.abbrevs, ar6.names):
        print(f"{abbr:5s} {num:3d} {name}")

    # 2D mask: each grid cell has an AR6 region number
    mask = ar6.mask(lon, lat)

    abbrev_to_number = dict(zip(ar6.abbrevs, ar6.numbers))

    time_values = ds["time"].values if "time" in ds.coords else np.arange(ds.sizes["time"])

    records = []

    for region_name, abbrevs in region_groups.items():
        print(f"\nProcessing region group: {region_name}")
        print("  AR6 subregions:", abbrevs)

        region_numbers = []
        for abbr in abbrevs:
            if abbr not in abbrev_to_number:
                raise ValueError(f"AR6 abbreviation not found: {abbr}")
            region_numbers.append(abbrev_to_number[abbr])

        region_mask = xr.zeros_like(mask, dtype=bool)
        for n in region_numbers:
            region_mask = region_mask | (mask == n)

        if region_mask.sum().item() == 0:
            print(f"  [Warning] No grid cells found for {region_name}")
            continue

        temp = pd.DataFrame()
        temp["time"] = time_values
        temp["region"] = region_name

        for var in vars_to_sum:
            if var not in ds:
                print(f"  [Skip] variable not found: {var}")
                continue

            da = ds[var]

            if lat_name not in da.dims or lon_name not in da.dims:
                print(f"  [Skip] no spatial dims: {var}")
                continue

            ts = da.where(region_mask).sum(dim=[lat_name, lon_name], skipna=True)
            temp[var] = ts.values

        records.append(temp)

    if not records:
        raise ValueError("No regional results were generated.")

    out = pd.concat(records, ignore_index=True)
    out.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"\nSaved: {out_csv}")


REGION_COUNTRY_NAMES = {
    "North America": [
        "United States of America",
        "Canada",
        "Greenland",
    ],

    "Latin America": [
        "Mexico",
        "Guatemala",
        "Belize",
        "Honduras",
        "El Salvador",
        "Nicaragua",
        "Costa Rica",
        "Panama",
        "Cuba",
        "Haiti",
        "Dominican Rep.",
        "Jamaica",
        "Bahamas",
        "Puerto Rico",
        "Trinidad and Tobago",
        "Barbados",
        "Dominica",
        "Grenada",
        "Curaçao",
        "Colombia",
        "Venezuela",
        "Guyana",
        "Suriname",
        "Ecuador",
        "Peru",
        "Brazil",
        "Bolivia",
        "Paraguay",
        "Uruguay",
        "Argentina",
        "Chile",
        "Falkland Is.",
    ],

    "Europe": [
        "United Kingdom",
        "Ireland",
        "France",
        "Germany",
        "Netherlands",
        "Belgium",
        "Luxembourg",
        "Switzerland",
        "Austria",
        "Spain",
        "Portugal",
        "Italy",
        "Malta",
        "Denmark",
        "Norway",
        "Sweden",
        "Finland",
        "Iceland",
        "Estonia",
        "Latvia",
        "Lithuania",
        "Poland",
        "Czechia",
        "Slovakia",
        "Hungary",
        "Slovenia",
        "Croatia",
        "Bosnia and Herz.",
        "Serbia",
        "Montenegro",
        "Albania",
        "North Macedonia",
        "Greece",
        "Bulgaria",
        "Romania",
        "Moldova",
        "Ukraine",
        "Belarus",
        "Russia",
        "Kosovo",
        "Isle of Man",
        "Faeroe Is.",
        "Åland",
    ],

    "Middle East": [
        # Middle East
        "Turkey",
        "Cyprus",
        "N. Cyprus",
        "Georgia",
        "Armenia",
        "Azerbaijan",
        "Iran",
        "Iraq",
        "Syria",
        "Lebanon",
        "Israel",
        "Palestine",
        "Jordan",
        "Saudi Arabia",
        "Yemen",
        "Oman",
        "United Arab Emirates",
        "Qatar",
        "Kuwait",

        # North Africa
        "Morocco",
        "W. Sahara",
        "Algeria",
        "Tunisia",
        "Libya",
        "Egypt",

        # Central Asia
        "Kazakhstan",
        "Uzbekistan",
        "Turkmenistan",
        "Kyrgyzstan",
        "Tajikistan",
    ],

    "Sub-Saharan Africa": [
        "Mauritania",
        "Mali",
        "Niger",
        "Chad",
        "Sudan",
        "S. Sudan",
        "Eritrea",
        "Djibouti",
        "Somaliland",
        "Somalia",
        "Ethiopia",
        "Senegal",
        "Gambia",
        "Guinea-Bissau",
        "Guinea",
        "Sierra Leone",
        "Liberia",
        "Côte d'Ivoire",
        "Ghana",
        "Togo",
        "Benin",
        "Burkina Faso",
        "Nigeria",
        "Cabo Verde",
        "Cameroon",
        "Central African Rep.",
        "Eq. Guinea",
        "Gabon",
        "Congo",
        "Dem. Rep. Congo",
        "São Tomé and Principe",
        "Uganda",
        "Kenya",
        "Rwanda",
        "Burundi",
        "Tanzania",
        "Angola",
        "Zambia",
        "Malawi",
        "Mozambique",
        "Zimbabwe",
        "Botswana",
        "Namibia",
        "South Africa",
        "Lesotho",
        "eSwatini",
        "Madagascar",
        "Comoros",
        "Mauritius",
    ],

    "East Asia": [
        "China",
        "Taiwan",
        "Hong Kong",
        "Mongolia",
        "North Korea",
        "South Korea",
        "Japan",
    ],

    "South Asia": [
        "Afghanistan",
        "Pakistan",
        "India",
        "Bangladesh",
        "Nepal",
        "Bhutan",
        "Sri Lanka",
        "Siachen Glacier",
    ],

    "Southeast Asia": [
        "Myanmar",
        "Thailand",
        "Laos",
        "Cambodia",
        "Vietnam",
        "Malaysia",
        "Singapore",
        "Brunei",
        "Philippines",
        "Indonesia",
        "Timor-Leste",
    ],

    "Oceania": [
        "Australia",
        "New Zealand",
        "Papua New Guinea",
        "Solomon Is.",
        "Vanuatu",
        "Fiji",
        "New Caledonia",
        "Fr. Polynesia",
        "Samoa",
        "Micronesia",
        "Marshall Is.",
        "Palau",
        "Kiribati",
    ],
}

def extract_9region_timeseries(
    nc_path,
    out_csv,
    vars_to_sum=None,
    time_offset=None,
    include_unassigned=True,
):
    """
    Aggregate gridded NetCDF variables to 9 regions using Natural Earth country masks.

    The region definition is based on REGION_COUNTRY_NAMES, where countries are
    matched by Natural Earth country names, not ISO3 codes.
    """

    if vars_to_sum is None:
        vars_to_sum = [
            "Net_Emissions",
        ]

    ds = xr.open_dataset(nc_path)

    lat_name = "lat" if "lat" in ds.coords else "latitude"
    lon_name = "lon" if "lon" in ds.coords else "longitude"

    # Convert longitude from 0–360 to -180–180 if needed
    if float(ds[lon_name].max()) > 180:
        ds = ds.assign_coords(
            {lon_name: (((ds[lon_name] + 180) % 360) - 180)}
        ).sortby(lon_name)

    lat = ds[lat_name]
    lon = ds[lon_name]

    # Natural Earth country boundaries
    try:
        countries = regionmask.defined_regions.natural_earth_v5_0_0.countries_50
    except Exception:
        countries = regionmask.defined_regions.natural_earth_v5_0_0.countries_110

    mask = countries.mask(lon, lat)

    name_to_number = dict(zip(countries.names, countries.numbers))
    number_to_abbrev = dict(zip(countries.numbers, countries.abbrevs))
    number_to_name = dict(zip(countries.numbers, countries.names))

    # Check country names not found in current Natural Earth version
    all_target_names = sorted(set(sum(REGION_COUNTRY_NAMES.values(), [])))
    not_found_names = [
        name for name in all_target_names
        if name not in name_to_number
    ]

    if not_found_names:
        print("\n[Warning] These country names were not found in Natural Earth:")
        for name in not_found_names:
            print("  ", name)
        print("They will be skipped. Check spelling against Natural Earth names.\n")

    time_values = ds["time"].values if "time" in ds.coords else np.arange(ds.sizes["time"])

    records = []

    # This mask tracks all grid cells already assigned to one of the 9 regions
    assigned_region_mask = xr.zeros_like(mask, dtype=bool)

    for region_name, country_name_list in REGION_COUNTRY_NAMES.items():

        region_numbers = [
            name_to_number[name]
            for name in country_name_list
            if name in name_to_number
        ]

        if len(region_numbers) == 0:
            print(f"[Warning] No valid countries found for region: {region_name}")
            continue

        region_mask = mask.isin(region_numbers)

        # Very important: update assigned mask
        assigned_region_mask = assigned_region_mask | region_mask

        n_cells = int(region_mask.sum().item())
        print(f"{region_name:25s}: {n_cells:8d} grid cells")

        temp = pd.DataFrame()
        temp["time"] = time_values

        if time_offset is not None:
            temp["year"] = time_values + time_offset

        temp["region"] = region_name

        for var in vars_to_sum:
            if var not in ds:
                print(f"[Skip] variable not found: {var}")
                continue

            da = ds[var]

            if lat_name not in da.dims or lon_name not in da.dims:
                print(f"[Skip] no spatial dims: {var}")
                continue

            ts = da.where(region_mask).sum(dim=[lat_name, lon_name], skipna=True)
            temp[var] = ts.values

        records.append(temp)

    # Optional: output unassigned land cells to check global coverage
    if include_unassigned:
        country_land_mask = mask.notnull()

        # Exclude non-target special territories from unassigned check
        exclude_names = [
            "Antarctica",
            "Fr. S. Antarctic Lands",
            "Heard I. and McDonald Is.",
            "S. Geo. and the Is.",
        ]

        exclude_numbers = [
            name_to_number[name]
            for name in exclude_names
            if name in name_to_number
        ]

        if len(exclude_numbers) > 0:
            excluded_mask = mask.isin(exclude_numbers)
        else:
            excluded_mask = xr.zeros_like(mask, dtype=bool)

        unassigned_mask = country_land_mask & (~assigned_region_mask) & (~excluded_mask)

        n_unassigned = int(unassigned_mask.sum().item())
        print(f"{'Unassigned':25s}: {n_unassigned:8d} grid cells")

        if n_unassigned > 0:
            temp = pd.DataFrame()
            temp["time"] = time_values

            if time_offset is not None:
                temp["year"] = time_values + time_offset

            temp["region"] = "Unassigned"

            for var in vars_to_sum:
                if var not in ds:
                    continue

                da = ds[var]

                if lat_name not in da.dims or lon_name not in da.dims:
                    continue

                ts = da.where(unassigned_mask).sum(dim=[lat_name, lon_name], skipna=True)
                temp[var] = ts.values

            records.append(temp)

            # Print which countries are still unassigned
            unassigned_numbers = np.unique(mask.where(unassigned_mask).values)
            unassigned_numbers = [
                int(x) for x in unassigned_numbers
                if np.isfinite(x)
            ]

            print("\nUnassigned countries in Natural Earth mask:")
            for n in unassigned_numbers:
                print(f"  {number_to_abbrev.get(n, 'NA'):5s} {number_to_name.get(n, 'Unknown')}")

    if len(records) == 0:
        raise ValueError("No regional records were generated. Check REGION_COUNTRY_NAMES.")

    out = pd.concat(records, ignore_index=True)
    out.to_csv(out_csv, index=False, encoding="utf-8-sig")

    print(f"\nSaved: {out_csv}")
    
# if __name__ == "__main__":
    # extract_ar6_region_groups(
        # nc_path=r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\summary_025deg.global_time_ge1000_7.9.nc",
        # out_csv="regional_AR6_ELUC_timeseries.csv",
    # )
    
if __name__ == "__main__":
    extract_9region_timeseries(
        nc_path=r"D:\Work\Research_doc\LULCC\BookKeeping\Out_ncfile\summary_025deg.global_time_ge1000_7.9.nc",
        out_csv="ELUC_9region_timeseries.csv",
        time_offset=850,   # 如果 time=1000 表示 1850，就保留850；如果time已经是年份，改成 None
        include_unassigned=True,
    )