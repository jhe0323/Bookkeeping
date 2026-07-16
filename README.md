# LULCC Carbon Bookkeeping Model README

## 1. Overview

This model is a Python-based land-use and land-cover change (LULCC) carbon bookkeeping model. It tracks annual carbon stock changes and LULCC-induced carbon fluxes caused by land-use transitions, including clearing, abandonment, wood harvest, and other land transitions.

The current version has been updated to support global simulations at multiple spatial resolutions, including 1°, 0.5°, and 0.25°. The model is designed to be broadly consistent with the bookkeeping logic used in BLUE, LUCE, and Global Carbon Budget (GCB)-style land-use emission accounting.

**Important requirement:** all input files must use the same spatial resolution, coordinate system, grid definition, and latitude/longitude orientation after preprocessing. For example, 0.25° LUH2 states, 0.25° LUH2 transitions, and 0.25° PFT maps should be used together.

---

## 2. Main Features of the Current Version

The current model version includes the following main features:

- Supports global or regional simulations using user-defined latitude and longitude bounding boxes.
- Supports 1°, 0.5°, and 0.25° input data, provided all input files are spatially consistent.
- Uses LUH2 state and transition files to derive annual land-cover fractions and transition areas.
- Uses static PFT maps to assign PFT-dependent biomass and soil carbon densities.
- Uses a YAML configuration file to store carbon densities, transition parameters, decay parameters, and allocation coefficients.
- Calculates carbon changes for biomass, soil, product pools, and atmosphere.
- Tracks LULCC fluxes using event-based bookkeeping processes.
- Produces annual gridded NetCDF outputs suitable for global budget checks, spatial mapping, and comparison with GCB/BLUE/LUCE results.
- Includes GCB-style outputs such as gross sources, gross sinks, and net emissions.
- Includes LUCE-style LULCC categories, including deforestation, non-forest conversion, reforestation, non-forest reconstruction, conversions between anthropogenic land, and wood harvest products.

---

## 3. Input Data

### 3.1 LUH2 Land-Cover State File

The LUH2 state file provides annual land-cover fractions for each grid cell.

Typical file:

```text
states.nc
```

The LUH2 sub-categories are aggregated into the model land-cover classes used for bookkeeping.

### 3.2 LUH2 Transition File

The LUH2 transition file provides annual transition fractions between land-cover classes.

Typical file:

```text
transitions.nc
```

The model converts transition fractions into transition areas by multiplying the transition fraction by grid-cell area.

### 3.3 PFT Map

The current workflow supports static PFT maps aligned to the LUH2 grid. Recent simulations use an IBIS-derived dominant PFT map converted to the target resolution, for example:

```text
IBIS_PFT_dominant_0.25deg.nc
```

The PFT map should be preprocessed so that:

- latitude and longitude match the internal LUH2 grid;
- longitude is converted to the model convention, normally `-180 to 180`;
- latitude orientation is handled correctly;
- PFT fractions sum to 1 over valid land cells;
- missing or invalid PFT cells are treated consistently.

Earlier versions used ORCHIDEE PFT maps. The current README assumes the newer IBIS/static dominant-PFT workflow, but the model can still use other PFT maps if they are converted into the expected format.

### 3.4 Parameter YAML File

Typical file:

```text
config.yml
```

The YAML file stores:

- PFT-specific biomass carbon density;
- PFT-specific soil carbon density;
- land-cover and PFT mapping information;
- clearing allocation parameters;
- abandonment/regrowth parameters;
- harvest allocation parameters;
- product-pool decay parameters;
- soil and slash decay parameters;
- model output variable settings.

Carbon density values in the YAML file are generally defined as areal carbon density, for example `t C ha^-1`. The model converts these values internally when calculating grid-cell total carbon stocks.

---

## 4. Land-Cover Classes

The model aggregates LUH2 sub-categories into the following major land-cover classes:

```text
Pr  = Primary vegetation
S   = Secondary vegetation
Pa  = Pasture/rangeland
C   = Cropland
U   = Urban land, read from LUH2 but normally not explicitly simulated as an active carbon-transition class
```

The active bookkeeping calculation mainly uses:

```text
Pr, S, Pa, C
```

Urban land is retained for completeness of land-cover accounting but is not currently treated as a full active LULCC event class in the main carbon-transition calculation.

---

## 5. PFT Classes

The model uses 15 PFT classes. Each PFT has corresponding biomass and soil carbon density parameters in `config.yml`.

The PFT map can be stored as either:

- a dominant PFT class per grid cell; or
- a PFT-fraction map with dimensions such as `(time, veget, lat, lon)`.

The current 0.25° workflow uses a preprocessed IBIS-derived PFT map with 15 PFT classes. The file should already be aligned to the target LUH2 grid before running the model.

---

## 6. Carbon Pools

The model tracks the following carbon pools:

```text
Biomass
Soil
Product 1-year pool
Product 10-year pool
Product 100-year pool
Atmosphere
```

In the annual NetCDF output, these pools may appear as variables such as:

```text
biomass_total
soil_total
P1
P10
P100
atmo or atmosphere
```

Depending on the summary script version, additional flux variables are also written for model diagnostics and GCB/LUCE comparison.

---

## 7. LULCC Event Types

The current model recognizes the following event groups:

### 7.1 Clearing

Clearing includes transitions from forest or natural vegetation to agricultural land, such as:

```text
Pr -> C
Pr -> Pa
S  -> C
S  -> Pa
```

Clearing carbon losses are allocated to:

- direct atmospheric emission;
- product pools;
- soil/slash pools;
- remaining biomass and soil changes after transition.

### 7.2 Abandonment / Regrowth

Abandonment and regrowth include transitions such as:

```text
C  -> S
Pa -> S
Pr -> S
S  -> S
```

These transitions are used to represent vegetation recovery, secondary forest regrowth, and related soil-carbon adjustment.

### 7.3 Harvest

Harvest is treated as a wood-removal process mainly associated with managed vegetation transitions such as:

```text
Pr -> S
S  -> S
```

The current harvest logic distinguishes between:

- requested harvest biomass;
- harvest biomass that can be met from available ecosystem biomass;
- unmet harvest biomass;
- biomass transferred to wood-product pools;
- biomass transferred to soil/slash pools;
- subsequent product and slash decay.

The latest harvest formulation is intended to avoid non-physical over-removal of biomass and to maintain carbon bookkeeping consistency when requested harvest exceeds available biomass.

### 7.4 Other Transitions

Other transitions include non-forest and anthropogenic land conversions such as:

```text
Pa -> C
C  -> Pa
```

These transitions are handled separately from forest-clearing and regrowth events.

---

## 8. LUCE-Style Output Categories

For comparison with LUCE-style LULCC accounting, the model can aggregate fluxes into the following six categories:

```text
Flux_FD   = Deforestation
Flux_NFC  = Non-forest conversion
Flux_FR   = Reforestation and harvest-related regrowth
Flux_NFR  = Non-forest reconstruction
Flux_CAL  = Conversions between anthropogenic land
Flux_WHp  = Wood harvest products
```

These categories are intended for high-level reporting and comparison with published LULCC emission products.

---

## 9. GCB-Style Budget Variables

For Global Carbon Budget-style analysis, the model can produce:

```text
Gross_Sources
Gross_Sinks
Net_Emissions
```

where:

```text
Net_Emissions = Gross_Sources + Gross_Sinks
```

By convention:

- positive values represent carbon emissions to the atmosphere;
- negative values represent carbon uptake or sink terms;
- net emissions combine both source and sink components.

These variables are useful for checking annual global budgets and comparing the model with BLUE/GCB-style outputs.

---

## 10. Important Diagnostic Variables

Depending on the current `summary_yearly.py` configuration, the model may also output diagnostic variables such as:

```text
unmet_harvest
area_clearing
area_abandonment
area_other
area_harvest
clearing_to_A
clearing_to_products
clearing_to_soil
harvest_to_products
harvest_to_soil
harvest_biomass_removed
abd_delta_b
abd_delta_ss
emit_clearing
emit_abandonment
emit_other
emit_harvest
emit_products
emit_abandonment_biomass
emit_abandonment_soil
area_deforestation
deforestation_to_A
deforestation_to_products
deforestation_to_soil
emit_products_clearing
emit_products_harvest
harvest_requested_biomass
harvest_met_biomass
harvest_unmet_biomass
harvest_beta_delta_sum
harvest_delta_B_h
harvest_R
harvest_delta_SS_h
harvest_unmet_to_products
harvest_unmet_to_soil
```

These variables are primarily used for debugging, process diagnosis, and closure checks. For final comparison with GCB/BLUE/LUCE products, the high-level variables such as `Gross_Sources`, `Gross_Sinks`, `Net_Emissions`, and LUCE category fluxes are usually more appropriate.

---

## 11. Units

### 11.1 Internal Units

The model internally calculates carbon stocks and fluxes using grid-cell area and carbon-density parameters.

If carbon densities in `config.yml` are given as:

```text
t C ha^-1
```

then the model converts them to total grid-cell carbon using:

```text
carbon_stock = area_m2 / 10000 * carbon_density_tC_ha
```

The resulting grid-cell carbon stock is in:

```text
t C
```

### 11.2 Output Units

For gridded NetCDF outputs, raw values are usually stored in native total-carbon units unless converted in a post-processing or plotting script.

Common conversions are:

```text
Tg C     = t C / 1e6
Gt C     = t C / 1e9
Tg C yr^-1 or Gt C yr^-1 = carbon-stock difference divided by number of years
```

For area-normalized flux maps:

```text
g C m^-2 yr^-1 = total_flux_tC_yr * 1e6 / area_m2
```

When plotting annual mean emissions or sinks over a period, the recommended workflow is:

```text
period_mean_flux = (carbon_stock_end - carbon_stock_start) / number_of_years
```

and then convert to either `Gt C yr^-1`, `Tg C yr^-1`, or `g C m^-2 yr^-1` depending on the purpose.

---

## 12. Recommended Simulation Strategy

For historical LULCC emission analysis, it is recommended to run the model from the base year (850) than only from the target analysis period.

For example, to analyze 1960-2020 emissions, a recommended setup is:

```text
base year: 850
analysis period: 1960-2020
```

Running from an early base year allows product pools, soil pools, slash pools, and regrowth cohorts to carry historical memory. This is important because LULCC carbon fluxes are not only determined by current-year transitions but also by legacy effects from past land-use changes.

Short runs starting too close to the target period may smooth or distort annual variability, especially for product-pool decay, soil adjustment, and regrowth effects.

---

## 13. Computation Workflow

The model workflow is:

```text
1. Read config.yml.
2. Read LUH2 states.nc.
3. Read LUH2 transitions.nc.
4. Build the internal latitude/longitude grid.
5. Calculate grid-cell area.
6. Read and align the PFT map to the internal grid.
7. Select simulation domain using a bounding box or grid slices.
8. Initialize biomass and soil pools for each grid cell.
9. Parse annual LUH2 transitions.
10. Allocate each transition to a LULCC event type.
11. Update biomass, soil, product, and atmosphere pools.
12. Apply product-pool and soil/slash decay.
13. Save annual gridded outputs to NetCDF.
14. Merge rank/band outputs if the model is run in parallel.
15. Use post-processing scripts to calculate annual means, period differences, maps, and global totals.
```

---

## 14. Main Code Structure

A typical project structure is:

```text
<ProjectRoot>/
├── In_ncfile/
│   ├── states_025deg.nc
│   ├── transitions_025deg.nc
│   └── IBIS_PFT_dominant_0.25deg.nc
│
├── Out_ncfile/
│   ├── summary_025deg.global.nc
│   ├── summary_025deg.rank000.nc
│   ├── summary_025deg.rank001.nc
│   └── ...
│
├── logs/
│   └── slurm or runtime logs
│
└── v2/
    ├── main.py
    ├── config.yml
    ├── LULCCSimulator.py
    ├── file_loader.py
    ├── parameter_loader.py
    ├── carbon_pools_init.py
    ├── transition.py
    ├── summary_yearly.py
    ├── events.py
    ├── merge025.py
    └── plotting / diagnostic scripts
```

The exact file names can differ between local Windows tests and Linux cluster runs.

---

## 15. Main Modules

### 15.1 `main.py`

Main entry point of the model. It defines:

- input paths;
- output paths;
- simulation period;
- spatial domain;
- model resolution;
- cluster rank/band settings if used;
- call to `LULCCSimulator`.

### 15.2 `LULCCSimulator.py`

Core simulation controller. Main responsibilities include:

- loading input data;
- building the internal grid;
- calculating grid-cell area;
- mapping internal grid indices to LUH2 dataset indices;
- reading initial land-cover areas;
- parsing annual transitions;
- running grid-cell simulations;
- writing NetCDF outputs.

### 15.3 `file_loader.py`

Handles LUH2 and PFT input files, including:

- reading NetCDF variables;
- checking latitude orientation;
- checking longitude convention;
- aligning PFT maps to the internal grid;
- supporting different spatial resolutions.

### 15.4 `parameter_loader.py`

Reads `config.yml` and exposes model parameters to the simulator.

### 15.5 `carbon_pools_init.py`

Initializes biomass and soil carbon pools based on:

- grid-cell area;
- land-cover fractions;
- PFT-specific carbon densities.

### 15.6 `events.py`

Defines event-level allocation logic for:

- clearing;
- abandonment;
- harvest;
- other land transitions.

### 15.7 `transition.py`

Defines process-level updates such as:

- biomass removal;
- biomass regrowth;
- soil adjustment;
- harvest residue transfer;
- product-pool transfer;
- product-pool decay;
- slash/soil decay.

### 15.8 `summary_yearly.py`

Aggregates annual outputs, including:

- pool totals;
- atmosphere fluxes;
- diagnostic variables;
- GCB-style budget variables;
- LUCE-style category fluxes.

---

## 16. Current Run Logic

A simplified call structure is:

```text
main.py
│
├─ set BASE_DIR / PROJECT_DIR / DATA_DIR / OUT_DIR
│
└─ if __name__ == "__main__":
   │
   ├─ sim = LULCCSimulator(config_path, state_path, trans_path, pft_path)
   │    │
   │    ├─ ParameterLoader(config.yml)
   │    ├─ FileLoader.load_luh2_dataset(states.nc)
   │    ├─ FileLoader.load_luh2_dataset(transitions.nc)
   │    ├─ LULCCSimulator._build_grid_and_area()
   │    └─ FileLoader.load_pft_map(...)
   │
   ├─ lat_slice, lon_slice = sim.indices_for_bbox(lat_min, lat_max, lon_min, lon_max)
   │
   └─ sim.run_simulation_grid(years, out_nc, lat_slice, lon_slice, start_year_idx)
        │
        ├─ create NetCDF output file
        ├─ write lat/lon/time dimensions
        ├─ create output variables
        └─ loop over grid cells
             │
             └─ yearly = sim.run_simulation(years, i, j, start_year_idx)
                  │
                  ├─ get cell PFT
                  ├─ read initial LULC areas
                  ├─ initialize carbon pools
                  ├─ for each simulation year:
                  │    ├─ parse LUH2 transitions
                  │    ├─ classify transitions into LULCC events
                  │    ├─ update biomass/soil/product/atmosphere pools
                  │    ├─ apply decay and regrowth
                  │    └─ summarize annual outputs
                  └─ return yearly outputs
```

---

## 17. Running Global 0.25° Simulations on a Cluster

For 0.25° global simulations, the full grid is large. The recommended workflow is to split the latitude dimension into multiple bands or ranks.

Example output pattern:

```text
summary_025deg.rank000.nc
summary_025deg.rank001.nc
summary_025deg.rank002.nc
...
summary_025deg.rank119.nc
```

After all ranks finish, merge them into a global file:

```text
summary_025deg.global.nc
```

When checking whether a rank finished correctly, compare file sizes and inspect whether the `done` variable is complete. Very small files usually indicate that the rank failed before writing valid results.

---

## 18. NetCDF Output

The output NetCDF file usually contains:

```text
time
lat
lon
```

and variables such as:

```text
biomass_total
soil_total
P1
P10
P100
atmo or atmosphere
done
Gross_Sources
Gross_Sinks
Net_Emissions
Flux_FD
Flux_NFC
Flux_FR
Flux_NFR
Flux_CAL
Flux_WHp
```

The exact variable list depends on the current `summary_yearly.py` and output configuration.

---

## 19. Post-Processing and Plotting

Typical post-processing tasks include:

### 19.1 Extracting a Period Mean Flux

For a variable stored as cumulative or pool-based annual values:

```text
mean_flux = (value_end - value_start) / number_of_years
```

For example, for 1960-2020:

```text
mean_flux_1960_2020 = (value_2020 - value_1960) / 60
```

### 19.2 Converting to Global Total Units

```text
Gt C yr^-1 = t C yr^-1 / 1e9
Tg C yr^-1 = t C yr^-1 / 1e6
```

### 19.3 Converting to Area-Normalized Flux

```text
g C m^-2 yr^-1 = t C yr^-1 * 1e6 / area_m2
```

### 19.4 Recommended Map Convention

For ELUC or atmosphere flux maps:

- positive values = emissions/source;
- negative values = sink/uptake;
- zero values can be shown as white;
- diverging color maps should only be used when both positive and negative values exist;
- if all values are positive, use a sequential emission color scale;
- if all values are negative, use a sequential sink color scale.

---

## 20. Common Checks Before Running

Before running the model, check:

```text
1. states.nc, transitions.nc, and PFT map have the same resolution.
2. Latitude orientation is handled correctly.
3. Longitude convention is consistent, normally -180 to 180.
4. PFT map is aligned with the internal LUH2 grid.
5. PFT fractions sum to 1 over valid land cells.
6. Carbon density units in config.yml are known and correctly converted.
7. start_year_idx matches the intended calendar year.
8. Output time dimension corresponds to the intended simulation period.
9. The run starts early enough to include legacy effects.
10. NetCDF variables are not all NaN or all zero after a test run.
```

---

## 21. Interpreting `start_year_idx`

LUH2 time indexing depends on the base year of the input file. For standard LUH2 files starting in year 850:

```text
calendar_year = 850 + time_index
```

For example:

```text
start_year_idx = 0     -> year 850
start_year_idx = 1100  -> year 1950
start_year_idx = 1152  -> year 2002
```

Always confirm the `time` variable in the NetCDF file before running production simulations.

---

## 22. Recommended Validation Workflow

Recommended validation steps are:

```text
1. Run a single-grid synthetic test.
2. Check hand-calculated carbon conservation for simple transitions.
3. Run a small regional test.
4. Check that no output variables are unexpectedly NaN.
5. Check annual global totals.
6. Compare deforestation-only emissions with BLUE/GCB-style references.
7. Check gross source, gross sink, and net emission closure.
8. Inspect spatial maps for unrealistic hotspots.
9. Diagnose clearing, harvest, product-pool, and soil-pool contributions separately.
10. Run full global simulations only after the above checks pass.
```

---

## 23. Known Notes and Cautions

- The model is sensitive to biomass and soil carbon density parameters.
- Clearing emissions can become very large if biomass carbon density is too high for certain PFT/land-cover combinations.
- Harvest calculations require special care when requested harvest exceeds available biomass.
- Short initialization periods can distort product-pool and soil-memory effects.
- PFT and LUH2 grid mismatch is one of the most common sources of spatial artifacts.
- Latitude flipping can cause north-south inversion in output maps if not handled consistently.
- For final publication-quality results, diagnostic variables should be used to verify process consistency before relying on aggregated budget variables.

---

## 24. Minimal Example

A typical Python entry point may look like:

```python
from LULCCSimulator import LULCCSimulator

config_path = "config.yml"
state_path = "In_ncfile/states_025deg.nc"
trans_path = "In_ncfile/transitions_025deg.nc"
pft_path = "In_ncfile/IBIS_PFT_dominant_0.25deg.nc"
out_nc = "Out_ncfile/summary_025deg.global.nc"

sim = LULCCSimulator(config_path, state_path, trans_path, pft_path)

lat_slice, lon_slice = sim.indices_for_bbox(
    lat_min=-90,
    lat_max=90,
    lon_min=-180,
    lon_max=180,
)

sim.run_simulation_grid(
    years=1171,
    out_nc=out_nc,
    lat_slice=lat_slice,
    lon_slice=lon_slice,
    start_year_idx=0,
)
```

This example assumes a LUH2 file starting in year 850 and running through year 2020. Adjust `years` and `start_year_idx` according to the actual time dimension of the input file.

---

## 25. Suggested Citation / Method Description

When describing the model in a manuscript or report, the method can be summarized as:

```text
We used an event-based LULCC carbon bookkeeping model driven by LUH2 annual land-cover states and transition matrices. The model tracks biomass, soil, product-pool, and atmospheric carbon changes for each grid cell. Land-use transitions were classified into clearing, abandonment/regrowth, harvest, and other transition events, and carbon transfers were calculated using PFT-specific biomass and soil carbon densities and event-specific allocation and decay parameters. Annual gridded outputs were aggregated into gross sources, gross sinks, net LULCC emissions, and LUCE-style transition categories for comparison with existing bookkeeping-based LULCC emission products.
```