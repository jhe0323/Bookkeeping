# LULCC Carbon Bookkeeping Model

A Python-based land-use and land-cover change (LULCC) carbon bookkeeping model for estimating annual carbon-stock changes and land-use change emissions. The model is driven by LUH2 land-cover states and transitions, PFT-dependent carbon densities, and event-specific allocation and decay parameters.

The code tracks carbon changes associated with:

- land clearing;
- agricultural abandonment and regrowth;
- wood harvest;
- cropland–pasture conversion;
- product-pool decay;
- biomass and soil-carbon relaxation.

The model is designed for comparison with BLUE-, LUCE-, and Global Carbon Budget-style bookkeeping estimates.

---

## 1. Repository branches

The repository contains three model configurations. They share the same core bookkeeping framework and differ mainly in the treatment of wood harvest.

| Branch | Harvest forcing | Harvest area constraint | Main interpretation |
|---|---|---|---|
| `main` | LUH2 `*_harv` area only | Strictly limited to the LUH2 harvest footprint | Pure area-driven harvest |
| `harvest-area-forced` | LUH2 `*_bioh` biomass demand | Strictly limited to the LUH2 harvest footprint | Biomass-demand-driven harvest with unmet demand |
| `harvest-bio-forced` | LUH2 `*_bioh` biomass demand | May expand beyond the LUH2 footprint within the same source cover and allowed PFT group | Biomass-demand-driven harvest with forced area expansion |

Switch between versions with:

```bash
git switch main
git switch harvest-area-forced
git switch harvest-bio-forced
```

### 1.1 `main`: pure LUH2-area-driven harvest

The default branch uses only LUH2 `*_harv` variables.

Main rules:

1. LUH2 `*_bioh` is ignored by the harvest calculation.
2. Harvested biomass is diagnosed from the model state on the LUH2 harvest area.
3. No harvest expansion outside the LUH2 footprint is allowed.
4. The area-implied biomass removal is allocated to wood-product and rapid-soil pools.

Conceptually:

```text
removed biomass
= current modeled biomass on source cover
× harvested share of the source-cover area
```

In this branch:

```text
harvest_requested_biomass = harvest_met_biomass
harvest_biomass_removed   = harvest_met_biomass
harvest_forced_biomass    = 0
harvest_unmet_biomass     = 0
harvest_extra_area        = 0
```

### 1.2 `harvest-area-forced`: strict-area biomass-demand harvest

This branch uses LUH2 `*_bioh` as the requested biomass removal but does not allow harvest outside the LUH2 `*_harv` footprint.

```text
requested biomass = met biomass + unmet biomass
removed biomass   = met biomass
forced biomass    = 0
extra area        = 0
```

Any biomass demand that cannot be supplied from the LUH2 harvest footprint is retained as `harvest_unmet_biomass` and is not added to product or soil pools.

### 1.3 `harvest-bio-forced`: biomass-demand harvest with area expansion

This branch uses LUH2 `*_bioh` as the requested biomass removal. The LUH2 harvest area is treated as the preferred footprint. If biomass within that footprint is insufficient, harvest may expand within the same source cover, harvest family, and allowed PFT group.

```text
requested biomass
= met biomass + forced biomass + unmet biomass

removed biomass
= met biomass + forced biomass
```

Additional removal outside the LUH2 footprint is recorded as:

```text
harvest_forced_biomass
harvest_extra_area
```

---

## 2. Core model features

- Global or regional simulations.
- Support for 1°, 0.5°, and 0.25° grids when all input files are spatially consistent.
- LUH2 state and transition forcing.
- Static or transient PFT-dependent biomass and soil carbon densities.
- Event-based carbon bookkeeping.
- Explicit biomass, slow-soil, rapid-soil, product, and atmospheric pools.
- Historical legacy effects from past land-use events.
- GCB-style gross sources, gross sinks, and net emissions.
- LUCE-style process categories.
- Gridded NetCDF output.
- Diagnostic variables for harvest, clearing, abandonment, and carbon closure.

---

## 3. Land-cover classes

The model aggregates LUH2 land-use classes into five internal categories:

| Internal code | Meaning | Typical LUH2 classes |
|---|---|---|
| `v` | Primary vegetation | `primf`, `primn` |
| `s` | Secondary vegetation | `secdf`, `secdn` |
| `p` | Pasture and rangeland | `pastr`, `range` |
| `c` | Cropland | `c3ann`, `c4ann`, `c3per`, `c4per`, `c3nfx` |
| `U` | Urban land | `urban` |

The active bookkeeping calculations mainly use `v`, `s`, `p`, and `c`. Urban land is retained for land-cover accounting but is not currently implemented as a complete active event class.

---

## 4. Event classification

### 4.1 Clearing

```text
v -> c
v -> p
s -> c
s -> p
```

Clearing transfers carbon among atmosphere, product pools, rapid soil/slash, slow soil, and biomass-deficit pools.

### 4.2 Abandonment and regrowth

```text
c -> s
p -> s
```

These transitions inherit the source carbon state into secondary vegetation and generate biomass and soil recovery trajectories.

### 4.3 Wood harvest

Wood harvest is processed separately from ordinary state-to-state transitions using LUH2 harvest families:

```text
primf
primn
secmf
secyf
secnf
```

Primary harvest is represented as a transfer from primary to secondary vegetation. Secondary harvest remains within secondary vegetation but creates a new harvest-history carbon deficit.

### 4.4 Other transitions

```text
c -> p
p -> c
```

These transitions are handled separately from clearing and abandonment.

---

## 5. Carbon pools

The equilibrium carbon array stores:

```text
B   biomass carbon
SS  slow soil carbon
```

The excess or legacy carbon array stores:

```text
B     biomass departure from equilibrium
SS    slow-soil departure from equilibrium
SR    rapid soil/slash carbon
P1    1-year product pool
P10   10-year product pool
P100  100-year product pool
A     direct atmospheric release
```

Historical departures are tracked separately for:

```text
l  clearing
h  harvest
a  abandonment
g  other transitions
```

This allows past land-use events to continue affecting present-day fluxes.

---

## 6. PFT and carbon-density configuration

The current model uses 15 PFT classes. PFT-dependent parameters are read from `config.yml`.

The PFT input may be:

- a 15-class PFT fraction map; or
- a dominant-PFT map converted to the expected grid and class convention.

All PFT files must match the LUH2 input grid in spatial resolution, coordinates, latitude orientation, longitude convention, land mask, and PFT numbering.

### Static and transient carbon density

The model supports both static and transient carbon densities through:

```yaml
dynamic_carbon_density:
  enabled: true
  path: dynamic_carbon_density.parquet
  base_year: 850
```

Set:

```yaml
enabled: false
```

to use the static values stored directly in `config.yml`.

When transient density is enabled, the Parquet file must contain year- and PFT-specific biomass and soil carbon densities for `v`, `s`, `p`, and `c`.

---

## 7. Required input data

Typical inputs are:

```text
states.nc
transitions.nc
PFT_map.nc
config.yml
dynamic_carbon_density.parquet   # required only when transient density is enabled
```

All gridded inputs must use the same resolution, coordinate system, grid-cell centers, latitude ordering, longitude convention, and time indexing.

The current `main.py` defines the input and output paths. Update those paths before running the model.

A common local layout is:

```text
<ProjectRoot>/
├── In_ncfile/
│   ├── states_1deg.nc
│   ├── transitions_1deg.nc
│   └── PFTmap_orchidee_1deg.nc
│
├── Out_ncfile/
│
└── Bookkeeping/
    ├── main.py
    ├── LULCCSimulator.py
    ├── events.py
    ├── transition.py
    ├── file_loader.py
    ├── parameter_loader.py
    ├── carbon_pools_init.py
    ├── summary_yearly.py
    ├── config.yml
    ├── dynamic_carbon_density.parquet
    ├── requirements.txt
    ├── README.md
    └── tools/                       # retained only on main
```

The `tools/` directory contains preprocessing, checking, merging, plotting, or diagnostic utilities. It is retained only on `main`; the two harvest experiment branches contain only the model files required for their respective simulations.

---

## 8. Main code modules

| File | Role |
|---|---|
| `main.py` | Defines paths, simulation period, spatial domain, and output file |
| `LULCCSimulator.py` | Controls loading, initialization, yearly simulation, diagnostics, and NetCDF writing |
| `events.py` | Implements clearing, abandonment, harvest, and other transition events |
| `transition.py` | Applies annual pool relaxation, decay, regrowth, and process-level flux accounting |
| `file_loader.py` | Loads and aligns LUH2 and PFT input files |
| `parameter_loader.py` | Reads `config.yml` and optional transient carbon-density data |
| `carbon_pools_init.py` | Defines indices and initializes equilibrium and excess carbon pools |
| `summary_yearly.py` | Summarizes annual stocks and atmospheric carbon |
| `config.yml` | Stores land-cover mapping, carbon densities, allocation coefficients, and response times |
| `tools/` | Utility scripts; available only on `main` |

The main call sequence is:

```text
main.py
└── LULCCSimulator
    ├── ParameterLoader
    ├── FileLoader
    ├── carbon_pools_init
    ├── events
    ├── transition
    └── summary_yearly
```

---

## 9. Installation

Python 3.8 or later is recommended.

```bash
pip install -r requirements.txt
```

Core dependencies include:

```text
numpy
netCDF4
PyYAML
pandas
pyarrow
```

`pandas` and `pyarrow` are required when transient carbon density is read from a Parquet file.

---

## 10. Running the model

The recommended entry point is:

```bash
python main.py
```

A minimal Python example is:

```python
from pathlib import Path
from LULCCSimulator import LULCCSimulator

base_dir = Path(__file__).resolve().parent
project_dir = base_dir.parent

sim = LULCCSimulator(
    config_path=str(base_dir / "config.yml"),
    LULC_path=str(project_dir / "In_ncfile" / "states_1deg.nc"),
    trans_path=str(project_dir / "In_ncfile" / "transitions_1deg.nc"),
    pft_path=str(project_dir / "In_ncfile" / "PFTmap_orchidee_1deg.nc"),
    lat_slice=None,
    lon_slice=None,
    area_unit="ha",
)

sim.run_simulation_grid(
    years=1172,
    out_nc=str(project_dir / "Out_ncfile" / "summary_global.nc"),
    lat_slice=None,
    lon_slice=None,
    start_year_idx=0,
    sync_every=200,
)
```

### Important meaning of `years`

`years` is the number of annual transitions simulated.

```text
number of stored states = years + 1
```

For example, a simulation containing states from 850 through 2022 has:

```text
1173 stored states
1172 annual transitions
```

Therefore:

```python
years = 1172
```

Always verify the actual time dimensions of `states.nc` and `transitions.nc` before a production run.

---

## 11. Spatial subsets and band runs

`lat_slice` and `lon_slice` can be passed to `LULCCSimulator` to load and simulate only part of the global grid.

```python
sim = LULCCSimulator(
    config_path="config.yml",
    LULC_path="states.nc",
    trans_path="transitions.nc",
    pft_path="PFT_map.nc",
    lat_slice=slice(80, 100),
    lon_slice=slice(120, 140),
)
```

For global 0.25° simulations, divide the domain into longitude or latitude bands, write one NetCDF file per band, and merge the outputs after all bands finish.

---

## 12. Main output variables

### Carbon stocks

```text
biomass_total
soil_total
P1
P10
P100
atmosphere
```

### GCB-style budget variables

```text
Gross_Sources
Gross_Sinks
Net_Emissions
Closure_Error
```

Sign convention:

```text
positive = atmospheric source
negative = atmospheric sink
```

The model uses:

```text
Net_Emissions = Gross_Sources + Gross_Sinks
```

### LUCE-style categories

```text
Flux_FD    deforestation
Flux_NFC   non-forest conversion
Flux_FR    forest regrowth
Flux_NFR   non-forest reconstruction
Flux_CAL   conversion between agricultural land classes
Flux_WHp   wood-product emissions
```

### Event-level fluxes

```text
Flux_Clearing
Flux_Abandonment
Flux_Harvest_Net
Flux_Harvest_SoilSlash
Flux_Harvest_Regrowth
Flux_Other
Flux_Products_Total
```

---

## 13. Harvest diagnostics

The branches use the same diagnostic names so that results can be compared directly.

### Demand and removal

```text
harvest_requested_biomass
harvest_met_biomass
harvest_forced_biomass
harvest_unmet_biomass
harvest_unmet_raw_biomass
harvest_biomass_removed
harvest_loss_biomass
harvest_patch_biomass_before
```

### Area

```text
harvest_luh2_area_frac
harvest_luh2_area
harvest_effective_area_frac
harvest_effective_area
harvest_extra_area
area_harvest
```

### Allocation

```text
harvest_to_products
harvest_to_soil
harvest_to_SR_from_biomass
harvest_to_SR_from_soil
harvest_forced_to_products
harvest_forced_to_soil
```

### Equation-level diagnostics

```text
harvest_beta_delta_sum
harvest_sigma_delta_sum
harvest_delta_B_h
harvest_delta_SS_h
harvest_R
```

Variables that are not applicable to a branch are retained and reported as zero, allowing NetCDF files from different schemes to be compared using the same analysis scripts.

---

## 14. Units

When `area_unit="ha"` and carbon densities are stored as:

```text
t C ha^-1
```

grid-cell carbon stocks and annual fluxes are stored in:

```text
t C
```

Common conversions are:

```text
Tg C = t C / 1e6
Pg C = t C / 1e9
```

For area-normalized maps:

```text
g C m^-2 yr^-1
= total flux in t C yr^-1 × 1e6 / grid-cell area in m²
```

Check the units of every input data source before changing `area_unit`.

---

## 15. Recommended validation workflow

Before running a full global simulation:

1. Check that states, transitions, and PFT maps use the same grid.
2. Verify latitude orientation and longitude convention.
3. Confirm that PFT fractions sum to 1 over valid land cells.
4. Confirm that `years` and `start_year_idx` remain within the forcing time dimension.
5. Run a single-grid synthetic test.
6. Test clearing, abandonment, harvest, and crop–pasture conversion separately.
7. Check biomass, soil, product, atmosphere, and total-system carbon closure.
8. Verify that no output variables are unexpectedly all zero or all `NaN`.
9. Compare global annual gross sources, gross sinks, and net emissions.
10. Inspect harvest requested, met, forced, unmet, and area diagnostics.
11. Inspect spatial maps for coordinate reversal or unrealistic hotspots.
12. Run the full global simulation only after these checks pass.

Recommended harvest identities:

```text
main:
requested = met = removed
forced = unmet = 0

harvest-area-forced:
requested = met + unmet
removed = met
forced = 0

harvest-bio-forced:
requested = met + forced + unmet
removed = met + forced
```

---

## 16. Known cautions

- Carbon-density assumptions strongly affect the magnitude of LULCC emissions.
- Static and transient carbon-density simulations should be clearly distinguished.
- LUH2 `*_harv` and `*_bioh` represent different constraints and should not be treated as interchangeable.
- Harvest biomass must not be added to product or soil pools unless it has been removed from the modeled ecosystem stock.
- PFT and LUH2 grid misalignment can create large spatial artifacts.
- A reversed latitude axis can invert the spatial result.
- Short simulations omit legacy emissions and sinks from earlier land-use events.
- Product pools, soil pools, and regrowth trajectories require an adequate initialization period.
- Urban land is not yet implemented as a complete active transition class.
- Gross source and sink diagnostics are internal bookkeeping outputs and are not always directly comparable with published process-level net components.

---

## 17. Development workflow

The default branch is the common base for model development.

```bash
git switch main
git pull --ff-only
git switch -c new-feature

# edit and test code

git add <files>
git commit -m "feat: describe the change"
git push -u origin new-feature
```

Repository-wide documentation should be maintained in this single `README.md`. Branch-specific behavior is documented in the branch table and harvest sections above rather than in separate README files.

Future tasks and planned changes can be recorded in:

```text
ROADMAP.md
```

Major completed changes can be recorded in:

```text
CHANGELOG.md
```

---

## 18. Suggested method description

> We used an event-based land-use and land-cover change carbon bookkeeping model driven by annual LUH2 land-cover states, transition matrices, and wood-harvest forcing. The model tracks biomass, slow- and rapid-soil carbon, wood-product pools, and atmospheric carbon for each grid cell while retaining legacy effects from historical clearing, abandonment, harvest, and agricultural land conversion. Carbon transfers were calculated using PFT-dependent carbon densities and event-specific allocation and response parameters. Annual gridded outputs were summarized as gross sources, gross sinks, net LULCC emissions, process-level fluxes, and LUCE-style transition categories.

---

## 19. Repository maintenance note

After merging this content into `README.md`, the separate file:

```text
README_area_driven_harvest.md
```

is no longer required and should be removed from the repository. Its content is incorporated into Sections 1.1 and 13 of this README.
