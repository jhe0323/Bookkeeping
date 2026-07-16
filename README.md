# LULCC Carbon Bookkeeping Model

A Python-based land-use and land-cover change (LULCC) carbon bookkeeping model for estimating annual carbon-stock changes and land-use change emissions.

The model is driven by:

- LUH2 land-cover states and transitions;
- LUH2 wood-harvest area and biomass-demand variables;
- PFT-dependent biomass and soil carbon densities;
- event-specific carbon-allocation and response parameters.

The model tracks carbon changes associated with:

- land clearing;
- agricultural abandonment and regrowth;
- wood harvest;
- cropland–pasture conversion;
- product-pool decay;
- biomass and soil-carbon relaxation.

The code is designed for comparison with BLUE-, LUCE-, and Global Carbon Budget-style bookkeeping estimates.

---

## 1. Repository structure

```text
Bookkeeping/
├── README.md
├── requirements.txt
├── .gitignore
│
├── src/
│   ├── __init__.py
│   ├── main.py
│   ├── LULCCSimulator.py
│   ├── events.py
│   ├── harvest.py
│   ├── transition.py
│   ├── file_loader.py
│   ├── parameter_loader.py
│   ├── carbon_pools_init.py
│   └── summary_yearly.py
│
├── config/
│   ├── config.yml
│   ├── dynamic_carbon_density.parquet
│   └── experiments/
│       ├── harvest_area.yml
│       ├── harvest_bio_strict.yml
│       └── harvest_bio_forced.yml
│
└── docs/
    ├── CHANGELOG.txt
    └── ROADMAP.md
```

The repository now uses a **single codebase on `main`**. Alternative harvest schemes are controlled by experiment configuration files rather than separate Git branches.

---

## 2. Harvest experiment configurations

Three harvest experiments are available.

| Experiment alias | Configuration file | Harvest forcing | Area treatment |
|---|---|---|---|
| `area` | `config/experiments/harvest_area.yml` | LUH2 `*_harv` area | Strictly limited to the LUH2 harvest footprint |
| `bio-strict` | `config/experiments/harvest_bio_strict.yml` | LUH2 `*_bioh` biomass demand | Demand can only be met within the LUH2 harvest footprint |
| `bio-forced` | `config/experiments/harvest_bio_forced.yml` | LUH2 `*_bioh` biomass demand | Harvest may expand beyond the LUH2 footprint within the same source cover and allowed PFT group |

### 2.1 Area-driven harvest

The area-driven experiment uses only LUH2 `*_harv` variables.

```text
removed biomass
= current modeled biomass on the source cover
× harvested share of the source-cover area
```

Expected diagnostic identities:

```text
harvest_requested_biomass = harvest_met_biomass
harvest_biomass_removed   = harvest_met_biomass
harvest_forced_biomass    = 0
harvest_unmet_biomass     = 0
harvest_extra_area        = 0
```

### 2.2 Strict-area biomass-demand harvest

The strict biomass-demand experiment uses LUH2 `*_bioh` as the requested biomass removal. Harvest is limited to the LUH2 `*_harv` footprint.

```text
requested biomass = met biomass + unmet biomass
removed biomass   = met biomass
forced biomass    = 0
extra area        = 0
```

Unmet demand is reported but is not added to product or rapid-soil pools.

### 2.3 Biomass-demand harvest with area expansion

The forced biomass-demand experiment first uses biomass within the LUH2 harvest footprint. If this is insufficient, harvest may expand within the same source cover, harvest family, and allowed PFT group.

```text
requested biomass
= met biomass + forced biomass + unmet biomass

removed biomass
= met biomass + forced biomass
```

Additional removal outside the LUH2 footprint is reported as:

```text
harvest_forced_biomass
harvest_extra_area
```

---

## 3. Main model modules

| File | Main role |
|---|---|
| `src/main.py` | Command-line entry point; resolves experiment configuration, input paths, output path, simulation period, and spatial domain |
| `src/LULCCSimulator.py` | Controls loading, initialization, yearly simulation, diagnostics, and NetCDF writing |
| `src/events.py` | Implements clearing, abandonment, and cropland–pasture transition events |
| `src/harvest.py` | Implements all area-driven and biomass-demand harvest schemes |
| `src/transition.py` | Applies annual pool relaxation, decay, regrowth, and process-level flux accounting |
| `src/file_loader.py` | Loads and aligns LUH2 and PFT input files |
| `src/parameter_loader.py` | Reads the base configuration, overlays experiment configuration, and loads optional transient carbon-density data |
| `src/carbon_pools_init.py` | Defines model indices and initializes equilibrium and legacy carbon pools |
| `src/summary_yearly.py` | Summarizes annual stocks and atmospheric carbon |
| `config/config.yml` | Stores land-cover mapping, carbon densities, allocation coefficients, and response times |
| `config/experiments/*.yml` | Selects the active harvest experiment |

The main call sequence is:

```text
src.main
  ↓
resolve experiment YAML
  ↓
LULCCSimulator
  ↓
ParameterLoader
  ├── config/config.yml
  └── config/experiments/<experiment>.yml
  ↓
FileLoader
  ├── states
  ├── transitions
  └── PFT map
  ↓
initialize C_bar, Delta, and frac_area
  ↓
annual simulation
  ├── clearing
  ├── abandonment
  ├── cropland–pasture conversion
  ├── wood harvest
  ├── pool relaxation and decay
  └── annual summary
  ↓
NetCDF output
```

---

## 4. Internal land-cover classes

The model aggregates LUH2 land-use classes into five internal categories.

| Internal code | Meaning | Typical LUH2 classes |
|---|---|---|
| `v` | Primary vegetation | `primf`, `primn` |
| `s` | Secondary vegetation | `secdf`, `secdn` |
| `p` | Pasture and rangeland | `pastr`, `range` |
| `c` | Cropland | `c3ann`, `c4ann`, `c3per`, `c4per`, `c3nfx` |
| `U` | Urban land | `urban` |

The active bookkeeping calculations mainly use `v`, `s`, `p`, and `c`. Urban land is retained for land-cover accounting but is not currently implemented as a complete active event class.

---

## 5. Event classification

### Clearing

```text
v → c
v → p
s → c
s → p
```

Clearing transfers carbon among direct atmospheric release, product pools, rapid soil/slash, slow soil, and biomass-deficit pools.

### Abandonment and regrowth

```text
c → s
p → s
```

These transitions inherit the source carbon state into secondary vegetation and generate biomass and soil recovery trajectories.

### Wood harvest

Wood harvest is processed separately from ordinary state-to-state transitions using the following LUH2 harvest families:

```text
primf
primn
secmf
secyf
secnf
```

Primary harvest is represented as a transfer from primary to secondary vegetation. Secondary harvest remains within secondary vegetation but creates a new harvest-history carbon deficit.

### Other transitions

```text
c → p
p → c
```

These transitions are handled separately from clearing and abandonment.

---

## 6. Carbon pools

The equilibrium carbon array stores:

```text
B   biomass carbon
SS  slow soil carbon
```

The legacy or excess carbon array stores:

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

## 7. PFT and carbon-density configuration

The current model uses 15 PFT classes. PFT-dependent parameters are stored in `config/config.yml`.

The PFT input may be:

- a 15-class PFT fraction map; or
- a dominant-PFT map converted to the expected grid and class convention.

All input files should use mutually consistent:

- spatial resolution;
- grid-cell centers;
- latitude orientation;
- longitude convention;
- land mask;
- PFT numbering.

### Static and transient carbon density

The model supports static and transient PFT-dependent biomass and soil carbon densities.

```yaml
dynamic_carbon_density:
  enabled: true
  path: dynamic_carbon_density.parquet
  base_year: 850
  min_year_policy: clip
  max_year_policy: clip
```

Set:

```yaml
enabled: false
```

to use the static carbon densities stored directly in `config/config.yml`.

When transient density is enabled, the Parquet file must contain year- and PFT-specific biomass and soil carbon densities for `v`, `s`, `p`, and `c`.

---

## 8. Required input data

The default command-line entry point expects:

```text
Bookkeeping/
├── In_ncfile/
│   ├── states_1deg.nc
│   ├── transitions_1deg.nc
│   └── PFTmap_orchidee_1deg.nc
│
├── Out_ncfile/
│
├── config/
└── src/
```

The paths are currently defined near the top of `src/main.py`:

```python
STATE_PATH = DATA_DIR / "states_1deg.nc"
TRANS_PATH = DATA_DIR / "transitions_1deg.nc"
PFT_PATH = DATA_DIR / "PFTmap_orchidee_1deg.nc"
CONFIG_PATH = CONFIG_DIR / "config.yml"
```

Update these values when using different input filenames.

The input data directories and generated outputs are intentionally excluded from Git through `.gitignore`.

---

## 9. Installation

Python 3.8 or later is recommended.

Create and activate a virtual environment, then install dependencies:

```bash
python -m venv .venv
```

Windows:

```bash
.venv\Scripts\activate
```

Linux or macOS:

```bash
source .venv/bin/activate
```

Install packages:

```bash
pip install -r requirements.txt
```

Required packages include:

```text
numpy
netCDF4
PyYAML
pandas
pyarrow
```

`pandas` and `pyarrow` are required when transient carbon density is enabled and the Parquet file is read.

---

## 10. Running the model

Run the model from the repository root using module mode.

### Area-driven experiment

```bash
python -m src.main --experiment area
```

### Strict-area biomass-demand experiment

```bash
python -m src.main --experiment bio-strict
```

### Biomass-demand experiment with area expansion

```bash
python -m src.main --experiment bio-forced
```

A custom experiment YAML can also be supplied:

```bash
python -m src.main \
  --experiment config/experiments/harvest_bio_forced.yml
```

Output files are named from the experiment configuration, for example:

```text
Out_ncfile/summary_global_harvest_area_driven.nc
Out_ncfile/summary_global_harvest_bio_strict.nc
Out_ncfile/summary_global_harvest_bio_forced.nc
```

The NetCDF file also stores experiment metadata:

```text
experiment_name
harvest_mode
harvest_use_bioh
harvest_allow_expansion
```

---

## 11. Simulation length and time indexing

`years` is the number of annual transitions simulated.

```text
number of stored states = years + 1
```

For example, states from 850 through 2022 contain:

```text
1173 stored states
1172 annual transitions
```

Always verify the actual time dimensions of `states.nc` and `transitions.nc` before a production run.

`start_year_idx` is an index offset into the external forcing data. It is not added to the internal carbon-pool array length.

---

## 12. Spatial subsets and band runs

`lat_slice` and `lon_slice` can be passed to `LULCCSimulator` to load and simulate only part of the global grid.

```python
simulator = LULCCSimulator(
    config_path="config/config.yml",
    experiment_path="config/experiments/harvest_area.yml",
    LULC_path="In_ncfile/states_1deg.nc",
    trans_path="In_ncfile/transitions_1deg.nc",
    pft_path="In_ncfile/PFTmap_orchidee_1deg.nc",
    lat_slice=slice(80, 100),
    lon_slice=slice(120, 140),
    area_unit="ha",
)
```

For high-resolution global simulations, the domain can be divided into longitude or latitude bands. Each band should write to a separate NetCDF file and the outputs can be merged after all bands finish.

---

## 13. Output variables

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

### LUCE-style process categories

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

## 14. Harvest diagnostics

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

`harvest_forced_to_products` and `harvest_forced_to_soil` are retained for backward-compatible output schemas. Forced harvest is already included in the physical totals `harvest_to_products` and `harvest_to_soil`.

### Equation-level diagnostics

```text
harvest_beta_delta_sum
harvest_sigma_delta_sum
harvest_delta_B_h
harvest_delta_SS_h
harvest_R
```

Variables that are not applicable to the selected experiment are retained and reported as zero so that outputs from different experiments can be compared using the same analysis scripts.

---

## 15. Units

When:

```text
area_unit = ha
carbon density = t C ha^-1
```

grid-cell carbon stocks and annual fluxes are stored in:

```text
t C
```

Common conversions:

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

## 16. Recommended validation workflow

Before running a full global simulation:

1. Check that states, transitions, and PFT maps use the same grid.
2. Verify latitude orientation and longitude convention.
3. Confirm that PFT fractions sum to 1 over valid land cells.
4. Confirm that `years` and `start_year_idx` remain within the forcing time dimension.
5. Run a single-grid or small-region test.
6. Test clearing, abandonment, harvest, and crop–pasture conversion separately.
7. Check biomass, soil, product, atmosphere, and total-system carbon closure.
8. Verify that no output variables are unexpectedly all zero or all `NaN`.
9. Compare annual gross sources, gross sinks, and net emissions.
10. Inspect harvest requested, met, forced, unmet, and area diagnostics.
11. Inspect spatial maps for coordinate reversal or unrealistic hotspots.
12. Run the full global simulation only after these checks pass.

Recommended harvest identities:

```text
area:
requested = met = removed
forced = unmet = 0
```

```text
bio-strict:
requested = met + unmet
removed = met
forced = 0
```

```text
bio-forced:
requested = met + forced + unmet
removed = met + forced
```

---

## 17. Known cautions

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
- The NetCDF `done` variable records completion during the current run; the current writer does not yet reopen an interrupted output file for automatic restart.

---

## 18. Development workflow

Use `main` as the common model codebase. Create temporary feature branches only for code development, not for experiment definitions.

```bash
git switch main
git pull --ff-only
git switch -c feature/<name>

# edit and test
git add <files>
git commit -m "feat: describe the change"
git push -u origin feature/<name>
```

After review and testing, merge the feature branch into `main`.

Alternative model experiments should normally be added as configuration files under:

```text
config/experiments/
```

Future tasks and planned changes are recorded in:

```text
docs/ROADMAP.md
```

Completed model changes are recorded in:

```text
docs/CHANGELOG.txt
```
