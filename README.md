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
├── server/
│   ├── __init__.py
│   ├── main_1deg.py
│   ├── main_025deg.py
│   ├── run_single_band.sbatch
│   └── clean_run.sh
│
└── docs/
    ├── CHANGELOG.txt
    └── ROADMAP.md
```

The repository uses two model branches for different land-use input formats:

| Branch | Land-use input | Main difference |
|---|---|---|
| `main` | Original LUH2 state and transition variables | Original LUH2 classes are aggregated inside the model |
| `vscp` | Pre-aggregated `v/s/c/p/U` state and transition variables | Aggregated variables are read directly |

Harvest schemes are not separated by Git branch. They are selected through experiment YAML files.

---

## 2. Input-format branches

### 2.1 `main`: original LUH2 input

The `main` branch reads original LUH2 state variables:

```text
primf, primn, secdf, secdn, urban,
c3ann, c4ann, c3per, c4per, c3nfx,
pastr, range
```

The model aggregates these variables internally:

```text
primf + primn                          → v
secdf + secdn                          → s
pastr + range                          → p
c3ann + c4ann + c3per + c4per + c3nfx → c
urban                                  → U
```

Original LUH2 state-to-state transition channels are also aggregated inside `LULCCSimulator`.

### 2.2 `vscp`: pre-aggregated input

The `vscp` branch reads the following state variables directly:

```text
v
s
c
p
U
```

The transition file contains variables named as:

```text
<src>_to_<dst>_frac
```

Examples:

```text
v_to_c_frac
v_to_p_frac
s_to_c_frac
s_to_p_frac
c_to_s_frac
p_to_s_frac
c_to_p_frac
p_to_c_frac
```

Same-class channels such as `v_to_v_frac` or `s_to_s_frac` may be retained for transition-matrix checks, but they are not treated as LULCC events by the simulator.

The aggregated transition file must still retain the original LUH2 wood-harvest variables:

```text
primf_harv     primf_bioh
primn_harv     primn_bioh
secmf_harv     secmf_bioh
secyf_harv     secyf_bioh
secnf_harv     secnf_bioh
```

Ordinary land-use transitions are therefore read as `v/s/c/p/U`, while wood harvest continues to use the original LUH2 harvest families.

---

## 3. Harvest experiment configurations

Three harvest experiments are available.

| Experiment alias | Configuration file | Harvest forcing | Area treatment |
|---|---|---|---|
| `area` | `config/experiments/harvest_area.yml` | LUH2 `*_harv` area | Strictly limited to the LUH2 harvest footprint |
| `bio-strict` | `config/experiments/harvest_bio_strict.yml` | LUH2 `*_bioh` biomass demand | Demand can only be met within the LUH2 harvest footprint |
| `bio-forced` | `config/experiments/harvest_bio_forced.yml` | LUH2 `*_bioh` biomass demand | Harvest may expand beyond the LUH2 footprint within the same source cover and allowed PFT group |

### 3.1 Area-driven harvest

The area-driven experiment uses LUH2 `*_harv` variables.

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

### 3.2 Strict-area biomass-demand harvest

The strict biomass-demand experiment uses LUH2 `*_bioh` as the requested biomass removal. Harvest is limited to the LUH2 `*_harv` footprint.

```text
requested biomass = met biomass + unmet biomass
removed biomass   = met biomass
forced biomass    = 0
extra area        = 0
```

Unmet demand is reported but is not added to product or rapid-soil pools.

### 3.3 Biomass-demand harvest with area expansion

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

## 4. Main model modules

| File | Main role |
|---|---|
| `src/main.py` | Local command-line entry point; resolves experiment configuration and local input/output paths |
| `src/LULCCSimulator.py` | Controls data parsing, initialization, annual simulation, diagnostics, and NetCDF writing |
| `src/events.py` | Implements clearing, abandonment, and cropland–pasture transition events |
| `src/harvest.py` | Implements area-driven and biomass-demand harvest schemes |
| `src/transition.py` | Applies annual pool relaxation, decay, regrowth, and process-level flux accounting |
| `src/file_loader.py` | Loads and aligns land-use and PFT input files |
| `src/parameter_loader.py` | Reads the base configuration, overlays experiment configuration, and loads optional transient carbon-density data |
| `src/carbon_pools_init.py` | Defines model indices and initializes equilibrium and legacy carbon pools |
| `src/summary_yearly.py` | Summarizes annual stocks and atmospheric carbon |
| `server/main_1deg.py` | Runs one 1° longitude band |
| `server/main_025deg.py` | Runs one 0.25° longitude band |
| `server/run_single_band.sbatch` | Selects resolution and experiment, then launches a SLURM array task |
| `server/clean_run.sh` | Removes selected model outputs and logs |

The main call sequence is:

```text
local:  src.main
server: server.main_1deg or server.main_025deg
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

## 5. Internal land-cover classes and events

The model uses five internal categories.

| Internal code | Meaning |
|---|---|
| `v` | Primary vegetation |
| `s` | Secondary vegetation |
| `p` | Pasture and rangeland |
| `c` | Cropland |
| `U` | Urban land |

The internal carbon-array order is:

```text
v, s, p, c, U
```

Input variables should always be matched by name rather than by their order in a NetCDF file.

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

Wood harvest is processed separately from ordinary state-to-state transitions using:

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

Urban land is retained for land-cover accounting, but `U`-related transitions are not yet implemented as a complete active bookkeeping process.

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

States, transitions, and the PFT map must use mutually consistent:

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
dynamic_carbon_density:
  enabled: false
```

to use the static carbon densities stored directly in `config/config.yml`.

When transient carbon density is enabled, the Parquet file must contain year- and PFT-specific biomass and soil carbon densities consistent with the model configuration.

---

## 8. Required input data and server layout

The server project layout is:

```text
/mnt/beegfs/product/lulc0120/
├── Code/                    # Git repository
├── In_ncfile/               # all input files are stored directly here
├── Out_ncfile/
├── logs/
└── site-packages/
```

There are no secondary input-data directories under `In_ncfile`. Original LUH2 inputs, aggregated VSCP inputs, PFT maps, and files at different resolutions may coexist in the same directory.

Example:

```text
/mnt/beegfs/product/lulc0120/In_ncfile/
├── states.nc
├── transitions.nc
├── states_1deg.nc
├── transitions_1deg.nc
├── states_vscp.nc
├── transitions_vscp.nc
├── states_vscp_1deg.nc
├── transitions_vscp_1deg.nc
├── IBIS_PFT_dominant_0.25deg.nc
└── IBIS_PFT_dominant_1deg.nc
```

Relative input filenames supplied to the server runners are resolved directly under:

```text
/mnt/beegfs/product/lulc0120/In_ncfile/
```

Current server defaults are:

| Resolution | State file | Transition file | PFT file |
|---|---|---|---|
| 1° | `states_1deg.nc` | `transitions_1deg.nc` | `IBIS_PFT_dominant_1deg.nc` |
| 0.25° | `states.nc` | `transitions.nc` | `IBIS_PFT_dominant_0.25deg.nc` |

The active input files must match the current code branch:

| Branch | Required state/transition content |
|---|---|
| `main` | Original LUH2 variables |
| `vscp` | Aggregated `v/s/c/p/U` variables |

When the actual filenames differ from the defaults, `server.main_1deg` and `server.main_025deg` accept:

```text
--state-file
--transition-file
--pft-file
```

The current `run_single_band.sbatch` wrapper forwards the experiment alias but does not forward custom input filenames. Therefore, SLURM array runs use the default filenames defined in the selected server runner unless those defaults are changed in that branch.

Local input and output directories are excluded from Git through `.gitignore`.

---

## 9. Installation

Python 3.8 or later is recommended.

### Local environment

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

### Server-local packages

The SLURM script adds the following directory to `PYTHONPATH`:

```text
/mnt/beegfs/product/lulc0120/site-packages
```

Packages can be installed there with:

```bash
python3 -m pip install \
  --target /mnt/beegfs/product/lulc0120/site-packages \
  -r requirements.txt
```

Required packages include:

```text
numpy
netCDF4
PyYAML
pandas
pyarrow
```

`pandas` and `pyarrow` are required when transient carbon density is enabled.

---

## 10. Running locally

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

The local paths and simulation length are currently defined in `src/main.py`. Update them before running with a different resolution, input filename, simulation period, or output location.

---

## 11. Running on the server

Run all commands from:

```bash
cd /mnt/beegfs/product/lulc0120/Code
```

The SLURM wrapper selects the runner through:

```text
RESOLUTION=1deg   → python3 -m server.main_1deg
RESOLUTION=025deg → python3 -m server.main_025deg
```

The experiment is selected independently:

```text
EXPERIMENT=area
EXPERIMENT=bio-strict
EXPERIMENT=bio-forced
```

### 11.1 1° global run

With the default `BAND_SIZE=10`, the globe is divided into 36 longitude bands:

```text
band IDs: 0–35
```

```bash
sbatch \
  --array=0-35%25 \
  --export=ALL,RESOLUTION=1deg,EXPERIMENT=area,BAND_SIZE=10,STAGGER_MAX=600 \
  server/run_single_band.sbatch
```

Replace `EXPERIMENT=area` with `bio-strict` or `bio-forced` as needed.

### 11.2 0.25° global run

With the default `BAND_SIZE=3`, the globe is divided into 120 longitude bands:

```text
band IDs: 0–119
```

```bash
sbatch \
  --array=0-119%35 \
  --export=ALL,RESOLUTION=025deg,EXPERIMENT=area,BAND_SIZE=3,STAGGER_MAX=600 \
  server/run_single_band.sbatch
```

### 11.3 Single-band run

1° band 0:

```bash
sbatch \
  --array=0-0 \
  --export=ALL,RESOLUTION=1deg,EXPERIMENT=area,BAND_SIZE=10,STAGGER_MAX=0 \
  server/run_single_band.sbatch
```

0.25° band 0:

```bash
sbatch \
  --array=0-0 \
  --export=ALL,RESOLUTION=025deg,EXPERIMENT=area,BAND_SIZE=3,STAGGER_MAX=0 \
  server/run_single_band.sbatch
```

To rerun a specific band, set the array range to that band ID. For example:

```bash
sbatch \
  --array=100-100 \
  --export=ALL,RESOLUTION=025deg,EXPERIMENT=area,BAND_SIZE=3,STAGGER_MAX=0 \
  server/run_single_band.sbatch
```

### 11.4 Direct server test with explicit filenames

The Python runner can be called directly when custom input filenames are needed:

```bash
BAND_ID=0 BAND_SIZE=10 STAGGER_MAX=0 \
python3 -m server.main_1deg \
  --experiment area \
  --state-file states_vscp_1deg.nc \
  --transition-file transitions_vscp_1deg.nc \
  --pft-file IBIS_PFT_dominant_1deg.nc \
  --start-year 850 \
  --end-year 2020
```

All three relative input filenames in this command are resolved under the same `In_ncfile` directory.

### 11.5 Output locations

1° outputs:

```text
/mnt/beegfs/product/lulc0120/Out_ncfile/
└── 1deg/
    └── <experiment_name>/
        ├── summary_1deg.rank000.nc
        ├── summary_1deg.rank001.nc
        └── ...
```

0.25° outputs:

```text
/mnt/beegfs/product/lulc0120/Out_ncfile/
└── 025deg/
    └── <experiment_name>/
        ├── summary_025deg.rank000.nc
        ├── summary_025deg.rank001.nc
        └── ...
```

Experiment directory names are read from the YAML files:

```text
harvest_area_driven
harvest_bio_strict
harvest_bio_forced
```

A band is skipped when its output has the expected dimensions and all `done` flags equal 1. An existing incomplete band file is removed and the complete band is run again.

### 11.6 Logs

SLURM output and error logs:

```text
/mnt/beegfs/product/lulc0120/logs/
├── slurm_band_<jobID>_<taskID>.out
└── slurm_band_<jobID>_<taskID>.err
```

Model logs:

```text
lulcc_<resolution>_<experiment>_<jobID>_<taskID>.log
```

### 11.7 Cleaning model outputs

Remove one resolution/experiment combination:

```bash
bash server/clean_run.sh 1deg area
bash server/clean_run.sh 025deg bio-strict
bash server/clean_run.sh 025deg bio-forced
```

Remove all model outputs and logs:

```bash
bash server/clean_run.sh all
```

The cleanup script asks for confirmation before deletion and does not remove files from `In_ncfile`.

---

## 12. Simulation length and spatial bands

`years` is the number of annual transitions simulated:

```text
number of stored states = years + 1
```

The server runners calculate:

```text
total_years = end_year - start_year
expected_time = total_years + 1
```

For example, a simulation from 850 through 2020 contains:

```text
1170 annual transitions
1171 stored states
```

`start_year_idx` is an index offset into the external forcing data. It is not added to the internal carbon-pool array length.

The server model is divided into longitude bands. The number of bands is:

```text
bands_total = 360 / BAND_SIZE
```

`BAND_SIZE` must divide 360 exactly.

---

## 13. Output variables and units

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
Flux_FD
Flux_NFC
Flux_FR
Flux_NFR
Flux_CAL
Flux_WHp
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

### Harvest demand and removal

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

### Harvest area

```text
harvest_luh2_area_frac
harvest_luh2_area
harvest_effective_area_frac
harvest_effective_area
harvest_extra_area
area_harvest
```

### Harvest allocation and equation diagnostics

```text
harvest_to_products
harvest_to_soil
harvest_to_SR_from_biomass
harvest_to_SR_from_soil
harvest_forced_to_products
harvest_forced_to_soil
harvest_beta_delta_sum
harvest_sigma_delta_sum
harvest_delta_B_h
harvest_delta_SS_h
harvest_R
```

Variables that are not applicable to the selected experiment are retained as zero so that the output schema remains consistent across experiments.

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

---

## 14. Known cautions

- The code branch and the land-use input format must be consistent.
- The `main` branch expects original LUH2 state and transition variables.
- The `vscp` branch expects pre-aggregated `v/s/c/p/U` variables.
- Aggregated transition files must retain the original LUH2 `*_harv` and `*_bioh` variables.
- Carbon-density assumptions strongly affect the magnitude of LULCC emissions.
- Static and transient carbon-density simulations should be clearly distinguished.
- LUH2 `*_harv` and `*_bioh` represent different constraints and should not be treated as interchangeable.
- Harvest biomass must not be added to product or soil pools unless it has been removed from the modeled ecosystem stock.
- PFT and land-use grid misalignment can create large spatial artifacts.
- A reversed latitude axis can invert the spatial result.
- Short simulations omit legacy emissions and sinks from earlier land-use events.
- Product pools, soil pools, and regrowth trajectories require an adequate initialization period.
- Urban land is not yet implemented as a complete active transition class.
- Gross source and sink diagnostics are internal bookkeeping outputs and are not always directly comparable with published process-level net components.
- The `done` variable records completion. The current server runner restarts an incomplete band from the beginning rather than continuing within the same band.
