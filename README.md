# LULCC Carbon Bookkeeping Model

A Python-based land-use and land-cover change (LULCC) carbon bookkeeping model for estimating annual carbon-stock changes, gross carbon sources and sinks, and net land-use change emissions.

The model follows the general bookkeeping framework used by BLUE- and LUCE-type models. Land-use events create departures from equilibrium biomass and soil carbon stocks. These departures are retained in legacy pools and relax over time, allowing historical land-use change to continue affecting present-day carbon fluxes.

The model is driven by:

- land-cover state data;
- annual land-use transition data;
- wood-harvest area or biomass-demand data;
- plant functional type (PFT) maps;
- PFT- and land-cover-specific biomass and soil carbon densities;
- event-specific allocation and response parameters.

The current `main` branch supports both original LUH2 variables and pre-aggregated `v/s/c/p/U` inputs through run configuration files. The same deterministic bookkeeping core is also used by the Monte Carlo uncertainty workflow.

---

## 1. Model framework

The model separates carbon storage into an equilibrium component and a legacy component.

### 1.1 Equilibrium carbon

The equilibrium carbon array, `C_bar`, contains the carbon stock expected for the current combination of:

- land-cover class;
- PFT;
- biomass or slow-soil pool;
- carbon-density parameter.

The equilibrium pools are:

```text
B   equilibrium biomass carbon
SS  equilibrium slow-soil carbon
```

### 1.2 Legacy carbon

Land-use events create differences between current carbon stocks and the equilibrium state. These differences are stored in `Delta`.

```text
B     biomass departure from equilibrium
SS    slow-soil departure from equilibrium
SR    rapid soil and slash carbon
P1    1-year product pool
P10   10-year product pool
P100  100-year product pool
A     direct atmospheric release
```

Legacy carbon is tracked separately by event history:

```text
l   clearing
h   harvest
a   abandonment
g   other land-use transitions
```

The annual model step relaxes these departures toward equilibrium and transfers released carbon to the atmospheric pool.

### 1.3 Annual carbon budget

The main atmospheric budget variables follow:

```text
Net_Emissions = Gross_Sources + Gross_Sinks
```

where:

```text
positive value   atmospheric carbon source
negative value   atmospheric carbon sink
```

The system-level carbon total is:

```text
system_carbon_total
= biomass_total
+ soil_total
+ P1
+ P10
+ P100
+ atmosphere
```

For static carbon density and conservative PFT remapping, the system closure error should remain close to numerical zero.

---

## 2. Internal land-cover classes

The model uses five internal land-cover classes.

| Code | Meaning |
|---|---|
| `v` | Primary vegetation |
| `s` | Secondary vegetation |
| `p` | Pasture and rangeland |
| `c` | Cropland |
| `U` | Urban land |

The array order is:

```text
v, s, p, c, U
```

The mapping from original LUH2 classes is controlled by `LUH2toLULC` in `config/config.yml`.

The current default mapping is:

```text
primf + primn                           → v
secdf + secdn                           → s
pastr + range                           → p
c3ann + c4ann + c3per + c4per + c3nfx → c
urban                                   → U
```

The mapping is configurable and should remain consistent with the PFT map and carbon-density parameterization.

---

## 3. Land-use processes

### 3.1 Clearing

Clearing includes:

```text
v → c
v → p
s → c
s → p
```

Clearing removes biomass from the source land-cover class and partitions it among:

- direct atmospheric release;
- 1-, 10-, and 100-year product pools;
- rapid soil and slash carbon;
- slow-soil carbon.

The target agricultural biomass stock is represented as an immediate biomass deficit and subsequently relaxes toward equilibrium.

### 3.2 Abandonment and regrowth

Abandonment includes:

```text
c → s
p → s
```

The source biomass and slow-soil state are inherited by secondary vegetation. The difference between the inherited state and the secondary-vegetation equilibrium state creates biomass and soil recovery trajectories.

### 3.3 Cropland-pasture conversion

Other agricultural transitions include:

```text
c → p
p → c
```

These transitions retain the source carbon state and create departures from the target equilibrium state.

### 3.4 Wood harvest

Wood harvest is processed separately from ordinary state-to-state transitions.

Supported LUH2 harvest families are:

```text
primf
primn
secmf
secyf
secnf
```

Primary harvest transfers affected area from primary vegetation to secondary vegetation. Secondary harvest remains within secondary vegetation but creates a new harvest-history biomass and soil departure.

### 3.5 Urban land

Urban land is retained in state accounting. Carbon-active transitions involving `U` are not yet implemented as a complete bookkeeping process.

---

## 4. Beginning-of-year transition accounting

All ordinary transitions in one model year are evaluated against the same beginning-of-year land-cover and carbon state.

For each source land-cover class:

1. all requested outgoing transitions are summed;
2. requested area is compared with beginning-of-year source area;
3. if requested area exceeds available area, all outgoing transitions are scaled proportionally;
4. individual transition effects are calculated from independent copies of the beginning-of-year state;
5. all event increments are combined and applied simultaneously.

This prevents reciprocal or multi-destination transitions from using area or carbon created earlier in the same annual loop.

The associated diagnostics are:

```text
transition_requested_area
transition_applied_area
transition_clipped_area
```

---

## 5. Land-use input formats

The `main` branch supports two input formats selected in the run YAML.

```yaml
run:
  input_format: original_luh2
```

or:

```yaml
run:
  input_format: vscp
```

### 5.1 Original LUH2 input

For `original_luh2`, the state file contains the original LUH2 classes configured under `LUH2toLULC`.

Typical state variables are:

```text
primf
primn
secdf
secdn
urban
c3ann
c4ann
c3per
c4per
c3nfx
pastr
range
```

Transition variables are detected from variables named:

```text
<src>_to_<dst>
```

or:

```text
<src>_to_<dst>_frac
```

Only variables that exist in the transition file and map to different internal land-cover classes are loaded and evaluated.

### 5.2 Pre-aggregated VSCP input

For `vscp`, the state file is read directly from:

```text
v
s
p
c
U
```

The required active state variables are:

```text
v
s
p
c
```

`U` is optional.

Transition variables are named:

```text
<src>_to_<dst>_frac
```

Examples include:

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

A VSCP transition file must retain the original LUH2 harvest variables when wood harvest is simulated.

---

## 6. Wood-harvest experiments

Three experiment configurations are provided.

| Configuration file | Mode | Demand source | Area behavior |
|---|---|---|---|
| `config/experiments/harvest_area.yml` | Area driven | LUH2 `*_harv` | Limited to the reported harvest footprint |
| `config/experiments/harvest_bio_strict.yml` | Biomass demand | LUH2 `*_bioh` | Demand can only be met within the reported footprint |
| `config/experiments/harvest_bio_forced.yml` | Biomass demand | LUH2 `*_bioh` | Additional source area may be used when the footprint is insufficient |

The experiment is selected in the run YAML:

```yaml
experiment:
  path: config/experiments/harvest_area.yml
```

### 6.1 Area-driven harvest

Area-driven harvest uses the LUH2 harvested-area fraction.

```text
requested biomass
= current modeled biomass within the effective harvested area
```

Expected relationships are:

```text
harvest_requested_biomass = harvest_met_biomass
harvest_biomass_removed   = harvest_met_biomass
harvest_forced_biomass    = 0
harvest_unmet_biomass     = 0
harvest_extra_area        = 0
```

### 6.2 Strict biomass-demand harvest

The requested biomass is read from `*_bioh`. Removal is restricted to biomass available within the LUH2 harvest footprint.

```text
requested biomass = met biomass + unmet biomass
removed biomass   = met biomass
forced biomass    = 0
extra area        = 0
```

Unmet biomass is reported but is not added to product or soil pools.

### 6.3 Biomass-demand harvest with area expansion

The model first uses biomass within the LUH2 harvest footprint. Remaining demand may be met from additional area within the same source-cover and forest/non-forest PFT group.

```text
requested biomass
= met biomass
+ forced biomass
+ unmet biomass

removed biomass
= met biomass
+ forced biomass
```

Additional removal is reported through:

```text
harvest_forced_biomass
harvest_extra_area
```

---

## 7. PFT maps

PFT names and the PFT count are determined from the parameter configuration.

The preferred explicit form is:

```yaml
pft:
  names:
    - PFT1
    - PFT2
    - PFT3
  forest_indices: [1, 2]
```

When the explicit `pft` section is absent, the model infers ordered `PFT<number>` entries from `carbon_density` and reads the legacy `forest_pft_indices` field.

The PFT map can be:

```text
lat × lon
```

Static dominant PFT codes.

```text
pft × lat × lon
```

Static fractional PFT composition.

```text
time × lat × lon
```

Dynamic dominant PFT codes.

```text
time × pft × lat × lon
```

Dynamic fractional PFT composition.

Dimension order is detected from dimension names.

### 7.1 Static PFT

A static PFT map is used for the complete simulation.

### 7.2 Dynamic PFT

Dynamic PFT behavior is controlled in the run YAML.

```yaml
pft:
  base_year: 850
  time_encoding: auto
  update_mode: annual_conservative
  min_year_policy: clip
  max_year_policy: clip
```

Available update modes are:

```text
fixed_initial
```

Use the PFT composition from the simulation start year for all years.

```text
annual_conservative
```

Update PFT composition annually while conserving total ecosystem carbon during remapping. Biomass and slow-soil equilibrium stocks are recalculated, and the corresponding mismatch is retained in the general legacy pool.

`PFT_Remap_Adjustment` reports numerical carbon change during remapping and should remain close to zero.

A PFT map finer than the model grid must be area-weighted to the model resolution before simulation. The loader does not silently aggregate a finer PFT map.

### 7.3 Missing PFT cells

Missing-PFT handling is controlled in the run YAML.

```yaml
pft:
  missing_cell_policy: skip
  nearest_search_radius: 8
  default_index: null
```

Supported policies are:

```text
error
nearest
default
skip
```

The current production 1° and 0.25° configurations use:

```text
missing_cell_policy: skip
```

---

## 8. Carbon density

The model supports static and transient PFT-dependent biomass and soil carbon densities.

### 8.1 Static carbon density

Static densities are read from:

```text
config/config.yml
```

Each PFT contains biomass and soil density for:

```text
v
s
p
c
```

With:

```text
area_unit = ha
carbon density = t C ha^-1
```

grid-cell carbon stocks and annual carbon fluxes are stored in:

```text
t C
```

### 8.2 Transient carbon density

Transient density is controlled by the effective parameter configuration.

```yaml
dynamic_carbon_density:
  enabled: true
  path: dynamic_carbon_density.parquet
  base_year: 850
  min_year_policy: clip
  max_year_policy: clip
```

The run YAML can override the base parameter file:

```yaml
model_overrides:
  dynamic_carbon_density:
    enabled: false
```

Configuration precedence is:

```text
config/config.yml
  ↓
experiment YAML
  ↓
model_overrides in the run YAML
  ↓
Monte Carlo parameter overrides, when MC mode is used
```

The current transient-density implementation refreshes equilibrium biomass and soil carbon for the existing land-cover-by-PFT area each year.

The resulting external stock adjustment is reported as:

```text
Carbon_Density_Adjustment
```

This adjustment is excluded from `Closure_Error` and is not included in `Net_Emissions`.

The first production Monte Carlo experiment uses static carbon density. Transient carbon density should be treated as a separate uncertainty experiment.

---

## 9. Time configuration

Simulation dates are expressed as calendar years.

```yaml
run:
  input_base_year: 850
  start_year: 850
  end_year: 2020

inputs:
  time_encoding: index
```

For an index-based LUH2 time axis:

```text
input index 0    → calendar year 850
input index 1170 → calendar year 2020
```

The number of annual transitions is:

```text
end_year - start_year
```

The number of stored annual states is:

```text
end_year - start_year + 1
```

The transition file must contain every annual transition year from `start_year` through `end_year - 1`.

Supported time encodings are:

```text
index
calendar_year
auto
```

---

## 10. Run configuration

A complete deterministic simulation is defined by a run YAML.

Current production server configurations are:

```text
config/run_1deg.yml
config/run_025deg.yml
```

Monte Carlo-specific run configurations are:

```text
config/run_mc_1deg.yml
config/run_mc_025deg.yml
```

A run configuration contains:

```text
run             simulation identity, resolution, input format and years
paths           input and output directories
inputs          state, transition and PFT files
parameters      base parameter file
experiment      harvest experiment file
model_overrides run-specific parameter switches
pft             PFT time handling and missing-cell behavior
validation      input validation settings
output          output prefix, compression and synchronization
server          longitude-band width and start staggering
monte_carlo     MC sample table and summary-output configuration, MC runs only
```

### 10.1 Current 1° deterministic configuration

The current server baseline uses pre-aggregated VSCP inputs:

```yaml
run:
  name: baseline_1deg_voidpft_v1
  resolution: 1deg
  input_format: vscp
  input_base_year: 850
  start_year: 850
  end_year: 2020
  area_unit: ha

paths:
  data_dir: /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/In_ncfile
  output_dir: /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/Out_ncfile

inputs:
  states: synthetic_states_1deg.nc
  transitions: synthetic_transitions_1deg.nc
  pft: IBIS_PFT_dominant_1deg.nc
  pft_variable: null
  time_encoding: index

parameters:
  path: config/config.yml

experiment:
  path: config/experiments/harvest_area.yml

model_overrides:
  dynamic_carbon_density:
    enabled: false

pft:
  base_year: 850
  time_encoding: auto
  update_mode: annual_conservative
  min_year_policy: clip
  max_year_policy: clip
  missing_cell_policy: skip
  nearest_search_radius: 8
  default_index: null

validation:
  run_before_simulation: false
  full_scan: true

output:
  filename_prefix: summary_1deg
  compression_level: 4
  sync_every: 200

server:
  band_size_deg: 10
  stagger_max: 0
```

### 10.2 Current 0.25° deterministic configuration

The corresponding production setup uses:

```text
resolution: 025deg
states: synthetic_states_025deg.nc
transitions: synthetic_transitions_025deg.nc
pft: IBIS_PFT_dominant_0.25deg.nc
missing_cell_policy: skip
band_size_deg: 3
```

Use a unique `run.name` whenever any input, experiment, model switch, parameter file, or scientifically relevant setup changes.

---

## 11. Input validation

The validator can check:

- required state variables;
- required transition and harvest variables;
- simulation-year coverage;
- state and transition grid consistency;
- negative or out-of-range fractions;
- land-cover state sums;
- PFT variable identification;
- dominant-PFT integer codes;
- fractional-PFT sums;
- configured PFT count;
- dynamic-PFT year coverage.

### 11.1 Validation before simulation

Enable automatic validation in a run YAML:

```yaml
validation:
  run_before_simulation: true
  full_scan: true
```

A validation report is written to:

```text
Out_ncfile/<resolution>/<run.name>/input_validation_report.json
```

The simulation stops when validation status is `FAIL`.

For server band jobs, `validation.run_before_simulation` should normally remain `false` so every band does not rescan the same global files.

### 11.2 Standalone validation

The current `tools/validate_inputs.py` uses a `CONFIG_PATH` constant.

Edit:

```text
tools/validate_inputs.py
```

and set, for example:

```python
CONFIG_PATH = Path("config/run_1deg.yml")
```

or:

```python
CONFIG_PATH = Path("config/run_mc_1deg.yml")
```

Then run:

```bash
python -m tools.validate_inputs
```

Use `full_scan: false` for a faster structural check and `full_scan: true` for a complete data scan.

---

## 12. Repository structure

The main production-relevant structure is:

```text
Bookkeeping/
├── README.md
├── requirements.txt
├── .gitignore
│
├── src/
│   ├── __init__.py
│   ├── main.py
│   ├── run_config.py
│   ├── run_manager.py
│   ├── input_validator.py
│   ├── LULCCSimulator.py
│   ├── file_loader.py
│   ├── parameter_loader.py
│   ├── carbon_pools_init.py
│   ├── events.py
│   ├── harvest.py
│   ├── transition.py
│   ├── summary_yearly.py
│   ├── mc_sampling.py
│   └── mc_runner.py
│
├── config/
│   ├── config.yml
│   ├── dynamic_carbon_density.parquet
│   ├── run_1deg.yml
│   ├── run_025deg.yml
│   ├── run_mc_1deg.yml
│   ├── run_mc_025deg.yml
│   ├── mc_parameters_1deg.yml
│   ├── mc_parameters_025deg.yml
│   ├── mc_parameters_smoke_1deg.yml
│   ├── mc_parameters_smoke_025deg.yml
│   └── experiments/
│       ├── harvest_area.yml
│       ├── harvest_bio_strict.yml
│       └── harvest_bio_forced.yml
│
├── server/
│   ├── main_1deg.py
│   ├── main_025deg.py
│   ├── main_mc_1deg.py
│   ├── main_mc_025deg.py
│   ├── run_1deg.sh
│   ├── run_025deg.sh
│   ├── submit_1deg.sh
│   ├── submit_025deg.sh
│   ├── run_mc_1deg.sh
│   ├── run_mc_025deg.sh
│   ├── submit_mc_1deg.sh
│   ├── submit_mc_025deg.sh
│   ├── submit_mc_batch.sh
│   ├── merge_1deg.sh
│   ├── merge_025deg_local.sh
│   └── clean_run.sh
│
└── tools/
    ├── validate_inputs.py
    ├── check_pft_land_overlap.py
    ├── generate_mc_samples.py
    ├── validate_mc_band.py
    ├── merge_mc_bands.py
    ├── merge_all_mc.py
    ├── summarize_mc.py
    └── check_mc_convergence.py
```

Legacy SLURM `.sbatch` files, if retained in the repository, are not part of the current DSUB production workflow.

### Main modules

| File | Role |
|---|---|
| `src/main.py` | Local deterministic entry point |
| `src/run_config.py` | Loads and resolves run YAML files |
| `src/run_manager.py` | Shared deterministic local/server orchestration, output checks and manifests |
| `src/input_validator.py` | One-time input and configuration validation |
| `src/LULCCSimulator.py` | Annual bookkeeping simulation |
| `src/file_loader.py` | NetCDF slicing, time conversion, grid alignment and PFT loading |
| `src/parameter_loader.py` | Base parameters, experiment overlays, overrides and dynamic density |
| `src/carbon_pools_init.py` | Pool definitions and initialization |
| `src/events.py` | Clearing, abandonment and cropland-pasture events |
| `src/harvest.py` | Area- and biomass-driven harvest |
| `src/transition.py` | Annual decay, regrowth and flux attribution |
| `src/summary_yearly.py` | Annual carbon-stock summaries |
| `src/mc_sampling.py` | Reproducible MC parameter sampling and physical validation |
| `src/mc_runner.py` | Summary-only MC band execution and stale-output checks |
| `server/main_1deg.py` | 1° deterministic server entry point |
| `server/main_025deg.py` | 0.25° deterministic server entry point |
| `server/main_mc_1deg.py` | 1° MC server entry point |
| `server/main_mc_025deg.py` | 0.25° MC server entry point |

---

## 13. Installation

The current source code requires Python 3.8 or later.

Create and activate a virtual environment:

```bash
python -m venv .venv
```

Linux or macOS:

```bash
source .venv/bin/activate
```

Windows:

```text
.venv\Scripts\activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Main dependencies include:

```text
numpy
netCDF4
PyYAML
pandas
pyarrow
xarray
scipy
```

`pandas` and `pyarrow` are required by the Monte Carlo workflow and by transient carbon-density input.

### 13.1 Current HPC Python environment

The current production server uses Python 3.8.9:

```bash
source /share/ccsuite/ENV/setenvpython389.sh
source /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/envs/LULCC/bin/activate
```

and:

```bash
export PYTHONPATH="/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/site-packages_py38:/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/Code:${PYTHONPATH:-}"
```

---

## 14. Input and output directories

The current production-server layout is:

```text
/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/
├── Code/                Git repository
├── In_ncfile/           states, transitions and PFT inputs
├── Out_ncfile/          deterministic and MC outputs
├── logs/                DSUB logs
├── envs/LULCC/          Python virtual environment
└── site-packages_py38/  server-local Python packages
```

Relative input filenames in a run YAML are resolved under `paths.data_dir`.

Current principal input names include:

```text
In_ncfile/
├── synthetic_states_1deg.nc
├── synthetic_transitions_1deg.nc
├── synthetic_states_025deg.nc
├── synthetic_transitions_025deg.nc
├── IBIS_PFT_dominant_1deg.nc
└── IBIS_PFT_dominant_0.25deg.nc
```

---

## 15. Running locally

The local entry point reads:

```text
config/run_local.yml
```

Create it from an appropriate run configuration and edit local paths as needed.

Linux or macOS:

```bash
cp config/run_1deg.yml config/run_local.yml
```

Windows:

```text
copy config\run_1deg.yml config\run_local.yml
```

Edit at least:

```text
run.name
run.resolution
run.input_format
run.start_year
run.end_year
paths.data_dir
paths.output_dir
inputs
experiment.path
model_overrides
pft
```

Run from the repository root:

```bash
python -m src.main
```

The local deterministic output is:

```text
Out_ncfile/<resolution>/<run.name>/<output.filename_prefix>.global.nc
```

---

## 16. Running deterministic simulations on the HPC server

The current production scheduler is DSUB.

Run from:

```bash
cd /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/Code
```

The run YAML controls the actual model setup. The shell scripts select the resolution, band ID and server environment.

### 16.1 1° global simulation

With:

```yaml
server:
  band_size_deg: 10
```

the globe is divided into 36 longitude bands:

```text
BAND_ID = 0 ... 35
```

Submit all bands:

```bash
bash server/submit_1deg.sh
```

Submit a subset, for example bands 0-5:

```bash
bash server/submit_1deg.sh 0 5
```

The submit helper generates one DSUB script per band under:

```text
server/generated_1deg/
```

and submits it with:

```text
dsub -s
```

### 16.2 0.25° global simulation

With:

```yaml
server:
  band_size_deg: 3
```

the globe is divided into 120 longitude bands:

```text
BAND_ID = 0 ... 119
```

Submit all bands:

```bash
bash server/submit_025deg.sh
```

Submit a subset:

```bash
bash server/submit_025deg.sh 0 11
```

### 16.3 Single-band test or rerun

A single 1° band:

```bash
bash server/submit_1deg.sh 0 0
```

A single 0.25° band:

```bash
bash server/submit_025deg.sh 0 0
```

### 16.4 Job status and logs

Use the DSUB job-status tools provided by the HPC platform, for example:

```text
djob
```

Logs are written under:

```text
/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/
```

### 16.5 Rerunning completed or incomplete bands

A deterministic band is considered reusable only when its dimensions, stored run identity and `done` array are complete and current.

An incomplete or stale output is removed and recomputed from the beginning.

To intentionally create a scientifically distinct simulation, use a new `run.name`.

### 16.6 Cleaning deterministic outputs

Examples:

```bash
bash server/clean_run.sh 1deg baseline_1deg_voidpft_v1
bash server/clean_run.sh 025deg baseline_025deg_v1
```

Delete all model outputs:

```bash
bash server/clean_run.sh all
```

Logs are retained by default. To intentionally remove matching logs:

```bash
CLEAN_LOGS=1 bash server/clean_run.sh ...
```

---

## 17. Deterministic output structure

Outputs are stored under:

```text
Out_ncfile/<resolution>/<run.name>/
```

### Local run

```text
<output.filename_prefix>.global.nc
```

### Server run

```text
<output.filename_prefix>.rank000.nc
<output.filename_prefix>.rank001.nc
...
```

The run directory can also contain:

```text
run_config.yml
run_manifest.json
run_manifest.rankXXX.json
input_validation_report.json
```

### 17.1 NetCDF global attributes

The NetCDF file records model provenance and runtime metadata, including:

- run name;
- resolution;
- input format;
- simulation years;
- state, transition and PFT paths;
- parameter and experiment paths;
- PFT source mode;
- PFT update mode;
- harvest mode;
- dynamic-density status;
- PFT names and count;
- area and carbon units;
- run-configuration hash;
- parameter and experiment hashes;
- Git commit;
- full run YAML;
- band ID and band size for server runs.

---

## 18. Output variables

### 18.1 Carbon stocks

```text
biomass_total
soil_total
P1
P10
P100
atmosphere
system_carbon_total
```

### 18.2 Carbon budget and closure

```text
Gross_Sources
Gross_Sinks
Net_Emissions
Closure_Error
Atmosphere_Closure_Error
Carbon_Density_Adjustment
PFT_Remap_Adjustment
```

Definitions:

```text
Closure_Error
= annual change in system_carbon_total
- Carbon_Density_Adjustment
- PFT_Remap_Adjustment
```

```text
Atmosphere_Closure_Error
= Net_Emissions
- annual atmosphere-pool increment
```

### 18.3 LUCE-style categories

```text
Flux_FD
Flux_NFC
Flux_FR
Flux_NFR
Flux_CAL
Flux_WHp
```

### 18.4 Event-level fluxes

```text
Flux_Clearing
Flux_Abandonment
Flux_Harvest_Net
Flux_Harvest_SoilSlash
Flux_Harvest_Regrowth
Flux_Other
Flux_Products_Total
```

The model is designed so that:

```text
Flux_Clearing
+ Flux_Abandonment
+ Flux_Harvest_Net
+ Flux_Other
= Net_Emissions
```

within numerical precision.

### 18.5 Transition-area diagnostics

```text
transition_requested_area
transition_applied_area
transition_clipped_area
area_clearing
area_abandonment
area_other
area_harvest
area_deforestation
```

### 18.6 Harvest diagnostics

```text
harvest_requested_biomass
harvest_met_biomass
harvest_unmet_biomass
harvest_unmet_raw_biomass
harvest_forced_biomass

harvest_luh2_area_frac
harvest_luh2_area
harvest_effective_area_frac
harvest_effective_area
harvest_extra_area

harvest_biomass_removed
harvest_patch_biomass_before
harvest_loss_biomass

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

Variables that do not apply to the selected experiment are retained as zero to preserve a consistent output schema.

---

## 19. Units

When:

```text
area_unit = ha
carbon density = t C ha^-1
```

carbon stocks and fluxes are:

```text
t C
```

Area variables are:

```text
ha
```

Fraction variables are dimensionless:

```text
1
```

Common conversions are:

```text
Tg C = t C / 1e6
Pg C = t C / 1e9
```

For area-normalized maps:

```text
g C m^-2 yr^-1
= grid-cell flux in t C yr^-1 × 1e6
  / grid-cell area in m²
```

---

## 20. Model cautions

- The land-use input format in the run YAML must match the actual NetCDF variables.
- The PFT map and parameter file must use the same PFT count and class order.
- Forest PFT indices determine forest/non-forest harvest allocation and LUCE process attribution.
- Original LUH2 and VSCP transition files must retain required harvest variables.
- Subsequent land-cover fractions are driven by transitions; the states file is used to initialize the selected start year.
- Dynamic carbon density introduces an external stock adjustment and should be interpreted separately from `Net_Emissions`.
- Dynamic PFT remapping changes equilibrium composition but is designed to conserve ecosystem carbon at the remapping step.
- Urban transitions are not yet a complete carbon-active process.
- A short simulation omits legacy emissions and sinks generated before the selected start year.
- Use a new `run.name` for every scientifically distinct set of inputs, parameters or model options.
- A complete deterministic band is not resumed grid by grid after interruption; incomplete output is recomputed from the beginning.
- Input validation should be run after replacing land-use, PFT or carbon-density data.
- Monte Carlo uncertainty currently represents parameter uncertainty conditional on the selected land-use forcing, PFT input and model structure; it is not a complete estimate of all possible structural and forcing uncertainty.

---

## 21. Recommended deterministic workflow

```text
1. Prepare states, transitions and PFT inputs.
2. Select original_luh2 or vscp input format.
3. Select a harvest experiment YAML.
4. Create or edit the run YAML.
5. Assign a unique run.name.
6. Validate the inputs.
7. Run a one-band test.
8. Inspect Closure_Error, transition_clipped_area and harvest diagnostics.
9. Submit all longitude bands.
10. Check all bands before merging outputs.
11. Merge or summarize the deterministic result.
```

---

## 22. Monte Carlo uncertainty analysis

The repository includes a summary-only Monte Carlo workflow for parameter uncertainty analysis at both 1° and 0.25°.

The deterministic bookkeeping core is not duplicated. Each Monte Carlo realization:

1. reads one reproducible parameter realization from a sample table;
2. merges the sampled overrides with the normal model configuration;
3. passes the effective parameters through the existing `ParameterLoader`;
4. calls the normal `LULCCSimulator`;
5. spatially sums annual additive outputs within each longitude band;
6. stores compact Parquet summaries rather than a full gridded NetCDF for every realization.

### 22.1 MC configurations

Run configurations:

```text
config/run_mc_1deg.yml
config/run_mc_025deg.yml
```

Scientific sampling configurations:

```text
config/mc_parameters_1deg.yml
config/mc_parameters_025deg.yml
```

Software smoke-test configurations:

```text
config/mc_parameters_smoke_1deg.yml
config/mc_parameters_smoke_025deg.yml
```

Core MC code:

```text
src/mc_sampling.py
src/mc_runner.py
```

Post-processing:

```text
tools/generate_mc_samples.py
tools/validate_mc_band.py
tools/merge_mc_bands.py
tools/merge_all_mc.py
tools/summarize_mc.py
tools/check_mc_convergence.py
```

The MC run YAMLs should remain synchronized with the corresponding deterministic run YAMLs except for `run.name` and the `monte_carlo` section.

### 22.2 Current MC baseline setup

#### 1°

```text
input_format: vscp
states: synthetic_states_1deg.nc
transitions: synthetic_transitions_1deg.nc
pft: IBIS_PFT_dominant_1deg.nc
missing_cell_policy: skip
dynamic_carbon_density.enabled: false
band_size_deg: 10
36 bands
```

#### 0.25°

```text
input_format: vscp
states: synthetic_states_025deg.nc
transitions: synthetic_transitions_025deg.nc
pft: IBIS_PFT_dominant_0.25deg.nc
missing_cell_policy: skip
dynamic_carbon_density.enabled: false
band_size_deg: 3
120 bands
```

### 22.3 Scientific parameter design

The first production MC design prioritizes carbon density and carbon allocation, followed by response and recovery times.

| Parameter family | Distribution | Working uncertainty |
|---|---|---|
| Natural biomass density (`v`,`s`) | shared lognormal multiplier | CV = 20% |
| Managed biomass density (`c`,`p`) | shared lognormal multiplier | CV = 20% |
| Natural soil density (`v`,`s`) + harvest `SOC_min_v/s` | shared lognormal multiplier | CV = 15% |
| Managed soil density (`c`,`p`) | shared lognormal multiplier | CV = 15% |
| Clearing response time | shared lognormal multiplier | CV = 25% |
| Abandonment biomass time | shared lognormal multiplier | CV = 25% |
| Abandonment soil time | shared lognormal multiplier | CV = 25% |
| Harvest response time | shared lognormal multiplier | CV = 25% |
| Fast/slow fractions | beta | concentration = 40 |
| Clearing allocation | Dirichlet/simplex | concentration = 60 |
| Harvest product split | Dirichlet/simplex | concentration = 60 |

These values are literature-informed working priors, not probability distributions directly reported by the cited bookkeeping studies.

Keep:

```yaml
scientific_ranges_confirmed: false
```

until the uncertainty design has been reviewed. Change it to `true` only for the approved production experiment.

Use the same priors at 1° and 0.25° so spatial resolution is not confounded with a different parameter distribution.

### 22.4 Coupling `SOC_min` to equilibrium soil density

`SOC_min_v` and `SOC_min_s` are absolute soil-carbon-density floors used during harvest disturbance.

If equilibrium soil density and `SOC_min` were sampled independently, some realizations could generate:

```text
SOC_min > equilibrium Soil
```

which would suppress modeled harvest soil loss.

Therefore the production parameter design applies the same shared random multiplier to:

```text
carbon_density.PFT*.Soil.v
carbon_density.PFT*.Soil.s
harvest_param.PFT*.SOC_min_v
harvest_param.PFT*.SOC_min_s
```

This preserves the configured ratio:

```text
SOC_min / equilibrium Soil
```

within every realization.

The MC parameter validator should additionally enforce:

```text
SOC_min_v <= Soil.v
SOC_min_s <= Soil.s
```

### 22.5 Syntax check

After changing MC code, run from the repository root:

```bash
python -m py_compile \
  src/mc_sampling.py \
  src/mc_runner.py \
  tools/generate_mc_samples.py \
  tools/validate_mc_band.py \
  tools/merge_mc_bands.py \
  tools/merge_all_mc.py \
  tools/summarize_mc.py \
  tools/check_mc_convergence.py \
  server/main_mc_1deg.py \
  server/main_mc_025deg.py
```

### 22.6 Validate MC inputs

Because the current `tools/validate_inputs.py` uses a `CONFIG_PATH` constant, set:

```python
CONFIG_PATH = Path("config/run_mc_1deg.yml")
```

and run:

```bash
python -m tools.validate_inputs
```

For 0.25° use:

```python
CONFIG_PATH = Path("config/run_mc_025deg.yml")
```

Validation should be run once before submitting many DSUB jobs.

### 22.7 Generate smoke-test samples

The smoke test checks the computational pipeline only. It is not a scientific uncertainty experiment.

For 1°:

```bash
python -m tools.generate_mc_samples \
  --spec config/mc_parameters_smoke_1deg.yml \
  --allow-unconfirmed-ranges \
  --force
```

For 0.25°:

```bash
python -m tools.generate_mc_samples \
  --spec config/mc_parameters_smoke_025deg.yml \
  --allow-unconfirmed-ranges \
  --force
```

The smoke table contains:

```text
sample_id 0   deterministic baseline
sample_id 1   small perturbation
sample_id 2   small perturbation
```

### 22.8 Mandatory baseline-equivalence test

Before a large ensemble, run the deterministic MC baseline for one band.

For 1° band 0:

```bash
bash server/submit_mc_1deg.sh 0 0 0 0
```

The MC output is located under:

```text
Out_ncfile/MC/1deg/mc_1deg_static_vscp_v1/
└── sample_000000/
    └── bands/
        ├── band_000.parquet
        └── band_000.json
```

Compare this band against the matching deterministic 1° `rank000.nc` using `tools.validate_mc_band`.

Check the current CLI if needed:

```bash
python -m tools.validate_mc_band --help
```

Do not launch the full scientific ensemble until all common additive variables pass within numerical tolerance.

### 22.9 Generate production samples

After the parameter ranges have been approved and:

```yaml
scientific_ranges_confirmed: true
```

generate 1° samples:

```bash
python -m tools.generate_mc_samples \
  --spec config/mc_parameters_1deg.yml \
  --force
```

or 0.25° samples:

```bash
python -m tools.generate_mc_samples \
  --spec config/mc_parameters_025deg.yml \
  --force
```

The current default design generates:

```text
500 random realizations
+ sample 0 deterministic baseline
= 501 rows
```

The sample table records, for each realization:

```text
sample_id
sample_kind
sample_seed
override_sha256
overrides_json
draws_json
```

The same master seed and unchanged specification produce reproducible samples.

### 22.10 DSUB MC submission

Each MC DSUB job currently represents:

```text
one sample_id × one longitude band
```

and uses one CPU.

#### 1°

Samples 0-2 across all 36 bands:

```bash
bash server/submit_mc_1deg.sh 0 2
```

Samples 1-10:

```bash
bash server/submit_mc_1deg.sh 1 10
```

Specific sample and band subset:

```bash
bash server/submit_mc_1deg.sh 1 10 0 5
```

#### 0.25°

Samples 0-2 across all 120 bands:

```bash
bash server/submit_mc_025deg.sh 0 2
```

Samples 1-5:

```bash
bash server/submit_mc_025deg.sh 1 5
```

Specific band subset:

```bash
bash server/submit_mc_025deg.sh 1 5 0 11
```

#### Generic wrapper

```bash
bash server/submit_mc_batch.sh 1deg   1 10
bash server/submit_mc_batch.sh 025deg 1 5
```

Re-submission is safe because completed bands are reused only when their result identity is current.

### 22.11 MC output structure

The MC output root is:

```text
Out_ncfile/MC/<resolution>/<mc_name>/
```

A realization is stored as:

```text
sample_000001/
├── bands/
│   ├── band_000.parquet
│   ├── band_000.json
│   ├── ...
│   └── band_NNN.parquet
├── global.parquet
└── sample_manifest.json
```

Band Parquet files contain:

```text
year
+ selected additive annual model outputs
```

The summary-only design avoids storing one full global gridded NetCDF for every Monte Carlo realization.

### 22.12 Reproducibility and stale-output protection

An existing MC band or merged realization is reused only if the current result-affecting identity matches the stored metadata.

The identity includes:

```text
sample-specific override hash
run-YAML hash
parameter-file hash
experiment-file hash
state-file fingerprint
transition-file fingerprint
PFT-file fingerprint
model-code hash
year range
output-variable list
```

The model-code hash is calculated directly from result-affecting Python files, so uncommitted local code changes also invalidate old results.

The sample-table SHA and Git commit are retained as provenance metadata but are not used alone to invalidate a realization. This allows an ensemble to be extended with additional samples without forcing unchanged earlier sample IDs to rerun.

If a band or merged output is incomplete, stale, corrupted, or inconsistent with the current identity, it is not silently reused.

### 22.13 Merge one realization

For 1° sample 1:

```bash
python -m tools.merge_mc_bands \
  --config config/run_mc_1deg.yml \
  --sample-id 1
```

For 0.25°, change the config:

```bash
python -m tools.merge_mc_bands \
  --config config/run_mc_025deg.yml \
  --sample-id 1
```

The merger checks all expected bands and rejects missing or stale band metadata.

### 22.14 Merge many realizations

While computation is still in progress:

```bash
python -m tools.merge_all_mc \
  --config config/run_mc_1deg.yml \
  --start 0 \
  --end 500 \
  --skip-incomplete
```

After all expected jobs finish, repeat without `--skip-incomplete`:

```bash
python -m tools.merge_all_mc \
  --config config/run_mc_1deg.yml \
  --start 0 \
  --end 500
```

The final strict run should fail if any expected realization cannot be merged.

### 22.15 Ensemble summary

For 1°:

```bash
python -m tools.summarize_mc \
  --config config/run_mc_1deg.yml
```

For 0.25°:

```bash
python -m tools.summarize_mc \
  --config config/run_mc_025deg.yml
```

The deterministic baseline is excluded from the MC distribution by default.

For every year and selected output variable, the ensemble summary contains:

```text
n
mean
sd
p05
p50
p95
```

Principal summary files include:

```text
ensemble_summary.parquet
ensemble_summary.json
net_emissions_samples.parquet
```

### 22.16 Convergence analysis

The largest completed ensemble is used as a reference distribution only.

It must not be allowed to demonstrate convergence by comparing against itself.

For example, with 500 completed random realizations, test only candidate sizes smaller than 500:

```bash
python -m tools.check_mc_convergence \
  --config config/run_mc_1deg.yml \
  --sizes 25,50,100,200,300,400 \
  --start-year 1850 \
  --end-year 2020 \
  --replicates 50
```

The convergence analysis evaluates:

```text
annual mean trajectory
annual P05 trajectory
annual P50 trajectory
annual P95 trajectory
cumulative mean
cumulative SD
cumulative P05
cumulative P50
cumulative P95
```

Current working pass thresholds are:

```text
cumulative mean error       <= 2%
cumulative SD/P05/P95 error <= 5%
annual statistic NRMSE      <= 5%
```

Outputs are:

```text
convergence_summary.parquet
convergence_summary.csv
convergence_summary.json
```

If at least one candidate size smaller than the full ensemble passes all criteria, the smallest passing value is reported as:

```text
recommended_minimum_n
```

If none passes:

```text
convergence_demonstrated: false
recommended_minimum_n: null
```

In that case, add more Monte Carlo realizations and repeat the convergence analysis. The full ensemble must not be considered converged merely because it matches itself.

### 22.17 Production strategy across resolutions

The current one-job-per-sample-per-band design implies:

```text
1°:
501 × 36 = 18,036 jobs

0.25°:
501 × 120 = 60,120 jobs
```

Therefore the recommended strategy is:

```text
1. Validate the 1° MC inputs.
2. Generate 1° smoke samples.
3. Run sample 0 / band 0.
4. Pass the baseline-equivalence test.
5. Review and approve scientific priors.
6. Generate the 1° production sample table.
7. Run the 1° ensemble in manageable DSUB batches.
8. Merge and summarize the 1° ensemble.
9. Evaluate MC convergence.
10. Use the convergence result to guide the initial 0.25° sample count.
11. Run the corresponding 0.25° sample IDs.
12. Extend the 0.25° ensemble if its own convergence check requires more samples.
```

The 1° convergence result is useful for computational planning but does not mathematically prove convergence at 0.25°; the 0.25° ensemble should also be checked once enough realizations are available.

### 22.18 Interpreting the uncertainty result

The current MC workflow quantifies parameter uncertainty conditional on:

```text
selected LUH2/VSCP land-use forcing
selected PFT map
selected static carbon-density framework
selected harvest experiment
current bookkeeping model structure
```

It does not automatically include:

```text
land-use forcing uncertainty
PFT-map uncertainty
alternative bookkeeping structures
dynamic/transient carbon-density uncertainty
other structural model uncertainty
```

Therefore reported MC intervals should be described as model-parameter uncertainty for the specified experiment rather than total uncertainty from all possible sources.

For publication-level reporting, record at minimum:

```text
sampling specification
parameter distributions and priors
master random seed
number of realizations
baseline configuration
mean and/or median
standard deviation
P05-P95 interval
convergence diagnostic
model and input-data version/provenance
```

---

## 23. Recommended complete workflow

### Deterministic model

```text
Prepare inputs
→ validate
→ one-band test
→ full-band simulation
→ inspect closure and process diagnostics
→ merge/output analysis
```

### Monte Carlo parameter uncertainty

```text
Validate MC configuration
→ smoke samples
→ baseline-equivalence test
→ approve priors
→ generate production sample table
→ 1° ensemble
→ merge
→ summarize
→ convergence
→ staged 0.25° ensemble
→ final uncertainty reporting
```
