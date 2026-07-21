# LULCC Carbon Bookkeeping Model

A Python-based land-use and land-cover change (LULCC) carbon bookkeeping model for estimating annual carbon-stock changes, gross sources and sinks, and net land-use change emissions.

The model follows the general bookkeeping framework used by BLUE- and LUCE-type models. Land-use events create departures from equilibrium biomass and soil carbon stocks. These departures are retained in legacy pools and relax over time, allowing historical land-use change to continue affecting present-day carbon fluxes.

The model is driven by:

- land-cover state data;
- annual land-use transition data;
- wood-harvest area or biomass-demand data;
- plant functional type (PFT) maps;
- PFT- and land-cover-specific biomass and soil carbon densities;
- event-specific allocation and response parameters.

The current `main` branch supports both original LUH2 variables and pre-aggregated `v/s/c/p/U` inputs through run configuration files. Separate code branches are not required.

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

For static carbon density and conservative PFT remapping, the system closure error should be close to numerical zero.

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

The mapping is configurable and should be kept consistent with the PFT map and carbon-density parameterization.

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

The target agricultural biomass stock is treated as an immediate biomass deficit and subsequently relaxes toward equilibrium.

### 3.2 Abandonment and regrowth

Abandonment includes:

```text
c → s
p → s
```

The source biomass and slow-soil state are inherited by secondary vegetation. The difference between the inherited state and the secondary-vegetation equilibrium state creates biomass and soil recovery trajectories.

### 3.3 Cropland–pasture conversion

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
```

The current transient-density implementation refreshes equilibrium biomass and soil carbon for the existing land-cover-by-PFT area each year.

The resulting external stock adjustment is reported as:

```text
Carbon_Density_Adjustment
```

This adjustment is excluded from `Closure_Error` and is not included in `Net_Emissions`.

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

A complete simulation is defined by a run YAML.

Current server configurations are:

```text
config/run_025deg.yml
config/run_1deg.yml
```

A run configuration contains:

```text
run             simulation identity, resolution, input format and years
paths           input and output directories
inputs          state, transition and PFT files
parameters      base parameter file
experiment      harvest experiment file
model_overrides run-specific parameter switches
pft             PFT time handling and update mode
validation      input validation settings
output          output prefix, compression and synchronization
server          longitude-band width and start staggering
```

Example:

```yaml
run:
  name: baseline_1deg_area
  resolution: 1deg
  input_format: original_luh2
  input_base_year: 850
  start_year: 850
  end_year: 2020
  area_unit: ha

paths:
  data_dir: ../In_ncfile
  output_dir: ../Out_ncfile

inputs:
  states: states_1deg.nc
  transitions: transitions_1deg.nc
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

Use a unique `run.name` whenever any input, experiment, model switch or parameter file changes. Completed outputs are skipped when their dimensions and `done` flags are complete.

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

### 11.2 Standalone validation

Edit:

```text
tools/validate_inputs.py
```

and set:

```python
CONFIG_PATH = Path("config/run_025deg.yml")
```

Then run:

```bash
python -m tools.validate_inputs
```

Use `full_scan: false` for a faster structural check and `full_scan: true` for a complete data scan.

---

## 12. Repository structure

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
│   └── summary_yearly.py
│
├── config/
│   ├── config.yml
│   ├── dynamic_carbon_density.parquet
│   ├── run_025deg.yml
│   ├── run_1deg.yml
│   └── experiments/
│       ├── harvest_area.yml
│       ├── harvest_bio_strict.yml
│       └── harvest_bio_forced.yml
│
├── server/
│   ├── __init__.py
│   ├── main_025deg.py
│   ├── main_1deg.py
│   ├── run_single_band.sbatch
│   └── clean_run.sh
│
├── tools/
│   └── validate_inputs.py
│
└── docs/
    ├── CHANGELOG.txt
    └── ROADMAP.md
```

### Main modules

| File | Role |
|---|---|
| `src/main.py` | Local entry point using `config/run_local.yml` |
| `src/run_config.py` | Loads and resolves run YAML files |
| `src/run_manager.py` | Shared local/server orchestration, output checks and manifests |
| `src/input_validator.py` | One-time input and configuration validation |
| `src/LULCCSimulator.py` | Annual bookkeeping simulation and NetCDF output |
| `src/file_loader.py` | NetCDF slicing, time conversion, grid alignment and PFT loading |
| `src/parameter_loader.py` | Base parameters, experiment overlays and dynamic density |
| `src/carbon_pools_init.py` | Pool definitions and initialization |
| `src/events.py` | Clearing, abandonment and cropland–pasture events |
| `src/harvest.py` | Area- and biomass-driven harvest |
| `src/transition.py` | Annual decay, regrowth and flux attribution |
| `src/summary_yearly.py` | Annual carbon-stock summaries |
| `server/main_025deg.py` | 0.25° server entry point |
| `server/main_1deg.py` | 1° server entry point |
| `server/run_single_band.sbatch` | SLURM longitude-band launcher |

---

## 13. Installation

The current source code requires Python 3.10 or later.

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

Main dependencies are:

```text
numpy
netCDF4
PyYAML
pandas
pyarrow
```

`pandas` and `pyarrow` are required when transient carbon density is enabled.

---

## 14. Input and output directories

The default server layout is:

```text
/mnt/beegfs/product/lulc0120/
├── Code/          Git repository
├── In_ncfile/     State, transition and PFT inputs
├── Out_ncfile/    Model outputs
├── logs/          SLURM logs
└── site-packages/ Server-local Python packages
```

Relative filenames in a run YAML are resolved under `paths.data_dir`.

Example:

```text
In_ncfile/
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

---

## 15. Running locally

The local entry point reads:

```text
config/run_local.yml
```

Create the file from one of the server configurations:

Linux or macOS:

```bash
cp config/run_1deg.yml config/run_local.yml
```

Windows:

```text
copy config\run_1deg.yml config\run_local.yml
```

Edit `config/run_local.yml`, especially:

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

The local output is:

```text
Out_ncfile/<resolution>/<run.name>/<output.filename_prefix>.global.nc
```

---

## 16. Running on the server

Edit the matching run configuration before submission:

```text
0.25° → config/run_025deg.yml
1°    → config/run_1deg.yml
```

The SLURM script selects only the resolution. The harvest experiment and all model options are selected inside the run YAML.

Run from:

```bash
cd /mnt/beegfs/product/lulc0120/Code
```

### 16.1 1° global simulation

With:

```yaml
server:
  band_size_deg: 10
```

the globe is divided into 36 bands:

```bash
sbatch \
  --array=0-35%25 \
  --export=ALL,RESOLUTION=1deg \
  server/run_single_band.sbatch
```

### 16.2 0.25° global simulation

With:

```yaml
server:
  band_size_deg: 3
```

the globe is divided into 120 bands:

```bash
sbatch \
  --array=0-119%35 \
  --export=ALL,RESOLUTION=025deg \
  server/run_single_band.sbatch
```

### 16.3 Single-band test

1° band 0:

```bash
sbatch \
  --array=0-0 \
  --export=ALL,RESOLUTION=1deg \
  server/run_single_band.sbatch
```

0.25° band 0:

```bash
sbatch \
  --array=0-0 \
  --export=ALL,RESOLUTION=025deg \
  server/run_single_band.sbatch
```

### 16.4 Rerunning a band

An output is skipped when:

- time, latitude and longitude dimensions match the current run;
- the complete `done` array equals 1.

An incomplete output is removed and recomputed from the beginning.

To intentionally rerun a completed simulation:

- use a new `run.name`; or
- remove the completed output file or run directory before submission.

---

## 17. Output structure

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

The run directory also contains:

```text
run_config.yml
run_manifest.json
input_validation_report.json   when validation is enabled
```

### 17.1 NetCDF global attributes

The NetCDF file records:

- run name;
- resolution;
- input format;
- input and simulation years;
- state, transition and PFT paths;
- parameter and experiment paths;
- PFT variable and source mode;
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

Fraction variables are:

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
- Use a new `run.name` for every distinct set of inputs, parameters or model options.
- A complete band is not resumed grid by grid after interruption; incomplete output is recomputed from the beginning.
- Input validation should be run after replacing land-use, PFT or carbon-density data.

---

## 21. Recommended workflow

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
10. Check all done variables before merging outputs.
```

