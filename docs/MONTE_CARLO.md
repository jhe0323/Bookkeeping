# Monte Carlo workflow for Bookkeeping

This upgrade is deliberately additive. It does **not** replace the existing
`LULCCSimulator.py`, `run_manager.py`, deterministic run YAMLs, or deterministic
SLURM script. The established model path remains unchanged.

## 1. Files added by the upgrade

- `src/mc_sampling.py` — generic reproducible parameter sampler, wildcard target
  expansion, Dirichlet/simplex allocation sampling, and physical constraints.
- `src/mc_runner.py` — summary-only Monte Carlo band runner. It calls the existing
  `LULCCSimulator.run_simulation()` and accumulates annual band totals in memory.
- `tools/generate_mc_samples.py` — generates `mc_samples.parquet` plus a sampling
  manifest containing seeds and hashes.
- `tools/merge_mc_bands.py` — merges the 36 one-degree longitude bands for one
  realization.
- `tools/merge_all_mc.py` — merges a range/all completed realizations.
- `tools/summarize_mc.py` — calculates mean, SD, P05, P50 and P95 across MC runs.
- `tools/validate_mc_band.py` — baseline equivalence test against the existing
  deterministic gridded NetCDF output.
- `server/main_mc_1deg.py` — one-degree MC entry point.
- `server/run_mc_band.sbatch` — maps each array task to `(sample_id, band_id)`.
- `server/submit_mc_batch.sh` — submits manageable sample batches.
- `config/run_mc_1deg.yml` — MC-specific run configuration; deterministic configs
  are untouched.
- `config/mc_parameters.yml` — scientific sampling template (disabled examples).
- `config/mc_parameters_smoke.yml` — tiny +/-1% SOFTWARE TEST only.

## 2. Merge into the current repository

From the repository root on the server:

```bash
cd /mnt/beegfs/product/lulc0120/Code
git checkout -b mc-upgrade
```

Copy the upgrade package into this directory while preserving `src/`, `tools/`,
`server/`, `config/`, and `docs/`. No existing file should need to be overwritten.
Then inspect:

```bash
git status
```

Expected result: only new untracked MC files.

## 3. Syntax and input checks

```bash
python3 -m py_compile \
  src/mc_sampling.py \
  src/mc_runner.py \
  tools/generate_mc_samples.py \
  tools/merge_mc_bands.py \
  tools/merge_all_mc.py \
  tools/summarize_mc.py \
  tools/validate_mc_band.py \
  server/main_mc_1deg.py
```

Validate the model inputs once before submitting an array:

```bash
python3 -m tools.validate_inputs --config config/run_mc_1deg.yml
```

`validation.run_before_simulation` stays false in the MC run YAML so every SLURM
band does not rescan the same global inputs.

## 4. Software smoke test

Generate two tiny perturbed realizations plus an unchanged baseline:

```bash
python3 -m tools.generate_mc_samples \
  --spec config/mc_parameters_smoke.yml \
  --allow-unconfirmed-ranges \
  --force
```

The output sample IDs are:

- `0`: deterministic baseline (no sample overrides)
- `1-2`: software-test perturbations

Run only baseline band 0 first:

```bash
MC_SAMPLE_ID=0 BAND_ID=0 \
RUN_CONFIG=/mnt/beegfs/product/lulc0120/Code/config/run_mc_1deg.yml \
python3 -m server.main_mc_1deg
```

The output is small and located at:

```text
../Out_ncfile/MC/1deg/mc_1deg_static_v1/
  sample_000000/
    bands/
      band_000.parquet
      band_000.json
```

## 5. Mandatory baseline-equivalence test

Compare the summary-only baseline band with the already-tested deterministic
NetCDF band. For the current `run_1deg.yml` naming convention an example is:

```bash
python3 -m tools.validate_mc_band \
  --mc-band ../Out_ncfile/MC/1deg/mc_1deg_static_v1/sample_000000/bands/band_000.parquet \
  --grid-nc ../Out_ncfile/1deg/baseline_1deg_area_pftfallback_v1/summary_1deg.rank000.nc
```

All common variables should report `PASS`. Do not start a large ensemble until
this equivalence test passes.

## 6. Define the scientific uncertainty experiment

Edit:

```text
config/mc_parameters.yml
```

The included uncertainty widths are examples only and are disabled. For every
scientific parameter group:

1. choose the parameter(s) supported by the uncertainty analysis;
2. replace the example distribution/range with literature- or data-supported
   values;
3. set `enabled: true` for the selected rules;
4. set `scientific_ranges_confirmed: true` only after reviewing all ranges.

Important parameter rules implemented in the sampler:

- positive quantities can use multiplicative uniform/normal/lognormal sampling;
- fractions can use beta sampling centered on the baseline and remain in [0,1];
- clearing allocation fractions can use a Dirichlet/simplex group, preserving
  their baseline total and preventing an unphysical negative residual;
- harvest product fractions can also use a simplex group;
- all generated samples are checked for non-negative carbon densities/timescales,
  valid [0,1] fractions, clearing allocation sum <= 1, and harvest product sum <= 1.

The first MC experiment is configured with dynamic carbon density disabled. Do
not sample the static `carbon_density` table while transient density is enabled.
A transient-density uncertainty experiment should instead perturb the time-varying
density series/multipliers explicitly.

## 7. Generate the production samples

For 500 random realizations plus one baseline (`sample_id 0`):

```bash
python3 -m tools.generate_mc_samples \
  --spec config/mc_parameters.yml \
  --force
```

Outputs:

```text
config/mc_samples.parquet
config/mc_samples.parquet.manifest.json
```

Each row stores `sample_id`, sample seed, parameter override JSON, draw JSON, and
an SHA-256 hash. The same sample table therefore reproduces the same ensemble.

## 8. Submit HPC jobs in batches

At 1 degree there are 36 longitude bands. The submission helper maps an array
index to:

```text
sample_id = sample_offset + task_id // 36
band_id   = task_id % 36
```

Example batches:

```bash
bash server/submit_mc_batch.sh 0 25 35
bash server/submit_mc_batch.sh 26 50 35
bash server/submit_mc_batch.sh 51 75 35
```

The third argument is the maximum number of concurrently running array tasks.
Continue the batches through sample 500. Using batches avoids one extremely large
18,000-task array and makes failures/restarts easier to manage.

A completed band is automatically skipped when its Parquet output and metadata
match the sample hash and run-config hash. Re-submitting a failed batch is safe.

## 9. Merge band summaries

Merge one sample:

```bash
python3 -m tools.merge_mc_bands \
  --config config/run_mc_1deg.yml \
  --sample-id 0
```

Merge all completed samples while jobs are still finishing:

```bash
python3 -m tools.merge_all_mc \
  --config config/run_mc_1deg.yml \
  --start 0 --end 500 \
  --skip-incomplete
```

After all jobs finish, run once without `--skip-incomplete`; this makes missing
bands an error and is a useful completeness check.

Each realization then has:

```text
sample_XXXXXX/global.parquet
sample_XXXXXX/sample_manifest.json
```

## 10. Calculate ensemble uncertainty

```bash
python3 -m tools.summarize_mc --config config/run_mc_1deg.yml
```

The deterministic baseline is excluded from the MC distribution by default.
Outputs:

```text
../Out_ncfile/MC/1deg/mc_1deg_static_v1/ensemble_summary.parquet
../Out_ncfile/MC/1deg/mc_1deg_static_v1/ensemble_summary.json
../Out_ncfile/MC/1deg/mc_1deg_static_v1/net_emissions_samples.parquet
```

`ensemble_summary.parquet` contains, for every year and selected variable:

- `n`
- `mean`
- `sd`
- `p05`
- `p50`
- `p95`

Carbon variables retain model units (`t C`). Convert to Pg C with `1e-9` only in
the downstream analysis/plotting step.

## 11. Git commit

After the baseline equivalence test succeeds:

```bash
git add src/mc_sampling.py src/mc_runner.py \
  tools/generate_mc_samples.py tools/merge_mc_bands.py \
  tools/merge_all_mc.py tools/summarize_mc.py tools/validate_mc_band.py \
  server/main_mc_1deg.py server/run_mc_band.sbatch server/submit_mc_batch.sh \
  config/run_mc_1deg.yml config/mc_parameters.yml config/mc_parameters_smoke.yml \
  docs/MONTE_CARLO.md

git commit -m "feat: add Monte Carlo sampling and summary-only HPC workflow"
```

## Design principle

The deterministic model core is intentionally not modified. The MC layer only:

1. creates parameter overrides;
2. passes them through the existing `ParameterLoader` override mechanism;
3. calls the existing cell simulator;
4. spatially accumulates additive annual outputs;
5. records sample/band hashes and seeds;
6. calculates ensemble statistics after the runs.

This isolates Monte Carlo engineering from the BLUE/LUCE carbon-process logic and
makes it much easier to verify that an unchanged baseline produces unchanged
model results.
