# Development and methodology changes

This document records material changes to the existing
`oelsmann24-verticallandmotion` FACTS2 module in plain language. Git remains the
authoritative record of exact line-by-line changes, while `CHANGELOG.md` summarizes
release-facing changes. The curated `diagnostics/` notebooks document the current
inputs and corrected outputs; ignored `data/discovery/` notebooks preserve historical
baseline comparisons. `notebooks/methodology_diagnostics.ipynb` is a lightweight
methodology summary.

For every substantial scientific or data-processing change, add an entry containing:

- the date and a concise title;
- the scientific motivation and expected effect;
- the source files and functions changed;
- the input datasets affected, including provenance and checksums when relevant;
- whether the CLI or NetCDF contract changed;
- the reference and candidate commands, using identical sites, years, sample count,
  and seed;
- quantitative reference-versus-candidate results;
- tests and container checks run; and
- limitations or unresolved questions.

Do not use this document or a notebook as a substitute for production code, automated
tests, commit messages, or public input-data documentation.

## 2026-11-03: Pre-push review and validation

At the user's request, prepared the existing local changes for commits and a push
to `correct-vlm-methodology`, not `main`. Reviewed production code, tests, notebook
source/outputs, dependencies, workflows and documentation. Checked pending files
for large/generated datasets, workstation paths and credential patterns; the local
agent guide and all raw/NetCDF diagnostic data remain ignored. The container
publication workflow runs only for `main` and version tags, not this feature branch.

All 63 tests pass with the optional diagnostics dependencies installed. Ruff checks,
the Just validation recipe, installed CLI help and `git diff --check` pass. The local
development Docker image builds, retains its non-root default UID 9876 and excludes
the local guide and plotting/notebook workspace. Host and Docker runs for convergence,
constant and GIA-only modes give exactly equal scientific arrays: 101 samples,
seed 1234, base year 2000, 50-year steps through 2300, six requested cities plus a
Pacific offshore point. The generated evidence is kept locally under ignored
`data/diagnostics/review_before_push_20261103/`.

Updated the README example to build/use the local development image, rather than
passing unreleased options to the published image. The example uses explicit host
mount paths, a tab-separated location file and the host user's non-root identity
so it can write to the dedicated output folder. Fixed the changelog comparison link.
No scientific calculations were changed during this final review. Maintainers
still need to review the scientific assumptions and agree on a versioned release
and any separate registry update; no container or release is published here.

## 2026-11-03: Optional GIA-only projections

Added `--gia-only` to the CLI. It bypasses reconstruction, GNSS, coastline and example
inputs and loads only GIA `rsl_mean` and `rsl_sterr`. Both projection paths use the
same quarter-degree GIA interpolation and nearest-site mapping. GIA-only mean and
error coefficients are the respective rates multiplied by elapsed time from the
selected base year, with no VLM trend, GRW or observational blend contribution.
The existing third seeded sampling stream supplies the GIA perturbations, matching
the full model's GIA fallback for identical sites, years and seed.

`vlm_postprocess` and `projection_parameters` accept an optional `gia_only=False`
argument. The default numerical behavior is unchanged. The CLI still requires all
observational inputs when the switch is absent. In GIA-only runs, `--vlm-rate-mode`
is ignored and the output records its effective `vlm_rate_mode` as `gia-only`;
variables, dimensions, units and fill values retain the FACTS local/RSL schema.

The README gives a complete local invocation. Automated tests compare synthetic
output samples against the analytic GIA projection, verify negative/zero/positive
elapsed years, confirm equality with existing GIA fallback, and check that both
rate-mode settings and different sample-block sizes give identical GIA-only results.
They also check missing observational inputs, input overwrite protection, and invalid
GIA inputs. Diagnostics and notebooks are unchanged at the user's request. The
registry, published image and input archive have not been changed; exposing the
new optional flag through FEB is a separate registry update at release.

Validation: all 63 tests pass, Ruff formatting/lint checks pass, and the development
Docker image builds. Host and container runs using the archived `GIA_stats.nc`,
101 samples, seed 1234, base year 2000 and 50-year steps through 2300 produced
exactly equal `(101, 7, 3)` scientific arrays for New York, Skellefteå and an
open-ocean location. All base-year samples were zero and all projections finite.

## 2026-11-02: Full global projection maps, including open-ocean GIA RSL

The first version of the shared projection notebook requested only coastal map
sites. Its blank open ocean was a diagnostic site-selection restriction, not a
zero GIA contribution. At the user's request, `diagnostics/helpers.py` now requests
all 64,800 one-degree cell centers across the globe plus the unchanged six cities.
Cell centers span 89.5°S–89.5°N and 179.5°W–179.5°E; cell edges reach both poles and
the date line. Mean, standard-deviation and difference maps use complete raster
cells rather than scattered coastal dots. No land/ocean/coastal mask hides any
of the diagnostic map cells. Inland RSL values are a model field, not an inland
inundation projection.

Production preprocessing already retains global GIA fields and gives observations
zero weight outside its coastal/distance selection. The production projection
then uses GIA RSL and its uncertainty alone at zero-weight sites. No scientific
formula, source input, production code, CLI or FACTS output schema changes here.
Both modes retain the same inputs, seed, sample count, dates and city coordinates.
The global run has its own ignored `data/diagnostics/shared_notebooks/global_grid/`
directory, preserving the earlier coastal-only products for comparison. Source
hashes and global-grid geometry invalidate incompatible caches. Small optional
tests cover complete globe/date-line/pole coverage, unchanged city locations,
correct value reshaping and rejection of missing/nonfinite map cells. Notebook
assertions also require the complete global grid.

Validation: the projection notebook was rerun successfully and its global and
difference maps checked visually. All 64,800 cells have finite means and standard
deviations at every saved decade. At the 60,013 zero-observation-weight cells,
summary means and spreads match GIA RSL rate times elapsed years and the sampled
GIA uncertainty to float32 storage precision; both modes are exactly identical
there. Direct production coefficients also match GIA exactly at offshore test
points in the Pacific, Atlantic, Indian and Arctic oceans. All six cities retain
exactly the same samples, means and standard deviations as the preceding coastal
run. Selected raw-sample reductions agree with the saved summaries. The production
source hashes are unchanged. All 50 tests pass with the optional diagnostics group,
and Ruff formatting/lint checks and `git diff --check` pass. No new Docker run was
needed for this notebook-only change; the diagnostics folder remains excluded
from the production image.

## 2026-11-02: Versioned input-data and two-mode projection notebooks

At the user's request, added the curated `diagnostics/` folder for inclusion in the
repository. This supersedes the earlier local-only preference for these two
notebooks, without moving or overwriting the historical discovery workspaces.
`01_input_data.ipynb` inventories every input variable and plots reconstruction
trends, trend/GRW uncertainties, spatial footprints, raw and CMR-corrected GNSS,
accepted source points, GIA fields, the coastal mask and spatial blend weights.
`02_projection_comparison.ipynb` uses only the corrected production model and
compares `converge-to-gia` with `constant` secular rates. The six city locations are
Venice, Skellefteå, London, Tokyo, Brisbane and New York. London replaces Valparaíso;
coordinates are explicitly listed as representative city locations, not exact
paper tide gauges. Embedded executed figures can be read without downloading data.

Supporting `helpers.py` calls the actual preprocessing and projection functions;
it does not copy or change the scientific implementation. Both modes use identical
inputs, sites, seed 1234, 500 samples, base year 2000 and ten-year output steps to
2300. The default maps compare 2150 mean and population standard deviation on shared
color scales, with a separate difference map. There are 14,242 integer-degree nodes
inside the corrected 150 km coastal mask, plus the six cities. The site selection
is derived from the current production mask rather than a historical cached mask.
These are VLM/GIA contributions to RSL, not total climate-driven sea-level projections.

Full samples, mean/std datasets at every decade, map PNGs, location files and run
manifests stay under ignored `data/diagnostics/shared_notebooks/`. Input hashes,
production-source hashes, helper hash, software versions and settings identify the
run; cached summaries are checksum-verified before reuse. The runner executes with
the current project interpreter. The locked optional `diagnostics` dependency group
provides plotting and Jupyter tools without adding them to the production image.
README and changelog link the new material. No production scientific source,
input archive, CLI, NetCDF schema, registry entry or release version changes in this
step; nothing is published to GitHub automatically.

Validation: both notebooks executed without errors and were checked visually.
Saved summaries are finite with nonnegative standard deviations, zero displacement
and spread at the base year, and the expected London coordinates. Means/stds at
ten selected sites (including every city) agree exactly with direct reductions of
the raw samples after float32 storage. Both modes are identical through 2050. The
47 existing tests pass, Ruff checks and the Just validation recipe pass, and
`git diff --check` is clean. Docker builds and its CLI help runs successfully;
the image still runs as UID 9876 and contains neither this diagnostics folder nor
the optional plotting/Jupyter packages. Full new global projection runs were made
on the host; prior host/container two-mode scientific equivalence is documented
in the preceding implementation entry.

## 2026-11-02: Optional convergence or constant secular rates

Added `--vlm-rate-mode` to the installed CLI and passed it through `cli.main`,
`vlm_postprocess`, and `projection_parameters`. The default, `converge-to-gia`,
retains the corrected cosine velocity transition from 2050 to 2150 and integrates
it to retain previously accumulated displacement. The alternative, `constant`,
uses elapsed time directly instead of the integrated taper. Thus the selected
spatial observation/GIA rate remains constant for the whole projection, and its
trend uncertainty continues accumulating. Unsupported mode values are rejected.

For constant mode the mean is `G*(t-b) + s*(t-b)*(U-V)`; for convergence it remains
`G*(t-b) + s*I(t,b)*(U-V)`. The VLM-trend and fallback-GIA uncertainty coefficients
use the same respective durations. GIA-only locations are identical in both modes.
Both modes preserve base-year referencing and fixed-seed/chunk reproducibility.

The switch applies to secular velocity and trend error only. The separate GRW
displacement-envelope model, calendar handling, and its taper remain unchanged;
constant mode does not invent a future GRW stochastic extrapolation. Input data,
CMR treatment, GIA bookkeeping, and the existing mean-rate spatial blend are unchanged.

The CLI addition is backward compatible: omitting it retains existing corrected
numerical outputs. An additive global NetCDF attribute, `vlm_rate_mode`, records the
choice. Existing variables, dimensions, units, fill values, and other metadata remain.
README help, Docker example, and changelog document the option. The separate module
registry is not edited: its existing commands continue to use the default. Maintainers
can expose the new choice through an optional string option sourced from
`module_inputs.options.vlm_rate_mode`, defaulting to `converge-to-gia`, at integration.

Regression coverage checks constant rates after 2150, accumulation of uncertainty,
GIA-only and partial-spatial-weight behavior, rebasing, default equivalence, invalid
mode rejection, CLI forwarding, output metadata, and deterministic chunking in both
modes. Local evidence under ignored `data/discovery/vlm_rate_modes_20261102/` keeps
the pre-switch source snapshot and six-city results for 500 samples, seed 1234,
base year 2000, and decadal years 2000–2300. The projection diagnostics notebook adds
a comparison section. The frozen original baseline and previous correction outputs
are preserved separately; no commit, push, registry update, or publication is made.

Validation: 47 tests pass, Ruff and the Just validation recipe pass, and the
development Docker image builds. Both modes were run locally and in Docker on the
six cities through 2300; host/container scientific arrays are exactly equal in both
modes. The default samples through 2150 are exactly equal to the previous candidate
outputs. A synthetic 2 mm/year non-GIA subsidence case with zero GIA reaches 200 mm
in 2150 and retains 200 mm through 2300 under convergence; constant mode reaches
300 mm in 2150 and 600 mm in 2300. The updated projection notebook executes without
errors. The new discovery outputs remain ignored, and `git diff --check` passes.

## 2026-11-02: Local diagnostics and candidate methodological correction

Status: implemented locally for comparison and scientific review; not committed,
published, or released. All pre-existing infrastructure and documentation changes
were preserved. The input datasets and CMR treatment are unchanged.

### Motivation and changes

- `vlm_preprocess.py`: replace the incorrect array `arctan` distance calculation
  with clipped `arctan2` and use a unit-sphere KD-tree for nearest-neighbor searches.
  GNSS is always queried against the reconstruction, regardless of relative dataset
  sizes. Normalize longitudes before applying the existing regional boxes.
- GNSS smoothing retains measurement uncertainty: with normalized weights `a`, its
  variance is `sum(a*(trend-mean)^2) + sum(a^2*sigma^2)`. Thus an isolated station
  retains its reported uncertainty. The scatter term represents unresolved spatial
  variation; independence of station measurement errors remains an approximation.
  Gaussian smoothing now searches its documented full three-scale radius.
- Normalize/sort longitude seams, interpolate by named dimensions, validate required
  fields/units and coordinates, reject negative standard deviations, and drop points
  missing a trend or its uncertainty with a logged count. Inputs are closed after
  loading needed variables; the example-file option remains accepted and validated.
- Preserve the static 0.25-degree coastal grid and existing region definitions, but
  retain the GRW array at source sites rather than expanding every source year onto
  the entire global grid. Extract source curves only for requested sites. If regional
  exclusions remove every observation, the calculation falls back to GIA.
- `vlm_postprocess.py`: compute a calendar-fixed cosine velocity taper from 2050 to
  2150 and integrate it analytically. If `s` is spatial weight, `I(t,b)` is the integral
  of the taper from base year `b`, `V` is observed upward VLM, `G` is GIA RSL, and `U`
  is GIA upward motion, mean displacement is
  `G*(t-b) + s*I(t,b)*(U-V)`. Accumulated non-GIA displacement is no longer erased
  in 2150, and rebasing is consistent with subtracting the same trajectory at `b`.
- VLM-trend uncertainty coefficient is `s*I*sigma_V`; fallback-GIA uncertainty
  coefficient is `((t-b)-s*I)*sigma_G`. They use independent seeded stratified normal
  streams instead of the same draw. The GRW envelope uses a third stream, the actual
  source calendar years, and subtracts its base-year envelope. It remains a marginal
  displacement-envelope approximation, not sampled GRW paths. Before the first
  finite GRW year, leading structural NaNs become zero; internal/trailing gaps and
  extrapolation during the active observational period are rejected.
- Retain the GIA sea-surface correction mean `G+U` without inventing missing joint
  covariance information. Shared spatial errors, GIA sea-surface uncertainty, and true
  GRW temporal covariance are limitations for subsequent scientific work.
- Write float32 samples directly in bounded sample/location blocks with an atomic
  final-file replacement. Preserve output names, dimensions/order, units, coordinate
  dtypes, fill values and global metadata; remove the intermediate int16 rounding
  and overflow. The output remains the local/RSL contribution and includes GIA.
- `cli.py` and `read_locationfile.py`: reject invalid parameters before preprocessing,
  use exact year steps, allow a single requested year, validate location coordinates
  and IDs with line-number errors, accept blank/comment lines, and prevent an output
  filename from overwriting an input. The chunk-size help now describes its existing
  sample-batch behavior; no option was removed or renamed.
- `.dockerignore`: exclude the ignored local agent guide from container images.

### Preserved baseline, local artifacts, and reproduction

The local workspace `data/discovery/vlm_review_20261102/` is ignored by Git and by the
Docker build. It contains exactly two new diagnostic notebooks, executed with their
plots embedded, a frozen original source copy, helper scripts, run logs and manifests,
input hashes, and separate baseline/candidate outputs. Exploratory data are not
part of the public package. The existing small methodology notebook contains pointers
and a synthetic verification example.

- `01_input_data_diagnostics.ipynb`: full variable inventory, VLM trend/uncertainty
  and footprint contours, input GRW spread, raw/corrected GNSS trends and uncertainty,
  GIA fields, coastline mask and spatial blending weights.
- `02_projection_diagnostics.ipynb`: original results first, candidate results and
  comparisons; city mean ± standard deviation, 2150 mean/std maps, and differences.
- Cities follow Figure 4 of Oelsmann et al. (2024), Nature Geoscience,
  <https://doi.org/10.1038/s41561-023-01357-2>: Venice, Skellefteå, Valparaíso, Tokyo,
  Brisbane and New York. The recorded representative city coordinates are not
  asserted to be the paper's exact tide-gauge coordinates. These plots show this
  module's VLM/GIA contribution, not the paper's combined climate projections.
- Both runs use the same five published input files, 500 samples, seed 1234, base
  year 2000, and years 2000–2150 every ten years. The six cities are accompanied by
  14,366 integer-degree coastal grid nodes within 150 km of a coast, selected once
  for both variants. They are not an area-weighted global statistic.
- `baseline_decadal_mean_std.nc` and `candidate_decadal_mean_std.nc` store actual
  finite-ensemble means and population standard deviations at every saved decade.
  Each is under 2 MB. Input/source SHA256 hashes and parameters are in the manifests;
  complete sample outputs remain available locally to audit rounding and statistics.

From the repository root, the reference command is
`uv run --locked python data/discovery/vlm_review_20261102/run_diagnostics.py baseline`;
the candidate command is identical except for its final `candidate` argument. The
baseline always imports its frozen source; the candidate imports the current package.
The local README also records notebook/environment reproduction commands.

### Quantitative comparison

The matched coastal-node RMS difference at 2150 is 66.60 mm for the mean and
50.88 mm for standard deviation. Largest absolute differences are 1409.44 mm and
619.97 mm respectively. These are sampled-node diagnostics, not global averages.

| City | Baseline mean ± std, 2150 (mm) | Candidate mean ± std, 2150 (mm) |
| --- | ---: | ---: |
| Venice | -14.06 ± 33.09 | 88.33 ± 43.37 |
| Skellefteå | -1230.49 ± 319.81 | -1342.39 ± 131.29 |
| Valparaíso | -81.61 ± 16.95 | -720.16 ± 218.04 |
| Tokyo | -64.35 ± 16.39 | -226.53 ± 140.76 |
| Brisbane | -40.84 ± 10.99 | 7.07 ± 13.04 |
| New York | 137.76 ± 55.83 | 187.62 ± 28.79 |

By 2150 the baseline is GIA-only. Most candidate differences then arise from retaining
the accumulated non-GIA mean and trend error and integrating a smaller historical
fallback-GIA uncertainty coefficient. Earlier GRW-envelope differences arise from
calendar alignment. Smaller spatial differences reflect corrected distances, longitude
seams, and GNSS smoothing; submillimetre differences reflect retaining float precision.
Separate uncertainty streams remove forced covariance but retain finite-sample
cross-correlations. No CMR or source-dataset difference explains these changes.

One local measurement gave baseline preprocessing/projection times of 249.87/7.40 s
and candidate times around 1/10 s. This is a single local comparison, with caching
and concurrent system activity, not a controlled performance benchmark. Projection
writing is slightly slower because it preserves float precision and streams output;
the large gain comes from preprocessing.

### Validation and release limitations

Scientific tests cover RSL signs and GIA cancellation, integrated taper/rebasing,
requested-range independence, exact year steps, actual GRW dates and missing data,
GNSS uncertainty, Gaussian support, spherical geometry and longitude conventions,
dataset sizes, GIA-only fallback, location validation, independent sampling,
fixed-seed chunk invariance, float precision and NetCDF schema, and a complete
synthetic CLI execution. No CI test requires the published archive.

The installed CLI, Ruff, pytest, and the Just validation recipe are checked locally.
The non-root development container is built and run on the six cities with the real
published inputs; its scientific contents are compared with the host candidate.
Detailed logs are in the ignored workspace. The system has no standalone `just`;
the same recipe is invoked with `uvx --from rust-just just validate`.

Final local verification: 39 tests passed; both notebooks executed without errors;
the Docker image built and ran with non-root execution; host and container scientific
arrays were exactly equal. A raw-file comparison confirmed the original NetCDF
variable set, dimension order, dtypes, fill values and attributes, except the expected
creation timestamp. Saved mean/std values were checked against the actual samples,
and source/input hashes match their manifests. `git diff --check` passed. None of the
discovery files or the agent guide is staged or included in the Docker image.

No commit, push, PR, registry update, container publication, version change, or release
is made. Maintainer review of the velocity-taper interpretation and remaining
covariance assumptions is still required before release. Registry help wording can
be coordinated then; the argument schema and production input archive are unchanged.

## 2026-09-14: Establish the local FACTS2 baseline

Status: infrastructure only; no scientific method changed.

Starting point:

- Branch: `correct-vlm-methodology`
- Commit: `1c4ab96fc29d9f1d174e99ed763a7eca02a72abc`
- Published input record: <https://doi.org/10.5281/zenodo.18199757>
- Archive MD5: `1de47917200232f9c68cf30759f88c24`

Changes:

- Corrected the installed command target to point directly to `cli.py`.
- Moved the Docker runtime user switch until after dependency installation so the
  image can create its virtual environment while still running the application as an
  unprivileged user.
- Added a CLI help smoke test and declared pytest as a development dependency.
- Added the test to the local validation recipe and pull-request workflow.
- Established an ignored local data layout for reference inputs, candidate inputs,
  outputs, and diagnostic artifacts.
- Added a diagnostic notebook template for dataset inventories and later numerical
  comparisons.

Compatibility:

- No CLI option, default, input schema, output schema, or numerical algorithm changed.
- The container entry point and non-root runtime are retained.

Validation evidence:

- `uv sync --locked --all-extras --dev` completed with Python 3.12.
- The installed `oelsmann24-verticallandmotion --help` command completed successfully.
- Before the Dockerfile repair, the image build failed while creating `.venv` with
  `Permission denied` under UID/GID 9876.
- After the repair, `oelsmann24-verticallandmotion:dev` built successfully, its help
  command completed, and the default runtime identity was confirmed as UID/GID
  9876:9876.
- Ruff formatting and lint checks passed.
- The CLI test passed with pytest.
- Two full container runs used all five published inputs, the New York test location,
  10 samples, years 2000--2100 at 10-year intervals, base year 2000, and seed 1234.
  Their scientific contents were equal. The byte-level NetCDF checksums differed only
  because each output records its creation time in the global `history` attribute.
- The reference output contained `sea_level_change(samples, years, locations)` with
  shape `(10, 11, 1)`, decoded dtype `float32`, units `mm`, and values from 0 to 157 mm.
- Both runs emitted an existing `RuntimeWarning` from `haversine` when evaluating
  `sqrt(1 - a)`. This warning is a baseline finding and was not corrected as part of
  the infrastructure repair.

Local generated evidence is ignored by Git and stored in:

- `data/diagnostics/summaries/reference_input_inventory.json`
- `data/diagnostics/summaries/baseline_output_summary.json`
- `data/outputs/reference/lslr.nc`
- `data/outputs/reference/lslr_repeat.nc`
