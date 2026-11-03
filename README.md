# oelsmann24-verticallandmotion

## Overview

This FACTS2 module uses data from the GlobalVLM workflow to generate vertical land
motion (VLM) contributions to local relative sea-level projections. It combines a
VLM reconstruction with GNSS estimates and glacial isostatic adjustment (GIA)
information, and writes samples for the locations specified by the user.

## Scientific background

The VLM reconstruction used as an input to this module is described in:

Oelsmann, J., Marcos, M., Passaro, M., Sanchez, L., Dettmering, D., Dangendorf, S.,
and Seitz, F. (2024). Regional variations in relative sea-level changes influenced
by nonlinear vertical land motion. *Nature Geoscience*, **17**, 137–144.
[https://doi.org/10.1038/s41561-023-01357-2](https://doi.org/10.1038/s41561-023-01357-2).

The paper provides the scientific background for the source reconstruction. This
repository documents its use in a FACTS2 module; the example below demonstrates a
module run for sea-level projections.

## Code reference

The module source is available in the
[oelsmann24-verticallandmotion repository](https://github.com/fact-sealevel/oelsmann24-verticallandmotion).
The published software release is
[v0.1.0](https://github.com/fact-sealevel/oelsmann24-verticallandmotion/releases/tag/v0.1.0).

## Data reference

### Input data

Oelsmann, J., and Malagón-Santos, V. (2026). *Vertical land motion data for FACTS*
(v1) [Data set]. Zenodo.
[https://doi.org/10.5281/zenodo.18199757](https://doi.org/10.5281/zenodo.18199757).

The bundled inputs and their sources are listed below, following the source
information in the [Zenodo record](https://zenodo.org/records/18199757).

| Input file | Source |
| --- | --- |
| `VLM_reconstruction.nc` | Oelsmann et al. [reconstruction data supplement](https://doi.org/10.5281/zenodo.8308347), associated with the paper cited above. |
| `NGL14_CMR_corrected.nc` | [Nevada Geodetic Laboratory](http://geodesy.unr.edu/); [MIDAS, Blewitt et al. (2016)](https://doi.org/10.1002/2015JB012552); [Blewitt et al. (2018)](https://doi.org/10.1029/2018EO104623). Contemporary mass redistribution (CMR) correction: [Frederikse et al. (2020)](https://doi.org/10.1038/s41586-020-2591-3). |
| `GIA_stats.nc` | [Caron et al. (2018)](https://doi.org/10.1002/2017GL076644), also distributed in the [Frederikse, Caron, and Landerer (2019) archive](https://doi.org/10.5281/zenodo.3485577). |
| `dist2coast_1deg_v2.grd` | Coastal-distance grid produced by NASA's Ocean Biology Processing Group (OBPG) with Generic Mapping Tools (GMT), distributed through PacIOOS. |
| `ssp585.2150.fair2.emuGLA.emulandice2.glaciers_quantiles.nc` | FACTS output sample included in the input archive. |

### Output data

The example writes `data/output/lslr.nc`, containing
`sea_level_change(samples, years, locations)` in millimetres relative to the chosen
base year, together with the site coordinates.

## Contributing modeling software

| Model | Version | Repository |
| --- | --- | --- |
| `oelsmann24-verticallandmotion` | `0.1.0` | [Source code](https://github.com/fact-sealevel/oelsmann24-verticallandmotion) |

## Reproduce a module run

Download and organize input data:
> [!WARNING]
> The following command downloads approximately 1.3 GB of compressed data. Extraction requires additional disk space.

```shell

#Download input data
mkdir -p ./data/input
curl -sL https://zenodo.org/records/18199757/files/oelsmann24_vlm_data.tar.gz | tar -zx -C ./data/input
mv ./data/input/data/* ./data/input/ && rmdir ./data/input/data

printf 'New_York\t12\t40.70\t-74.01\n' > ./data/input/location.lst

#Make dir for output data 
mkdir -p ./data/output
```

Build the development image from this checkout, then run it. The new options below
are not yet available in the published release image.

```shell
docker build -t oelsmann24-verticallandmotion:dev .

docker run --rm \
--user "$(id -u):$(id -g)" \
-v "$(pwd)/data/input:/mnt/oelsmann24_in:ro" \
-v "$(pwd)/data/output:/mnt/oelsmann24_out" \
oelsmann24-verticallandmotion:dev \
--pipeline-id "1234" \
--nsamps 500 \
--input-vlm-reconstruction-file "/mnt/oelsmann24_in/VLM_reconstruction.nc" \
--input-gps-vlm-file "/mnt/oelsmann24_in/NGL14_CMR_corrected.nc" \
--input-gia-vlm-file "/mnt/oelsmann24_in/GIA_stats.nc" \
--input-dist2coast-file "/mnt/oelsmann24_in/dist2coast_1deg_v2.grd" \
--input-example-file "/mnt/oelsmann24_in/ssp585.2150.fair2.emuGLA.emulandice2.glaciers_quantiles.nc" \
--seed 1234 \
--baseyear 2000 \
--pyear-start 2000 \
--pyear-step 10 \
--pyear-end 2100 \
--location-file "/mnt/oelsmann24_in/location.lst" \
--chunk-size 50 \
--vlm-rate-mode "converge-to-gia" \
--output-lslr-file "/mnt/oelsmann24_out/lslr.nc"
```

## Command-line options
Several options and configurations are available when running the container.
```shell
Usage: oelsmann24-verticallandmotion [OPTIONS]

  Command-line interface for the oelsmann24-verticallandmotion package.

Options:
  --pipeline-id TEXT              Unique identifier for the processing pipeline.
  --input-vlm-reconstruction-file TEXT
                                  Path to the input VLM reconstruction data
                                  file. Required unless --gia-only.
  --input-gps-vlm-file TEXT       Path to the input GPS VLM data file. Required
                                  unless --gia-only.
  --input-gia-vlm-file TEXT       Path to the input GIA VLM data file.
                                  [required]
  --input-dist2coast-file TEXT    Path to the input distance to coast GRD file.
                                  Required unless --gia-only.
  --input-example-file TEXT       Path to the input example data file. Required
                                  unless --gia-only.
  --nsamps INTEGER                Number of samples to generate.  [required]
  --seed INTEGER                  Seed value for random number generator.
                                  [default: 1234]
  --baseyear INTEGER              Base or reference year for projections
  --pyear-start INTEGER           Year for which projections start  [default:
                                  2000]
  --pyear-step INTEGER            Step sizes in years between pyear-start and
                                  pyear-stop  [default: 10]
  --pyear-end INTEGER             Year for which projections end  [default:
                                  2100]
  --location-file TEXT            File that contains name, id, lat, and lon of
                                  points for localization  [required]
  --chunk-size INTEGER            Number of samples per computation block
                                  (locations are streamed separately)  [default:
                                  50]
  --vlm-rate-mode [converge-to-gia|constant]
                                  Converge secular rates to GIA during
                                  2050-2150, or keep the observed rates
                                  constant. GRW variability is unchanged.
                                  [default: converge-to-gia]
  --gia-only                      Project only GIA RSL and its uncertainty; skip
                                  observational inputs and ignore --vlm-rate-
                                  mode.
  --output-lslr-file TEXT         Path to the output local sea-level rise
                                  projections file.  [required]
  --help                          Show this message and exit.
```

See this help documentation by running:

```shell
docker run --rm oelsmann24-verticallandmotion:dev --help
```

## Building the container locally

You can build the container with Docker by cloning the repository locally and then running

```shell
docker build -t oelsmann24-verticallandmotion:dev .
```

from the repository root.

### Local data and diagnostic workspace

Keep large inputs and generated results under the ignored `data/` directory. A useful
development layout is:

```text
data/
├── downloads/             # Downloaded archives
├── inputs/
│   ├── reference/         # Published reference inputs
│   └── candidate/         # Replacement inputs under evaluation
├── shared/                # Location files and other small run inputs
├── outputs/
│   ├── reference/
│   └── candidate/
└── diagnostics/           # Generated figures, tables, and summaries
```

The container does not contain these data. Bind-mount the selected host input directory
read-only and mount a dedicated output directory as writable. For example, adapt the
earlier run command to use the locally built image and these mounts:

```shell
-v "$(pwd)/data/inputs/reference:/mnt/oelsmann24_in:ro" \
-v "$(pwd)/data/shared:/mnt/shared_in:ro" \
-v "$(pwd)/data/outputs/reference:/mnt/oelsmann24_out" \
oelsmann24-verticallandmotion:dev \
```

Use `/mnt/shared_in/location.lst` for `--location-file`. On Linux, if the fixed non-root
container user cannot write to the host output directory, add
`--user "$(id -u):$(id -g)"` before the mount arguments. This still runs the application
without root privileges.

The versioned [diagnostics folder](diagnostics/README.md) contains two executed
notebooks: [input/source data](diagnostics/01_input_data.ipynb) and
[corrected projections](diagnostics/02_projection_comparison.ipynb). The latter
compares convergence to GIA with constant secular rates for Venice, Skellefteå,
London, Tokyo, Brisbane and New York, and maps both modes' mean and standard
deviation in 2150 over a complete global 1° grid, including open-ocean GIA RSL.
Supporting Python source imports the production model rather than duplicating its
scientific logic. Reproduce both notebooks with:

```shell
uv sync --locked --group diagnostics --dev
uv run --locked --group diagnostics python diagnostics/execute_notebooks.py
```

Large inputs and generated NetCDFs remain ignored under `data/`; only the notebook
source, embedded plots and helpers are versioned. The optional `diagnostics`
dependencies are not installed in the production Docker image. The existing
`notebooks/methodology_diagnostics.ipynb` remains a lightweight methodology summary.
Other exploratory notebooks, historical baseline comparisons, source snapshots
and generated data stay local under ignored `data/discovery/`. The historical
workspace is described in `data/discovery/vlm_review_20261102/README.md` when available.

### Projection method

#### Combining observations and GIA

The model includes the relative sea-level (RSL) contribution from GIA everywhere.
In the open ocean, and where coastal observations are unavailable, it uses GIA only.

At the coast, we use land-motion observations from GNSS and the Oelsmann et al.
(2024) reconstruction to adjust the GIA estimate. We first subtract the vertical
land motion predicted by GIA from the observed land motion. The difference is the
land motion not explained by GIA. We use this difference to adjust the GIA RSL
prediction: extra subsidence increases relative sea level, while extra uplift
reduces it.

A coastal mask limits this correction to locations within 150 km of the coast.
Within that band, a cosine blend reduces its spatial weight with distance from
the nearest accepted observation: full weight at the observation, half weight
at 50 km, and zero weight at 100 km or farther away. Outside these limits, the
model uses GIA only. This keeps the GIA RSL contribution everywhere without
counting GIA land motion twice. Do not add a second complete GIA projection to
this output.

#### Projection options

Choose how the secular VLM rates evolve by passing either option to the local or
Docker command:

- `--vlm-rate-mode constant`: retain the observation-based secular rates throughout
  the projection, with the spatial observation/GIA blend unchanged. Trend
  uncertainty continues to accumulate with elapsed time.
- `--vlm-rate-mode converge-to-gia` (default): gradually bring secular rates toward
  GIA. The temporal observation weight stays at one through 2050, decreases with a
  cosine taper to zero in 2150, and remains zero thereafter.

The temporal taper acts on velocities, which are integrated to obtain displacement.
Displacement accumulated before 2150 is therefore retained afterward. Both modes
reference displacement to the selected base year. The option controls secular
rates and trend uncertainty; the separate random-walk displacement envelope keeps
its taper to zero in 2150 in both modes.

#### GIA-only projections

Add `--gia-only` to produce only the GIA contribution to RSL and its uncertainty
at every requested location. The GIA RSL rate (`rsl_mean`) and its standard error
(`rsl_sterr`) accumulate linearly from the selected base year. Observed VLM and
random-walk variability are excluded, and `--vlm-rate-mode` has no effect.

Only the GIA data and location file are needed; reconstruction, GNSS, distance-to-coast
and example inputs can be omitted. From the local checkout, run:

```shell
uv run --locked oelsmann24-verticallandmotion \
  --gia-only \
  --input-gia-vlm-file ./data/input/GIA_stats.nc \
  --location-file ./data/input/location.lst \
  --nsamps 500 \
  --seed 1234 \
  --baseyear 2000 \
  --pyear-start 2000 \
  --pyear-end 2100 \
  --pyear-step 10 \
  --output-lslr-file ./data/output/gia_only.nc
```

For Docker, use the same switch with the locally built development image and the
mounted input/output paths. The output keeps the usual samples, years and locations
and records `vlm_rate_mode="gia-only"`. Omitting the switch retains the existing
observation/GIA projection behavior.

#### Uncertainty and sampling

Projection uncertainty combines observation-based VLM trend errors, Gaussian
random-walk displacement spreads from the Oelsmann et al. (2024) reconstruction
representing temporal VLM variability, and GIA RSL errors for the fallback GIA
contribution. These three error components are perturbed using independent seeded
sampling streams; their joint covariance is not modeled. The GIA sea-surface
correction is applied as a fixed mean, without sampling its uncertainty.

The module writes samples; quantiles can be computed downstream. Detailed
method changes and assumptions are recorded in the
[development notes](docs/development_changes.md).

## Support

Source code is available online at https://github.com/fact-sealevel/oelsmann24-verticallandmotion. This software is open source, available under the MIT license.

Please file issues in the issue tracker at https://github.com/fact-sealevel/oelsmann24-verticallandmotion/issues.
