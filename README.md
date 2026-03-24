# oelsmann24-verticallandmotion

This module contains the vertical land motion module from the GlobalVLM workflow.

> [!CAUTION]
> This is a prototype. It is likely to change in breaking ways. It might delete all your data. Don't use it in production.

## Example

Download and organize input data:
> [!WARNING]
> The following command will download 8 GB of data.

```shell

#Download input data
mkdir -p ./data/input
curl -sL https://zenodo.org/records/18199757/files/oelsmann24_vlm_data.tar.gz | tar -zx -C ./data/input
mv ./data/input/data/* ./data/input/ && rmdir ./data/input/data

echo "New_York	12	40.70	-74.01" > ./data/input/location.lst

#Make dir for output data 
mkdir -p ./data/output
```

Run the application in a container, such as: 
```shell
docker run --rm \
-v ./data/input:/mnt/oelsmann24_in:ro \
-v ./data/output:/mnt/oelsmann24_out \
ghcr.io/fact-sealevel/oelsmann24-verticallandmotion:latest \
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
--output-lslr-file "/mnt/oelsmann24_out/lslr.nc"
```

## Features
Several options and configurations are available when running the container.
```shell
Usage: oelsmann24-verticallandmotion [OPTIONS]

  Command-line interface for the oelsmann24-verticallandmotion package.

Options:
  --pipeline-id TEXT              Unique identifier for the processing
                                  pipeline.
  --input-vlm-reconstruction-file TEXT
                                  Path to the input VLM reconstruction data
                                  file.  [required]
  --input-gps-vlm-file TEXT       Path to the input GPS VLM data file.
                                  [required]
  --input-gia-vlm-file TEXT       Path to the input GIA VLM data file.
                                  [required]
  --input-dist2coast-file TEXT    Path to the input distance to coast GRD
                                  file.  [required]
  --input-example-file TEXT       Path to the input example data file.
                                  [required]
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
  --chunk-size INTEGER            Number of locations to process at a time
                                  [default: 50]
  --output-lslr-file TEXT         Path to the output local sea-level rise
                                  projections file.  [required]
  --help                          Show this message and exit.
```

See this help documentation by running:

```shell
docker run --rm ghcr.io/fact-sealevel/oelsmann24-verticallandmotion:latest --help
```

## Building the container locally

You can build the container with Docker by cloning the repository locally and then running

```shell
docker build -t oelsmann24-verticallandmotion .
```

from the repository root.

## Support

Source code is available online at https://github.com/fact-sealevel/oelsmann24-verticallandmotion. This software is open source, available under the MIT license.

Please file issues in the issue tracker at https://github.com/fact-sealevel/oelsmann24-verticallandmotion/issues.
