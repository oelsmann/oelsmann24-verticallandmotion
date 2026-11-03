import click
import logging

from oelsmann24_verticallandmotion.vlm_preprocess import (
    global_preprocess_verticallandmotion,
    preprocess_gia_only,
)
from oelsmann24_verticallandmotion.vlm_postprocess import (
    vlm_postprocess,
    projection_years,
    VLM_RATE_MODES,
)
from oelsmann24_verticallandmotion.read_locationfile import ReadLocationFile
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "--pipeline-id",
    type=str,
    default="1234",
    help="Unique identifier for the processing pipeline.",
)
@click.option(
    "--input-vlm-reconstruction-file",
    type=str,
    help="Path to the input VLM reconstruction data file. Required unless --gia-only.",
)
@click.option(
    "--input-gps-vlm-file",
    type=str,
    help="Path to the input GPS VLM data file. Required unless --gia-only.",
)
@click.option(
    "--input-gia-vlm-file",
    type=str,
    required=True,
    help="Path to the input GIA VLM data file.",
)
@click.option(
    "--input-dist2coast-file",
    type=str,
    help="Path to the input distance to coast GRD file. Required unless --gia-only.",
)
@click.option(
    "--input-example-file",
    type=str,
    help="Path to the input example data file. Required unless --gia-only.",
)
@click.option(
    "--nsamps",
    type=int,
    required=True,
    help="Number of samples to generate.",
)
@click.option(
    "--seed",
    type=int,
    default=1234,
    show_default=True,
    help="Seed value for random number generator.",
)
@click.option(
    "--baseyear",
    type=int,
    default=2005,
    help="Base or reference year for projections",
)
@click.option(
    "--pyear-start",
    type=int,
    default=2000,
    help="Year for which projections start",
    show_default=True,
)
@click.option(
    "--pyear-step",
    default=10,
    type=int,
    show_default=True,
    help="Step sizes in years between pyear-start and pyear-stop",
)
@click.option(
    "--pyear-end",
    type=int,
    default=2100,
    show_default=True,
    help="Year for which projections end",
)
@click.option(
    "--location-file",
    type=str,
    help="File that contains name, id, lat, and lon of points for localization",
    required=True,
)
@click.option(
    "--chunk-size",
    type=int,
    default=50,
    show_default=True,
    help="Number of samples per computation block (locations are streamed separately)",
)
@click.option(
    "--vlm-rate-mode",
    type=click.Choice(VLM_RATE_MODES),
    default="converge-to-gia",
    show_default=True,
    help="Converge secular rates to GIA during 2050-2150, or keep the observed rates constant. GRW variability is unchanged.",
)
@click.option(
    "--gia-only",
    is_flag=True,
    help="Project only GIA RSL and its uncertainty; skip observational inputs and ignore --vlm-rate-mode.",
)
@click.option(
    "--output-lslr-file",
    type=str,
    required=True,
    help="Path to the output local sea-level rise projections file.",
)
def main(
    pipeline_id: str,
    input_vlm_reconstruction_file: str | None,
    input_gps_vlm_file: str | None,
    input_gia_vlm_file: str,
    input_dist2coast_file: str | None,
    input_example_file: str | None,
    nsamps: int,
    seed: int,
    baseyear: int,
    pyear_start: int,
    pyear_step: int,
    pyear_end: int,
    location_file: str,
    chunk_size: int,
    output_lslr_file: str,
    vlm_rate_mode: str,
    gia_only: bool,
):
    """Command-line interface for the oelsmann24-verticallandmotion package."""
    try:
        projection_years(baseyear, pyear_start, pyear_end, pyear_step)
        if nsamps < 1 or chunk_size < 1 or seed < 0:
            raise ValueError(
                "nsamps and chunk-size must be positive; seed must be nonnegative"
            )
        observation_inputs = {
            "--input-vlm-reconstruction-file": input_vlm_reconstruction_file,
            "--input-gps-vlm-file": input_gps_vlm_file,
            "--input-dist2coast-file": input_dist2coast_file,
            "--input-example-file": input_example_file,
        }
        inputs = [input_gia_vlm_file, location_file]
        if not gia_only:
            for flag, path in observation_inputs.items():
                if path is None:
                    raise ValueError(
                        f"Missing option '{flag}' (required unless --gia-only)"
                    )
                inputs.append(path)
        for path in inputs:
            if not Path(path).is_file():
                raise ValueError(f"Input file does not exist: {path}")
        protected_inputs = inputs + [p for p in observation_inputs.values() if p]
        if Path(output_lslr_file).resolve() in [
            Path(p).resolve() for p in protected_inputs
        ]:
            raise ValueError("Output path must not overwrite an input file")
        ReadLocationFile(location_file)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc
    logger.info("Starting preprocessing...")
    if gia_only:
        preprocess_dict = preprocess_gia_only(input_gia_vlm_file)
    else:
        preprocess_dict = global_preprocess_verticallandmotion(
            pipeline_id=pipeline_id,
            input_vlm_reconstruction_file=input_vlm_reconstruction_file,
            input_gps_vlm_file=input_gps_vlm_file,
            input_gia_vlm_file=input_gia_vlm_file,
            input_dist2coast_file=input_dist2coast_file,
            input_example_file=input_example_file,
        )
    logger.info("Preprocessing completed.")
    logger.info("Starting postprocessing...")

    vlm_postprocess(
        preprocess_dict=preprocess_dict,
        nsamps=nsamps,
        rng_seed=seed,
        location_file=location_file,
        baseyear=baseyear,
        pyear_start=pyear_start,
        pyear_end=pyear_end,
        pyear_step=pyear_step,
        chunksize=chunk_size,
        pipeline_id=pipeline_id,
        output_lslr_file=output_lslr_file,
        rate_mode=vlm_rate_mode,
        gia_only=gia_only,
    )
    logger.info("Postprocessing completed.")
