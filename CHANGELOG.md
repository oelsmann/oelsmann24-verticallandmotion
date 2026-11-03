# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Optional `--gia-only` switch to project GIA RSL means and uncertainties without
  observational inputs. The usual output schema and GIA sampling stream are retained.
- Versioned input-data and corrected-projection notebooks in `diagnostics/`, with
  embedded plots and reproducible Python helpers. Compare both rate modes on
  full global maps (including open-ocean GIA RSL) and at six cities, including
  London instead of Valparaíso; retain large data locally and keep plotting/Jupyter
  dependencies in an optional group.
- Optional `--vlm-rate-mode converge-to-gia|constant` switch for converging secular
  rates to GIA or continuing the selected rates. Convergence remains the default;
  the chosen mode is recorded in output metadata. GRW treatment is unchanged.
- CLI help smoke test and a reproducible methodology-diagnostics notebook.
- Plain-language development and methodology change record.
- Scientific, input-validation, deterministic-sampling, NetCDF-schema, and synthetic
  end-to-end regression tests with no large-data downloads.

### Fixed

- Correct the console-script target to reference `oelsmann24_verticallandmotion.cli`.
- Install container dependencies before switching to the non-root runtime user.
- Correct spherical distances, longitude seam handling, GNSS matching when dataset
  sizes change, and propagation of GNSS measurement uncertainty during smoothing.
- Use exact requested year steps and calendar-based decay/GRW interpolation.
- Integrate the observational velocity taper and reference the resulting trajectory
  to the chosen base year, retaining accumulated non-GIA displacement.
- Separate the VLM-trend and fallback-GIA uncertainty sampling streams.
- Preserve floating-point precision and stream sample blocks to NetCDF, avoiding
  integer overflow and the full global space-by-time preprocessing cube.
- Validate parameters, locations, input units and dimensions before computation;
  protect input files from accidental output overwrites.

These numerical changes are a candidate methodological correction and require
scientific review before release. The source input datasets and CMR treatment are
unchanged; CLI option names and the FACTS local/RSL NetCDF schema are retained.


## [0.1.0] - 2026-03-24

-Initial release

[Unreleased]: https://github.com/fact-sealevel/oelsmann24-verticallandmotion/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/fact-sealevel/oelsmann24-verticallandmotion/releases/tag/v0.1.0
