"""Project the observation/GIA contribution to RSL using integrated velocity weights."""

import os
from pathlib import Path
import tempfile
import time

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial import cKDTree
from scipy.stats import norm

from .read_locationfile import ReadLocationFile
from .vlm_preprocess import _to_xyz, _arc_km_from_chord

VLM_RATE_MODES = ("converge-to-gia", "constant")


def match_coord_vectors(site_lats, site_lons, lat_flat, lon_flat):
    """Nearest grid node on a sphere, accepting either longitude convention."""
    if len(lat_flat) == 0:
        raise ValueError("Cannot map onto an empty grid")
    tree = cKDTree(_to_xyz(lon_flat, lat_flat))
    chord, idx = tree.query(_to_xyz(site_lons, site_lats))
    return idx, _arc_km_from_chord(chord)


def projection_years(baseyear, start, end, step):
    """Exact integer step; end is included only when it falls on that step."""
    if any(
        isinstance(v, bool) or not isinstance(v, (int, np.integer))
        for v in [baseyear, start, end, step]
    ):
        raise ValueError("Projection years and step must be integers")
    if not 2000 <= baseyear <= 2300 or not 2000 <= start <= end <= 2300:
        raise ValueError(
            "Base and projection years must be in 2000..2300; start <= end"
        )
    if step < 1:
        raise ValueError("Projection year step must be positive")
    return np.arange(start, end + 1, step, dtype=np.int64)


def temp_decay(years, start_year=None, period1=2050, period2=2150):
    """Calendar-fixed cosine velocity weight, independent of requested year range."""
    if period2 <= period1:
        raise ValueError("Decay end must follow decay start")
    years = np.asarray(years, dtype=float)
    offset = np.clip(years - period1, 0, period2 - period1)
    weight = (1 + np.cos(np.pi * offset / (period2 - period1))) / 2
    inside = (years >= period1) & (years < period2)
    return weight[inside], years[inside] - period1, weight


def integrated_decay(years, baseyear, period1=2050, period2=2150):
    """Exact integral of the cosine velocity taper from baseyear to each year."""
    if period2 <= period1:
        raise ValueError("Decay end must follow decay start")
    duration = period2 - period1

    def primitive(y):
        y = np.asarray(y, dtype=float)
        q = np.clip(y - period1, 0, duration)
        return (
            np.minimum(y, period1)
            + q / 2
            + duration / (2 * np.pi) * np.sin(np.pi * q / duration)
        )

    return primitive(years) - primitive(baseyear)


def _grw_at_years(grw, source_ids, years):
    """Use actual years; leading pre-projection NaNs are structural zeros.

    Internal/trailing missing values are rejected. No extrapolation is allowed
    while the observational contribution is active; outside it the result is zero.
    """
    if not np.any(source_ids >= 0):
        return np.zeros((len(source_ids), len(years)))
    if set(grw.dims) != {"x", "time"}:
        raise ValueError("GRW_ensemble_std must have x and time dimensions")
    t = grw.time
    source_years = (
        t.dt.year.values
        if np.issubdtype(t.dtype, np.datetime64) or t.dtype.kind == "O"
        else np.asarray(t, dtype=float)
    )
    if (
        len(source_years) < 2
        or not np.all(np.isfinite(source_years))
        or np.any(np.diff(source_years) <= 0)
        or np.min(source_years) < 1800
    ):
        raise ValueError(
            "GRW time must contain increasing, unique absolute calendar years"
        )
    valid_sites = source_ids >= 0
    unique, inverse = np.unique(source_ids[valid_sites], return_inverse=True)
    values = grw.transpose("x", "time").isel(x=unique).values.astype(float)
    if np.any(values[np.isfinite(values)] < 0):
        raise ValueError("GRW standard deviations cannot be negative")
    for row in values:
        finite = np.flatnonzero(np.isfinite(row))
        if not len(finite):
            raise ValueError("GRW source has no finite standard deviations")
        row[: finite[0]] = 0
        if not np.all(np.isfinite(row)):
            raise ValueError("GRW has internal or trailing missing values")
    years = np.asarray(years, dtype=float)
    active = temp_decay(years)[2] > 0
    if np.any((years[active] < source_years[0]) | (years[active] > source_years[-1])):
        raise ValueError("Requested active GRW years exceed source time coverage")
    result = np.zeros((len(source_ids), len(years)))
    if np.any(active):
        f = interp1d(source_years, values, axis=1, bounds_error=True)
        result[np.ix_(valid_sites, active)] = f(years[active])[inverse]
    return result * 1000  # Input displacement spread in m; output mm.


def projection_parameters(
    preprocess_dict,
    site_lats,
    site_lons,
    years,
    baseyear,
    rate_mode="converge-to-gia",
    gia_only=False,
):
    """Return site/year mean and three independent uncertainty coefficients.

    GIA sea-surface mean is retained in the observational branch. Its covariance
    with GIA radial deformation is unavailable; no new uncertainty model for that
    correction is invented here. Spatial and GRW path covariance remain simplified.
    rate_mode controls secular velocity and trend uncertainty only; the existing
    GRW displacement-envelope taper is unchanged in both modes.
    gia_only bypasses observational fields and uses the full GIA RSL rate/error.
    """
    if rate_mode not in VLM_RATE_MODES:
        raise ValueError(f"rate_mode must be one of {VLM_RATE_MODES}")
    gia = preprocess_dict["GIA_VLM_int"]
    fields = gia if gia_only else preprocess_dict["VLM_REC_MERGED"]
    xx, yy = np.meshgrid(fields.lon, fields.lat)
    idx, _ = match_coord_vectors(site_lats, site_lons, yy.ravel(), xx.ravel())

    def take(da):
        return da.transpose("lat", "lon").values.ravel()[idx]

    dt = np.asarray(years, dtype=float) - baseyear
    if gia_only:
        mean = dt[:, None] * take(gia.rsl_mean)
        gia_sigma = dt[:, None] * take(gia.rsl_sterr)
        zeros = np.zeros_like(mean)
        return mean, zeros, zeros.copy(), gia_sigma

    spatial_weight = take(preprocess_dict["weights"])
    trend = take(fields.VLM_trend_coefficient_mean) * 1000
    sigma = take(fields.VLM_trend_coefficient_uncertainty) * 1000
    if np.any(
        (spatial_weight > 0) & (~np.isfinite(trend) | ~np.isfinite(sigma) | (sigma < 0))
    ):
        raise ValueError("Observation trend or uncertainty missing at a weighted site")
    trend = np.where(spatial_weight > 0, trend, 0)
    sigma = np.where(spatial_weight > 0, sigma, 0)
    duration = (
        integrated_decay(years, baseyear) if rate_mode == "converge-to-gia" else dt
    )
    obs_duration = duration[:, None] * spatial_weight[None, :]
    gia_duration = dt[:, None] - obs_duration
    mean = obs_duration * (-trend + take(gia.asl)) + gia_duration * take(gia.rsl_mean)
    obs_sigma = obs_duration * sigma
    gia_sigma = gia_duration * take(gia.rsl_sterr)
    ids = take(fields.source_index).astype(int)
    eval_years = np.concatenate([np.asarray(years), [baseyear]])
    grw = _grw_at_years(preprocess_dict["VLM_REC_MERGED_GRW"], ids, eval_years)
    # This displacement envelope, not a velocity, is referenced consistently to b.
    envelope = grw * temp_decay(eval_years)[2][None, :]
    grw_sigma = spatial_weight[None, :] * (envelope[:, :-1] - envelope[:, -1:]).T
    return mean, obs_sigma, grw_sigma, gia_sigma


def sample_draws(nsamps, seed):
    """Independent stratified Gaussian streams, stable across chunk sizes."""
    if nsamps < 1 or seed < 0:
        raise ValueError("nsamps must be positive and seed nonnegative")
    q = norm.ppf(np.arange(1, nsamps + 1) / (nsamps + 1))
    return tuple(np.random.default_rng(seed + i).permutation(q) for i in range(3))


def vlm_postprocess(
    preprocess_dict: dict,
    nsamps,
    rng_seed,
    location_file,
    baseyear,
    pyear_start,
    pyear_end,
    pyear_step,
    chunksize,
    pipeline_id,
    output_lslr_file: str,
    rate_mode="converge-to-gia",
    gia_only=False,
):
    """Stream floating-point RSL samples to the established FACTS NetCDF schema."""
    years = projection_years(baseyear, pyear_start, pyear_end, pyear_step)
    if chunksize < 1:
        raise ValueError("chunk-size must be positive")
    draws = sample_draws(nsamps, rng_seed)
    _, ids, lats, lons = ReadLocationFile(location_file)
    coefficients = projection_parameters(
        preprocess_dict,
        lats,
        lons,
        years,
        baseyear,
        rate_mode=rate_mode,
        gia_only=gia_only,
    )
    if not all(np.all(np.isfinite(a)) for a in coefficients):
        raise ValueError("Projection contains missing or nonfinite coefficients")
    output = Path(output_lslr_file)
    location_block = min(len(ids), 1024)
    # Bound arithmetic memory even when the requested sample batch is very large.
    sample_block = min(
        chunksize, nsamps, max(1, 64 * 1024**2 // (len(years) * location_block * 8 * 5))
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".vlm-", suffix=".nc", dir=output.parent)
    os.close(fd)
    try:
        with Dataset(temporary, "w") as nc:
            for name, size in [
                ("samples", nsamps),
                ("years", len(years)),
                ("locations", len(ids)),
            ]:
                nc.createDimension(name, size)
            for name, data, dtype, dim in [
                ("samples", np.arange(nsamps), "i8", "samples"),
                ("years", years, "f8", "years"),
                ("locations", ids, "i8", "locations"),
                ("lat", lats, "f8", "locations"),
                ("lon", lons, "f8", "locations"),
            ]:
                options = {"fill_value": np.nan} if dtype == "f8" else {}
                nc.createVariable(name, dtype, (dim,), **options)[:] = data
            nc.setncatts(
                {
                    "description": (
                        "GIA contributions to RSL change"
                        if gia_only
                        else "Global VLM contributions to RSL change according to Global vlm workflow"
                    ),
                    "history": "Created " + time.ctime(),
                    "source": "SLR Framework: Global vlm workflow",
                    "scenario": "NA",
                    "baseyear": baseyear,
                    "vlm_rate_mode": "gia-only" if gia_only else rate_mode,
                }
            )
            var = nc.createVariable(
                "sea_level_change",
                "f4",
                ("samples", "years", "locations"),
                zlib=True,
                complevel=4,
                fill_value=np.nan,
                chunksizes=(sample_block, len(years), location_block),
            )
            var.setncatts({"units": "mm", "missing_value": np.float32(np.nan)})
            mean, obs, grw, gia = coefficients
            for loc in range(0, len(ids), location_block):
                sl = slice(loc, loc + location_block)
                for start in range(0, nsamps, sample_block):
                    stop = min(nsamps, start + sample_block)
                    block = mean[None, :, sl] + sum(
                        coefficient[None, :, sl] * draw[start:stop, None, None]
                        for coefficient, draw in zip([obs, grw, gia], draws)
                    )
                    var[start:stop, :, sl] = block.astype(np.float32)
        os.replace(temporary, output)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
