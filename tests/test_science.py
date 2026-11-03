import numpy as np
import pytest
import xarray as xr

from oelsmann24_verticallandmotion.vlm_preprocess import (
    EARTH_RADIUS_KM,
    haversine,
    compute_dist_between_datasets,
    radial_smoothing,
    merge_gps_and_vlm_reconstruction,
    periodic_interp,
    flatten_lon_lat,
)
from oelsmann24_verticallandmotion.vlm_postprocess import (
    temp_decay,
    integrated_decay,
    projection_years,
    projection_parameters,
    sample_draws,
    vlm_postprocess,
    _grw_at_years,
)
from oelsmann24_verticallandmotion.read_locationfile import ReadLocationFile


def point_data(lons, lats):
    return xr.Dataset(
        coords={"x": np.arange(len(lons)), "lon": ("x", lons), "lat": ("x", lats)}
    )


@pytest.mark.parametrize("lat,lon", [(0, 0), (0, 90), (0, 180), (90, 0)])
def test_haversine_scalar_array(lat, lon):
    scalar = haversine((0.0, 0.0), (float(lat), float(lon)))
    array = haversine((0.0, 0.0), (np.array([lat]), np.array([lon])))[0]
    assert np.isfinite(array)
    assert scalar == pytest.approx(array)
    if lon == 180:
        assert array == pytest.approx(np.pi * EARTH_RADIUS_KM * 1000)


def test_spherical_matching_seam_and_bruteforce():
    sources = point_data([-179.9, -20.0, 50.0], [0.0, 40.0, -20.0])
    queries = point_data([180.1, -18.0, 51.0], [0.0, 41.0, -19.0])
    idx, d = compute_dist_between_datasets(queries, sources)
    brute = np.array(
        [
            haversine((la, lo), (sources.lat.values, sources.lon.values))
            for la, lo in zip(queries.lat.values, queries.lon.values)
        ]
    )
    np.testing.assert_array_equal(idx, brute.argmin(axis=1))
    np.testing.assert_allclose(d, brute.min(axis=1), atol=1e-7)


def test_smoothing_measurement_error_and_scatter():
    ds = point_data([0.0], [0.0]).assign(trend=("x", [2.0]), trend_un=("x", [1.0]))
    out = radial_smoothing(ds)
    assert out.trend_s_un.item() == 1
    with pytest.raises(ValueError):
        radial_smoothing(ds, L_km=0)
    ds = point_data([0.0, 1.5], [0.0, 0.0]).assign(
        trend=("x", [0.0, 10.0]), trend_un=("x", [1.0, 1.0])
    )
    assert radial_smoothing(ds, L_km=100, kernel="gaussian").trend_s[0] > 0


def test_valid_neighbor_count():
    ds = point_data([0.0, 0.01], [0.0, 0.0]).assign(
        trend=("x", [2.0, np.nan]), trend_un=("x", [1.0, 1.0])
    )
    assert radial_smoothing(ds, min_neighbors=2).trend_s.isnull().all()


def test_merge_size_and_longitudes():
    rec = point_data([0.0], [0.0]).assign(
        VLM_trend_coefficient_mean=("x", [0.0]),
        VLM_trend_coefficient_uncertainty=("x", [0.001]),
        GRW_ensemble_std=(("x", "time"), np.zeros((1, 2))),
    )
    rec = rec.assign_coords(time=[2020, 2150])
    gps = point_data([-155.0, -154.0, -153.0], [20.0, 20.0, 20.0]).assign(
        trend=("x", [1.0, 2.0, 3.0]),
        trend_un=("x", [1.0, 1.0, 1.0]),
        trend_CMR_corrected=("x", [1.0, 2.0, 3.0]),
    )
    a, _ = merge_gps_and_vlm_reconstruction(rec, gps)
    b, _ = merge_gps_and_vlm_reconstruction(rec, gps.assign_coords(lon=gps.lon % 360))
    xr.testing.assert_allclose(a, b)
    assert (a.kind == 1).sum() == 3


@pytest.mark.parametrize(
    "years", [[2000, 2040], [2200, 2300], [2100], [2000, 2060, 2100]]
)
def test_decay_calendar_and_empty_window(years):
    y = np.array(years)
    weights = temp_decay(y)[2]
    assert np.all((weights >= 0) & (weights <= 1))
    if 2100 in years:
        assert weights[list(years).index(2100)] == pytest.approx(0.5)


def test_integral_persists_and_rebases():
    np.testing.assert_allclose(
        integrated_decay([2050, 2100, 2150, 2300], 2000), [50, 90.915494309, 100, 100]
    )
    y = np.array([2000, 2100, 2150, 2300])
    np.testing.assert_allclose(
        integrated_decay(y, 2100),
        integrated_decay(y, 2000) - integrated_decay([2100], 2000),
    )


def test_projection_years_exact():
    np.testing.assert_array_equal(
        projection_years(2005, 2000, 2100, 30), [2000, 2030, 2060, 2090]
    )
    np.testing.assert_array_equal(projection_years(2000, 2100, 2100, 1), [2100])
    for args in [(1999, 2000, 2100, 10), (2000, 2100, 2000, 10), (2000, 2000, 2100, 0)]:
        with pytest.raises(ValueError):
            projection_years(*args)


def prepared(v=-0.002, sig=0.001, gia=0.0, gia_un=0.5, rad=0.0, grw=None, w=1.0):
    def field(value):
        return xr.DataArray(
            [[value]], dims=("lat", "lon"), coords={"lat": [0.0], "lon": [0.0]}
        )

    if grw is None:
        grw = np.zeros((1, 4))
    return {
        "VLM_REC_MERGED": xr.Dataset(
            {
                "VLM_trend_coefficient_mean": field(v),
                "VLM_trend_coefficient_uncertainty": field(sig),
                "source_index": field(0),
            }
        ),
        "weights": field(w),
        "GIA_VLM_int": xr.Dataset(
            {
                "asl": field(gia + rad),
                "rsl_mean": field(gia),
                "rsl_sterr": field(gia_un),
            }
        ),
        "VLM_REC_MERGED_GRW": xr.DataArray(
            grw, dims=("x", "time"), coords={"time": [2000, 2020, 2100, 2150]}
        ),
    }


def test_rsl_sign_integrated_and_gia_cancellation():
    mean, *_ = projection_parameters(
        prepared(sig=0), [0], [0], [2000, 2050, 2150, 2300], 2000
    )
    np.testing.assert_allclose(mean[:, 0], [0, 100, 200, 200])
    mean, *_ = projection_parameters(
        prepared(v=0.003, gia=1, rad=3), [0], [0], [2100, 2150], 2000
    )
    np.testing.assert_allclose(mean[:, 0], [100, 150])


@pytest.mark.parametrize("mode", ["converge-to-gia", "constant"])
def test_rebase_full_trajectory_and_zero(mode):
    pp = prepared(grw=np.array([[0.0, 0.0, 0.01, 0.015]]))
    y = [2000, 2050, 2100, 2150]
    a = projection_parameters(pp, [0], [0], y, 2000, rate_mode=mode)
    b = projection_parameters(pp, [0], [0], y, 2100, rate_mode=mode)
    for x, z in zip(a, b):
        np.testing.assert_allclose(z, x - x[2], atol=1e-10)


def test_grw_dates_and_missing():
    grw = prepared(grw=np.array([[np.nan, 0, 0.01, 0.015]]))["VLM_REC_MERGED_GRW"]
    np.testing.assert_allclose(
        _grw_at_years(grw, np.array([0]), [2000, 2060, 2300]), [[0, 5, 0]]
    )
    bad = grw.copy()
    bad.values[0, 2] = np.nan
    with pytest.raises(ValueError, match="missing"):
        _grw_at_years(bad, np.array([0]), [2100])
    with pytest.raises(ValueError, match="coverage"):
        _grw_at_years(grw, np.array([0]), [1990])


def test_independent_draws_and_variance():
    draws = sample_draws(10000, 1234)
    assert abs(np.corrcoef(draws[0], draws[2])[0, 1]) < 0.03
    assert not np.array_equal(draws[0], draws[2])
    a = 2 * draws[0] + 3 * draws[2]
    assert a.var() == pytest.approx(13, rel=0.04)


@pytest.mark.parametrize("mode", ["converge-to-gia", "constant"])
def test_gia_fallback_without_observation_series(mode):
    pp = prepared(w=0.0, gia=1.0)
    pp["VLM_REC_MERGED"]["source_index"][:] = -1
    pp["VLM_REC_MERGED_GRW"] = pp["VLM_REC_MERGED_GRW"].isel(x=slice(0, 0))
    mean, obs, grw, gia = projection_parameters(
        pp, [0], [0], [2000, 2100], 2000, rate_mode=mode
    )
    np.testing.assert_allclose(mean[:, 0], [0, 100])
    assert not np.any(obs) and not np.any(grw)
    np.testing.assert_allclose(gia[:, 0], [0, 50])


@pytest.mark.parametrize("mode", ["converge-to-gia", "constant"])
def test_gia_only_matches_fallback_and_rebases_without_observations(mode):
    pp = prepared(w=0.0, gia=2.0, gia_un=0.4)
    pp["VLM_REC_MERGED"]["source_index"][:] = -1
    years = [2000, 2100, 2150, 2300]
    expected = projection_parameters(pp, [0], [0], years, 2100, rate_mode=mode)
    # No observational fields, weights, GRW data, or radial correction are needed.
    gia = pp["GIA_VLM_int"][["rsl_mean", "rsl_sterr"]]
    actual = projection_parameters(
        {"GIA_VLM_int": gia},
        [0],
        [0],
        years,
        2100,
        rate_mode=mode,
        gia_only=True,
    )
    for a, b in zip(actual, expected):
        np.testing.assert_array_equal(a, b)
    np.testing.assert_allclose(actual[0][:, 0], [-200, 0, 100, 400])
    np.testing.assert_allclose(actual[3][:, 0], [-40, 0, 20, 80])


def test_request_range_independence():
    pp = prepared()
    full = projection_parameters(pp, [0], [0], [2000, 2050, 2100, 2150, 2300], 2000)
    subset = projection_parameters(pp, [0], [0], [2100, 2300], 2000)
    for a, b in zip(full, subset):
        np.testing.assert_allclose(a[[2, 4]], b)


def test_grw_cf_calendar_and_unused_missing_source():
    import cftime

    grw = xr.DataArray(
        [[np.nan, np.nan], [0.0, 0.01]],
        dims=("x", "time"),
        coords={
            "time": [
                cftime.DatetimeNoLeap(2000, 1, 1),
                cftime.DatetimeNoLeap(2100, 1, 1),
            ]
        },
    )
    np.testing.assert_allclose(
        _grw_at_years(grw, np.array([-1, 1]), [2050]), [[0.0], [5.0]]
    )


def test_periodic_grid_dim_order_and_no_mutation():
    da = xr.DataArray(
        np.tile([1.0, 0.0, -1.0, 0.0, 1.0], (3, 1)).T,
        dims=("lon", "lat"),
        coords={"lon": [0, 90, 180, 270, 360], "lat": [-90, 0, 90]},
    )
    original = da.copy(deep=True)
    interp = periodic_interp(da, [0], [-180, 180, -1, 359])
    np.testing.assert_allclose(interp.values[0], interp.values[1])
    np.testing.assert_allclose(interp.values[2], interp.values[3])
    xr.testing.assert_identical(da, original)


def test_flatten_time():
    da = xr.DataArray(
        np.arange(12).reshape(2, 2, 3),
        dims=("time", "lat", "lon"),
        coords={"time": [2000, 2010], "lat": [0, 1], "lon": [0, 1, 2]},
    )
    out = flatten_lon_lat(da)
    np.testing.assert_array_equal(out.sla.values, da.values.reshape(2, 6))


def test_location_parser(tmp_path):
    p = tmp_path / "sites.lst"
    p.write_text("\n  # comment\nTest\t1\t0\t359\n")
    assert ReadLocationFile(p)[3][0] == 359
    for text in ["", "T\t1\t100\t0\n", "T\t1\t0\tnan\n", "T\t1\t0\t0\nX\t1\t0\t0\n"]:
        p.write_text(text)
        with pytest.raises(ValueError):
            ReadLocationFile(p)


@pytest.mark.parametrize("mode", ["converge-to-gia", "constant"])
def test_netcdf_streaming_precision_determinism_schema(tmp_path, mode):
    loc = tmp_path / "site.lst"
    loc.write_text("Test\t1\t0\t0\n")
    paths = [tmp_path / "a.nc", tmp_path / "b.nc"]
    for chunk, path in zip([1, 17], paths):
        vlm_postprocess(
            prepared(v=-0.5),
            51,
            1234,
            loc,
            2000,
            2000,
            2300,
            10,
            chunk,
            "test",
            str(path),
            rate_mode=mode,
        )
    with xr.open_dataset(paths[0]) as a, xr.open_dataset(paths[1]) as b:
        xr.testing.assert_equal(a, b)
        assert a.sea_level_change.dims == ("samples", "years", "locations")
        assert a.sea_level_change.dtype == np.float32
        assert a.sea_level_change.attrs["units"] == "mm"
        assert a.baseyear == 2000
        assert a.vlm_rate_mode == mode
        assert a.sea_level_change.sel(years=2000).max() == 0
        assert a.sea_level_change.sel(years=2300).mean() > 40000
        assert np.any(a.sea_level_change.values % 1 != 0)
        assert np.isnan(a.sea_level_change.encoding["_FillValue"])


def test_rate_mode_default_and_constant_secular_rates():
    pp = prepared(
        v=-0.002,
        sig=0.001,
        gia=1.0,
        gia_un=0.5,
        rad=3.0,
        grw=np.array([[0.0, 0.0, 0.01, 0.015]]),
    )
    years = [2000, 2050, 2100, 2150, 2200, 2300]
    default = projection_parameters(pp, [0], [0], years, 2000)
    taper = projection_parameters(
        pp, [0], [0], years, 2000, rate_mode="converge-to-gia"
    )
    constant = projection_parameters(pp, [0], [0], years, 2000, rate_mode="constant")
    for a, b in zip(default, taper):
        np.testing.assert_array_equal(a, b)
    # Observational RSL rate = -(-2) + (1+3) = 6 mm/year.
    np.testing.assert_allclose(constant[0][:, 0], (np.array(years) - 2000) * 6)
    assert taper[0][-1, 0] - taper[0][-2, 0] == pytest.approx(100.0)  # GIA 1 mm/year
    assert constant[0][-1, 0] - constant[0][-2, 0] == pytest.approx(600.0)
    assert taper[1][-1, 0] == pytest.approx(100.0)
    assert constant[1][-1, 0] == pytest.approx(300.0)
    assert constant[3][-1, 0] == 0  # no fallback-GIA error where s=1
    np.testing.assert_array_equal(constant[2], taper[2])  # GRW model unchanged


def test_constant_rate_retains_spatial_blend():
    pp = prepared(v=-0.002, gia=1.0, rad=3.0, w=0.3)
    mean, obs, grw, gia = projection_parameters(
        pp, [0], [0], [2000, 2300], 2000, rate_mode="constant"
    )
    assert mean[-1, 0] == pytest.approx(300 * (0.3 * 6 + 0.7 * 1))
    assert obs[-1, 0] == pytest.approx(90.0)
    assert gia[-1, 0] == pytest.approx(105.0)


def test_invalid_rate_mode():
    with pytest.raises(ValueError, match="rate_mode"):
        projection_parameters({}, [0], [0], [2100], 2000, rate_mode="unknown")


@pytest.mark.parametrize(
    "nsamps,chunk,seed", [(0, 1, 0), (1, 0, 0), (1, -1, 0), (1, 1, -1)]
)
def test_invalid_sampling(tmp_path, nsamps, chunk, seed):
    with pytest.raises(ValueError):
        vlm_postprocess(
            prepared(),
            nsamps,
            seed,
            "unused",
            2000,
            2000,
            2100,
            10,
            chunk,
            "test",
            str(tmp_path / "out.nc"),
        )
