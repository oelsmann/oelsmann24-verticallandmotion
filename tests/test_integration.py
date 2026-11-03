from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner
import numpy as np
import pytest
import xarray as xr

from oelsmann24_verticallandmotion.cli import main
from oelsmann24_verticallandmotion.vlm_preprocess import (
    load_datasets,
    preprocess_gia_only,
)
from oelsmann24_verticallandmotion.vlm_postprocess import sample_draws


@pytest.fixture
def inputs(tmp_path):
    rec = xr.Dataset(
        {
            "VLM_trend_coefficient_mean": ("x", [-0.002]),
            "VLM_trend_coefficient_uncertainty": ("x", [0.001]),
            "GRW_ensemble_std": (
                ("time", "x"),
                np.array([[np.nan], [0.0], [0.01], [0.015]]),
            ),
        },
        coords={
            "x": [5],
            "lon": ("x", [0.0]),
            "lat": ("x", [0.0]),
            "time": np.array(
                ["2000-12-31", "2020-12-31", "2100-12-31", "2150-12-31"],
                dtype="datetime64[ns]",
            ),
        },
    )
    gps = xr.Dataset(
        {k: ("x", [1.0, 2.0]) for k in ["trend", "trend_un", "trend_CMR_corrected"]},
        coords={"x": [0, 1], "lon": ("x", [205.0, 206.0]), "lat": ("x", [20.0, 20.0])},
    )
    gridcoords = {"lon": [0.0, 90.0, 180.0, 270.0, 360.0], "lat": [-90.0, 0.0, 90.0]}
    gia = xr.Dataset(
        {
            k: (("lon", "lat"), np.full((5, 3), v))
            for k, v in [("rsl_mean", 1.0), ("rsl_sterr", 0.2), ("rad_mean", 0.0)]
        },
        coords=gridcoords,
    )
    coast = xr.Dataset(
        {"z": (("y", "x"), np.zeros((3, 5)))},
        coords={"x": gridcoords["lon"], "y": gridcoords["lat"]},
    )
    datasets = [rec, gps, gia, coast, xr.Dataset()]
    paths = []
    for i, ds in enumerate(datasets):
        p = tmp_path / f"in{i}.nc"
        ds.to_netcdf(p)
        paths.append(str(p))
    return paths


def arguments(inputs, tmp_path):
    loc = tmp_path / "locations.lst"
    loc.write_text("Test\t1\t0\t0\n")
    keys = [
        "input-vlm-reconstruction-file",
        "input-gps-vlm-file",
        "input-gia-vlm-file",
        "input-dist2coast-file",
        "input-example-file",
    ]
    args = [word for key, value in zip(keys, inputs) for word in ["--" + key, value]]
    return args + [
        "--location-file",
        str(loc),
        "--nsamps",
        "21",
        "--baseyear",
        "2000",
        "--pyear-start",
        "2000",
        "--pyear-end",
        "2150",
        "--pyear-step",
        "10",
        "--output-lslr-file",
        str(tmp_path / "output.nc"),
    ]


@pytest.mark.parametrize(
    "mode,expected", [("converge-to-gia", 350.0), ("constant", 450.0)]
)
def test_full_cli_synthetic(inputs, tmp_path, mode, expected):
    result = CliRunner().invoke(
        main, arguments(inputs, tmp_path) + ["--vlm-rate-mode", mode]
    )
    assert result.exit_code == 0, repr(result.exception)
    with xr.open_dataset(tmp_path / "output.nc") as out:
        assert out.sea_level_change.shape == (21, 16, 1)
        assert out.sea_level_change.sel(years=2150).mean() == pytest.approx(
            expected, abs=1e-3
        )
        assert out.vlm_rate_mode == mode


def test_all_observations_excluded_falls_back_to_gia(inputs, tmp_path):
    with xr.open_dataset(inputs[0]) as ds:
        rec = ds.load().assign_coords(lon=("x", [-100.0]), lat=("x", [60.0]))
    rec.to_netcdf(inputs[0])
    with xr.open_dataset(inputs[1]) as ds:
        gps = ds.load().assign_coords(lon=("x", [0.0, 1.0]), lat=("x", [0.0, 0.0]))
    gps.to_netcdf(inputs[1])
    result = CliRunner().invoke(main, arguments(inputs, tmp_path))
    assert result.exit_code == 0, repr(result.exception)
    with xr.open_dataset(tmp_path / "output.nc") as out:
        assert out.sea_level_change.sel(years=2150).mean() == pytest.approx(
            150.0, abs=1e-3
        )


@pytest.mark.parametrize(
    "flag,value",
    [
        ("--nsamps", "0"),
        ("--chunk-size", "-2"),
        ("--seed", "-1"),
        ("--pyear-step", "0"),
        ("--pyear-end", "1990"),
        ("--vlm-rate-mode", "unknown"),
    ],
)
def test_cli_rejects_before_preprocess(inputs, tmp_path, flag, value):
    args = arguments(inputs, tmp_path) + [flag, value]
    with patch(
        "oelsmann24_verticallandmotion.cli.global_preprocess_verticallandmotion"
    ) as pp:
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 2
        pp.assert_not_called()


def test_cli_protects_inputs(inputs, tmp_path):
    args = arguments(inputs, tmp_path) + ["--output-lslr-file", inputs[0]]
    with patch(
        "oelsmann24_verticallandmotion.cli.global_preprocess_verticallandmotion"
    ) as pp:
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 2
        pp.assert_not_called()


def test_loader_rejects_wrong_units(inputs):
    path = Path(inputs[0])
    with xr.open_dataset(path) as source:
        ds = source.load()
    ds.VLM_trend_coefficient_mean.attrs["units"] = "mm/year"
    ds.to_netcdf(path)
    with pytest.raises(ValueError, match="expected m/year"):
        load_datasets(*inputs)


@pytest.mark.parametrize("mode", ["converge-to-gia", "constant"])
def test_gia_only_cli_analytic_samples_without_other_inputs(inputs, tmp_path, mode):
    args = arguments(inputs, tmp_path)
    for flag in [
        "--input-vlm-reconstruction-file",
        "--input-gps-vlm-file",
        "--input-dist2coast-file",
        "--input-example-file",
    ]:
        i = args.index(flag)
        Path(args[i + 1]).unlink()
        del args[i : i + 2]
    # A GIA-only file need not contain radial motion, already included in RSL.
    with xr.open_dataset(inputs[2]) as source:
        gia = source[["rsl_mean", "rsl_sterr"]].load()
    gia.rsl_mean[0, 1] = 2.0
    gia.rsl_sterr[0, 1] = 0.4
    gia.to_netcdf(inputs[2])
    (tmp_path / "locations.lst").write_text("Coast\t1\t0\t0\nOcean\t2\t0\t90\n")
    args += [
        "--gia-only",
        "--vlm-rate-mode",
        mode,
        "--baseyear",
        "2100",
        "--pyear-end",
        "2300",
        "--pyear-step",
        "100",
        "--input-example-file",
        str(tmp_path / "ignored-missing-example.nc"),
    ]
    outputs = []
    for chunk in [1, 17]:
        with patch(
            "oelsmann24_verticallandmotion.cli.global_preprocess_verticallandmotion",
            side_effect=AssertionError("GIA-only mode must bypass observations"),
        ):
            result = CliRunner().invoke(main, args + ["--chunk-size", str(chunk)])
        assert result.exit_code == 0, repr(result.exception)
        with xr.open_dataset(tmp_path / "output.nc") as out:
            outputs.append(out.load())
    xr.testing.assert_equal(*outputs)
    out = outputs[0]
    z = sample_draws(21, 1234)[2]
    dt = np.array([-100, 0, 100, 200])
    expected = (
        dt[None, :, None] * np.array([2.0, 1.0])[None, None, :]
        + z[:, None, None] * dt[None, :, None] * np.array([0.4, 0.2])[None, None, :]
    )
    np.testing.assert_allclose(out.sea_level_change.values, expected, rtol=1e-6)
    assert out.vlm_rate_mode == "gia-only"
    assert out.baseyear == 2100
    assert out.sea_level_change.dims == ("samples", "years", "locations")
    assert out.sea_level_change.dtype == np.float32
    assert out.sea_level_change.attrs["units"] == "mm"
    assert np.isnan(out.sea_level_change.encoding["_FillValue"])


@pytest.mark.parametrize(
    "flag",
    [
        "--input-vlm-reconstruction-file",
        "--input-gps-vlm-file",
        "--input-dist2coast-file",
        "--input-example-file",
    ],
)
def test_blended_cli_still_requires_observation_inputs(inputs, tmp_path, flag):
    args = arguments(inputs, tmp_path)
    i = args.index(flag)
    del args[i : i + 2]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 2
    assert f"Missing option '{flag}'" in result.output
    assert not (tmp_path / "output.nc").exists()


def test_gia_only_cli_protects_input(inputs, tmp_path):
    args = arguments(inputs, tmp_path) + ["--gia-only", "--output-lslr-file", inputs[2]]
    with patch("oelsmann24_verticallandmotion.cli.preprocess_gia_only") as pp:
        result = CliRunner().invoke(main, args)
    assert result.exit_code == 2
    assert "must not overwrite" in result.output
    pp.assert_not_called()


@pytest.mark.parametrize(
    "problem,match",
    [
        ("missing", "missing variables"),
        ("units", "must be in mm/year"),
        ("negative", "cannot be negative"),
        ("nonfinite", "contains missing values"),
    ],
)
def test_gia_only_validates_rates_and_errors(inputs, problem, match):
    with xr.open_dataset(inputs[2]) as source:
        gia = source.load()
    if problem == "missing":
        gia = gia.drop_vars("rsl_sterr")
    elif problem == "units":
        gia.rsl_mean.attrs["units"] = "m/year"
    elif problem == "negative":
        gia.rsl_sterr[:] = -1
    else:
        gia.rsl_mean[:] = np.nan
    gia.to_netcdf(inputs[2])
    with pytest.raises(ValueError, match=match):
        preprocess_gia_only(inputs[2])
