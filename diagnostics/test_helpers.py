"""Small optional diagnostics tests; no public input files are needed."""

import numpy as np
import pytest
import xarray as xr

pytest.importorskip("cartopy", reason="Install the optional diagnostics group")
from diagnostics import helpers as dh  # noqa: E402


def diagnostic_dataset():
    sites = dh.location_table()
    return xr.Dataset(
        {
            "mean": ("locations", sites.lat.to_numpy() + sites.lon.to_numpy()),
            "lat": ("locations", sites.lat.to_numpy()),
            "lon": ("locations", sites.lon.to_numpy()),
        },
        coords={"locations": sites.id.to_numpy()},
    )


def test_complete_global_sites_with_unchanged_cities():
    sites = dh.location_table()
    assert len(sites) == 64806
    assert sites.id.is_unique
    assert tuple(sites.name[:6]) == tuple(city[0] for city in dh.CITIES)
    np.testing.assert_array_equal(sites.lat[:6], [city[1] for city in dh.CITIES])
    np.testing.assert_array_equal(sites.lon[:6], [city[2] for city in dh.CITIES])
    nodes = sites.iloc[6:]
    assert nodes[["lat", "lon"]].drop_duplicates().shape[0] == 64800
    assert (nodes.lat.min() - 0.5, nodes.lat.max() + 0.5) == (-90, 90)
    assert (nodes.lon.min() - 0.5, nodes.lon.max() + 0.5) == (-180, 180)
    # This central Pacific cell must not be removed by a coastal/ocean mask.
    assert ((nodes.lat == 0.5) & (nodes.lon == -149.5)).any()


def test_map_cells_preserve_values_and_reject_incomplete_or_unordered_grid():
    ds = diagnostic_dataset()
    field = dh.global_map_field(ds, "mean")
    assert field.shape == (180, 360)
    np.testing.assert_array_equal(field, dh.MAP_LATS[:, None] + dh.MAP_LONS[None, :])
    for invalid in [
        ds.isel(locations=slice(None, -1)),
        ds.isel(locations=slice(None, None, -1)),
    ]:
        with pytest.raises(ValueError, match="complete, ordered global"):
            dh.global_map_field(invalid, "mean")


def test_map_rejects_nonfinite_cells_instead_of_plotting_silent_holes():
    ds = diagnostic_dataset()
    ds["mean"].values[-1] = np.nan
    with pytest.raises(ValueError, match="nonfinite"):
        dh.global_map_field(ds, "mean")
