"""Reproducible diagnostics around the production FACTS2 package.

Scientific transformations stay in src/oelsmann24_verticallandmotion. This module
only selects diagnostic sites, runs that package, summarizes samples, and plots.
"""

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from oelsmann24_verticallandmotion.vlm_postprocess import vlm_postprocess
from oelsmann24_verticallandmotion.vlm_preprocess import (
    global_preprocess_verticallandmotion,
)

ROOT = Path(__file__).resolve().parents[1]
FILES = {
    "input_vlm_reconstruction_file": "VLM_reconstruction.nc",
    "input_gps_vlm_file": "NGL14_CMR_corrected.nc",
    "input_gia_vlm_file": "GIA_stats.nc",
    "input_dist2coast_file": "dist2coast_1deg_v2.grd",
    "input_example_file": "ssp585.2150.fair2.emuGLA.emulandice2.glaciers_quantiles.nc",
}
CITIES = (
    ("Venice", 45.43, 12.33),
    ("Skellefteå", 64.75, 21.00),
    ("London", 51.5074, -0.1278),
    ("Tokyo", 35.65, 139.77),
    ("Brisbane", -27.47, 153.03),
    ("New York", 40.70, -74.01),
)
MODES = {"converge-to-gia": "Converge to GIA", "constant": "Constant secular rates"}
COLORS = {"converge-to-gia": "#176782", "constant": "#a54f22"}
MAP_LATS = np.arange(-89.5, 90, 1.0)
MAP_LONS = np.arange(-179.5, 180, 1.0)


@dataclass(frozen=True)
class Config:
    input_dir: Path = ROOT / "data/inputs/reference"
    output_dir: Path = ROOT / "data/diagnostics/shared_notebooks/global_grid"
    baseyear: int = 2000
    start: int = 2000
    end: int = 2300
    step: int = 10
    nsamps: int = 500
    seed: int = 1234

    def paths(self):
        missing = [
            name for name in FILES.values() if not (self.input_dir / name).is_file()
        ]
        if missing:
            raise FileNotFoundError(
                "Missing public inputs: "
                + ", ".join(missing)
                + ". Obtain https://doi.org/10.5281/zenodo.18199757 manually and "
                "set INPUT_DIR in the notebook. No data are downloaded automatically."
            )
        return {key: str(self.input_dir / name) for key, name in FILES.items()}


def preprocess(config):
    return global_preprocess_verticallandmotion(
        "diagnostic_notebooks", **config.paths()
    )


def inventory(config):
    """Inspect every source variable without loading the multi-GB coastal grid."""
    config.paths()
    rows = []
    for filename in FILES.values():
        with xr.open_dataset(config.input_dir / filename) as ds:
            for name, variable in ds.data_vars.items():
                rows.append(
                    {
                        "file": filename,
                        "variable": name,
                        "dimensions": str(dict(variable.sizes)),
                        "dtype": str(variable.dtype),
                        "units": variable.attrs.get(
                            "units", variable.attrs.get("unit", "not declared")
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def provenance(config):
    config.paths()
    settings = {
        k: v for k, v in asdict(config).items() if k not in {"input_dir", "output_dir"}
    }
    return {
        "settings": settings,
        "cities": CITIES,
        "map_sites": {
            "selection": "all global one-degree cell centers; no coastal/ocean/land mask",
            "latitude_centers": [-89.5, 89.5],
            "longitude_centers": [-179.5, 179.5],
            "spacing_degrees": 1,
        },
        "input_sha256": {
            name: _sha256(config.input_dir / name) for name in FILES.values()
        },
        "source_sha256": {
            p.name: _sha256(p)
            for p in sorted((ROOT / "src/oelsmann24_verticallandmotion").glob("*.py"))
        },
        "helper_sha256": _sha256(Path(__file__).resolve()),
        "versions": {
            name: version(name)
            for name in ["numpy", "scipy", "xarray", "netCDF4", "matplotlib", "cartopy"]
        },
    }


def location_table():
    """Six cities plus every global 1° cell center, including the open ocean."""
    lon, lat = np.meshgrid(MAP_LONS, MAP_LATS)
    table = pd.DataFrame(CITIES, columns=["name", "lat", "lon"])
    table["id"] = np.arange(1, len(CITIES) + 1)
    maps = pd.DataFrame({"lat": lat.ravel(), "lon": lon.ravel()})
    maps["id"] = np.arange(1000, 1000 + len(maps))
    maps["name"] = [f"global_{n}" for n in maps.id]
    return pd.concat([table, maps], ignore_index=True)[["name", "id", "lat", "lon"]]


def global_map_field(dataset, name):
    """Reshape actual projected global nodes into complete map cells, not dots.

    City rows are excluded. Coordinates are checked so an old coastal-only cache
    cannot silently produce a map with missing ocean cells.
    """
    nodes = dataset.isel(locations=slice(len(CITIES), None))
    lat, lon = np.meshgrid(MAP_LATS, MAP_LONS, indexing="ij")
    if (
        nodes.sizes["locations"] != lat.size
        or not np.array_equal(nodes.lat.values, lat.ravel())
        or not np.array_equal(nodes.lon.values, lon.ravel())
    ):
        raise ValueError("Expected the complete, ordered global one-degree map grid")
    values = nodes[name].transpose("locations").values
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Global map field {name} contains nonfinite values")
    return xr.DataArray(
        values.reshape(lat.shape),
        dims=("lat", "lon"),
        coords={"lat": MAP_LATS, "lon": MAP_LONS},
        attrs=nodes[name].attrs,
    )


def summarize_samples(raw_path, summary_path):
    """Actual finite-ensemble means/std in bounded location blocks (ddof=0)."""
    with xr.open_dataset(raw_path) as ds:
        means = np.empty((ds.sizes["years"], ds.sizes["locations"]), dtype=np.float32)
        stds = np.empty_like(means)
        for start in range(0, ds.sizes["locations"], 512):
            section = slice(start, start + 512)
            block = ds.sea_level_change.isel(locations=section).values.astype(
                np.float64
            )
            means[:, section] = block.mean(axis=0)
            stds[:, section] = block.std(axis=0, ddof=0)
        result = xr.Dataset(
            {
                "mean": (("years", "locations"), means),
                "std": (("years", "locations"), stds),
                "lat": ds.lat,
                "lon": ds.lon,
            },
            coords={"years": ds.years, "locations": ds.locations},
            attrs=dict(ds.attrs),
        )
        result.attrs["std_definition"] = (
            "Population standard deviation of generated samples, ddof=0"
        )
        result.attrs["projection_scope"] = (
            "VLM/GIA contribution to RSL, not total climate-driven RSL"
        )
        for name in ["mean", "std"]:
            result[name].attrs["units"] = "mm"
        result.to_netcdf(
            summary_path,
            encoding={name: {"dtype": "f4", "zlib": True} for name in ["mean", "std"]},
        )


def ensure_projections(config, force=False):
    """Recompute when inputs, production code, helper, versions or settings change."""
    expected = provenance(config)
    fingerprint = hashlib.sha256(
        json.dumps(expected, sort_keys=True).encode()
    ).hexdigest()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = config.output_dir / "manifest.json"
    summary_paths = {mode: config.output_dir / f"{mode}_mean_std.nc" for mode in MODES}
    if (
        not force
        and manifest_path.is_file()
        and all(p.is_file() for p in summary_paths.values())
    ):
        previous = json.loads(manifest_path.read_text())
        if previous.get("fingerprint") == fingerprint and all(
            previous.get("summary_sha256", {}).get(mode) == _sha256(path)
            for mode, path in summary_paths.items()
        ):
            return {
                mode: xr.load_dataset(path) for mode, path in summary_paths.items()
            }, previous

    pp = preprocess(config)
    table = location_table()
    location_path = config.output_dir / "locations.lst"
    table.to_csv(location_path, sep="\t", header=False, index=False)
    for mode, path in summary_paths.items():
        raw = config.output_dir / f"{mode}_samples.nc"
        vlm_postprocess(
            pp,
            config.nsamps,
            config.seed,
            str(location_path),
            config.baseyear,
            config.start,
            config.end,
            config.step,
            50,
            "diagnostic_notebooks",
            str(raw),
            rate_mode=mode,
        )
        summarize_samples(raw, path)
    manifest = {
        "fingerprint": fingerprint,
        "provenance": expected,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "git_note": "Source hashes identify the actual code including uncommitted corrections.",
        "locations": len(table),
        "global_map_nodes": len(table) - len(CITIES),
        "summary_sha256": {mode: _sha256(path) for mode, path in summary_paths.items()},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return {
        mode: xr.load_dataset(path) for mode, path in summary_paths.items()
    }, manifest


def map_axes(nrows=1, ncols=1, figsize=None):
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize or (8 * ncols, 4.6 * nrows),
        subplot_kw={"projection": ccrs.Robinson()},
        squeeze=False,
        layout="constrained",
    )
    for ax in axes.flat:
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor="#eeeeeb", zorder=0)
        ax.coastlines(resolution="110m", linewidth=0.4, color="#666666")
    return fig, axes


def map_points(ax, lon, lat, values, title, units, cmap="RdBu_r", limits=None):
    values = np.asarray(values)
    if limits is None:
        upper = max(float(np.nanpercentile(abs(values), 98)), 0.01)
        limits = (-upper, upper) if cmap == "RdBu_r" else (0, upper)
    im = ax.scatter(
        (np.asarray(lon) + 180) % 360 - 180,
        lat,
        c=values,
        s=3,
        cmap=cmap,
        vmin=limits[0],
        vmax=limits[1],
        transform=ccrs.PlateCarree(),
        rasterized=True,
    )
    ax.set_title(title)
    ax.get_figure().colorbar(
        im,
        ax=ax,
        orientation="horizontal",
        pad=0.025,
        shrink=0.85,
        label=units,
        extend="both" if limits[0] < 0 else "max",
    )


def map_grid(ax, field, title, units, cmap="RdBu_r", limits=None):
    field = field.transpose("lat", "lon").sortby("lat").sortby("lon")
    if limits is None:
        upper = max(float(np.nanpercentile(abs(field.values), 98)), 0.01)
        limits = (-upper, upper) if cmap == "RdBu_r" else (0, upper)
    im = ax.pcolormesh(
        field.lon,
        field.lat,
        field,
        cmap=cmap,
        vmin=limits[0],
        vmax=limits[1],
        transform=ccrs.PlateCarree(),
        rasterized=True,
    )
    ax.set_title(title)
    ax.get_figure().colorbar(
        im,
        ax=ax,
        orientation="horizontal",
        pad=0.025,
        shrink=0.85,
        label=units,
        extend="both" if limits[0] < 0 else "max",
    )


def mark_cities(ax, labels=False):
    for index, (name, lat, lon) in enumerate(CITIES, 1):
        ax.plot(
            lon, lat, "o", ms=3, mec="black", mfc="white", transform=ccrs.PlateCarree()
        )
        if labels:
            ax.text(
                lon + 4, lat + 2, str(index), fontsize=8, transform=ccrs.PlateCarree()
            )


def city_plot(results, baseyear):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, layout="constrained")
    for index, (ax, (name, _, _)) in enumerate(zip(axes.flat, CITIES)):
        for mode, ds in results.items():
            mean, std = (
                ds["mean"].isel(locations=index),
                ds["std"].isel(locations=index),
            )
            ax.plot(ds.years, mean, color=COLORS[mode], lw=2, label=MODES[mode])
            ax.fill_between(
                ds.years, mean - std, mean + std, color=COLORS[mode], alpha=0.14
            )
        ax.axvspan(2050, 2150, color="#666666", alpha=0.045)
        ax.axhline(0, color=".65", lw=0.7)
        ax.grid(alpha=0.2)
        ax.set(
            title=name,
            xlabel="Year",
            ylabel=f"RSL contribution (mm; relative to {baseyear})",
        )
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Corrected model: mean ± one ensemble standard deviation", fontsize=15)
    return fig


def projection_maps(results, year):
    selected = {mode: ds.sel(years=year) for mode, ds in results.items()}
    grids = {
        mode: {name: global_map_field(ds, name) for name in ["mean", "std"]}
        for mode, ds in selected.items()
    }
    mean_limit = max(
        float(
            np.percentile(
                np.concatenate([abs(g["mean"].values).ravel() for g in grids.values()]),
                98,
            )
        ),
        1.0,
    )
    std_limit = max(
        float(
            np.percentile(
                np.concatenate([g["std"].values.ravel() for g in grids.values()]), 98
            )
        ),
        1.0,
    )
    fig, axes = map_axes(2, 2, figsize=(15, 9))
    for row, (mode, grid) in enumerate(grids.items()):
        map_grid(
            axes[row, 0],
            grid["mean"],
            f"{MODES[mode]} — {year} mean",
            "mm",
            limits=(-mean_limit, mean_limit),
        )
        map_grid(
            axes[row, 1],
            grid["std"],
            f"{MODES[mode]} — {year} standard deviation",
            "mm (1σ)",
            cmap="magma",
            limits=(0, std_limit),
        )
        for ax in axes[row]:
            mark_cities(ax)
    return fig


def difference_maps(results, year):
    fig, axes = map_axes(1, 2, figsize=(15, 5))
    for ax, name in zip(axes.flat, ["mean", "std"]):
        delta = global_map_field(
            results["constant"].sel(years=year), name
        ) - global_map_field(results["converge-to-gia"].sel(years=year), name)
        map_grid(
            ax,
            delta,
            f"Constant − convergence: {year} {name}",
            "mm",
        )
    return fig


def save_figure(fig, config, name):
    directory = config.output_dir / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / name, dpi=140, bbox_inches="tight")
    plt.show()


def main():
    """Run model summaries without Jupyter: uv run --group diagnostics python diagnostics/helpers.py."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Config.input_dir)
    parser.add_argument("--output-dir", type=Path, default=Config.output_dir)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    _, manifest = ensure_projections(
        Config(input_dir=args.input_dir, output_dir=args.output_dir), force=args.force
    )
    print(
        json.dumps(
            {
                key: manifest[key]
                for key in ["fingerprint", "locations", "global_map_nodes"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
