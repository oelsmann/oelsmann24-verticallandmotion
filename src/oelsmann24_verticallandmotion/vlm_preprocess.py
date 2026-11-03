import numpy as np
import xarray as xr
from scipy.spatial import cKDTree
from typing import Optional, Tuple
import copy
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

""" vlm_preprocess.py

This runs the post-processing stage for the global vertical land motion component.

Parameters:
nsamps = Number of samples to draw
rng_seed = Seed value for the random number generator
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

"""


EARTH_RADIUS_KM = 6371.0088


# ── Helper functions: spherical geometry & kernels ─────────────────────────────
def _to_xyz(lon_deg, lat_deg):
    """Convert lon/lat (deg) to unit-sphere 3D coordinates for fast KDTree queries.

    Parameters
    ----------
    lon_deg, lat_deg : array-like
        Longitudes and latitudes in degrees.

    Returns
    -------
    np.ndarray of shape (N, 3)
        Unit vectors (x, y, z) on the sphere.
    """
    lon_deg, lat_deg = np.broadcast_arrays(lon_deg, lat_deg)
    if not np.all(np.isfinite(lon_deg) & np.isfinite(lat_deg)):
        raise ValueError("Coordinates must be finite")
    if np.any(np.abs(lat_deg) > 90):
        raise ValueError("Latitudes must be between -90 and 90 degrees")
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    cl = np.cos(lat)
    x = cl * np.cos(lon)
    y = cl * np.sin(lon)
    z = np.sin(lat)
    return np.column_stack((x, y, z))


def _arc_km_from_chord(chord):
    """Convert chord length on the unit sphere to great-circle arc length (km)."""
    arc_rad = 2.0 * np.arcsin(np.clip(chord * 0.5, 0.0, 1.0))
    return EARTH_RADIUS_KM * arc_rad


def _kernel_weights(dist_km, L_km, kernel="cos", exponent=1.0):
    """
    Compute radial distance weights on [0, L_km]; outside the radius -> 0.

    Parameters
    ----------
    dist_km : array-like
        Distances in km from the evaluation point.
    L_km : float
        Kernel radius/scale in km.
    kernel : {'cos', 'tri', 'gaussian'}
        - 'cos'     : raised cosine (Hann) on [0, L]
        - 'tri'     : linear taper to 0 at L
        - 'gaussian': exp(-(d/L)^2), truncated at 3*L for speed
    exponent : float
        Optional exponent applied to weights (w ** exponent).

    Returns
    -------
    np.ndarray
        Weights aligned with `dist_km`.
    """
    if not np.isfinite(L_km) or L_km <= 0 or not np.isfinite(exponent) or exponent <= 0:
        raise ValueError("Kernel scale and exponent must be finite and positive")
    d = np.asarray(dist_km)
    w = np.zeros_like(d, dtype=float)

    if kernel == "cos":  # raised cosine (Hann) on [0,L]
        inside = d <= L_km
        w[inside] = 0.5 * (1.0 + np.cos(np.pi * d[inside] / L_km))
    elif kernel == "tri":  # linear taper
        inside = d <= L_km
        w[inside] = 1.0 - d[inside] / L_km
    elif kernel == "gaussian":
        w = np.exp(-((d / L_km) ** 2))
        inside = d <= (3 * L_km)  # cut tail for speed
        w[~inside] = 0.0
    else:
        raise ValueError("kernel must be one of {'cos','tri','gaussian'}")

    if exponent != 1.0:
        w = w ** float(exponent)
    return w


def radial_smoothing(
    ds: xr.Dataset,
    *,
    var_name: str = "trend",
    uncert_name: str = "trend_un",  # if missing, uncertainty weighting is skipped
    dim: str = "x",
    lon_name: str = "lon",
    lat_name: str = "lat",
    L_km: float = 150.0,  # radius/scale (km)
    kernel: str = "cos",
    exponent: float = 1.0,
    use_uncert: bool = True,
    min_neighbors: int = 1,
) -> xr.Dataset:
    """
    Radially smooth `ds[var_name]` using distance-decay weights and optional 1/σ² weights.

    Output variables added to a copy of `ds`:
      - `<var_name>_s`    : uncertainty- and distance-weighted mean
      - `<var_name>_s_un` : weighted standard deviation among neighbors

    Notes
    -----
    - Neighbor search uses a KDTree on unit-sphere 3D vectors with a chord-length cutoff
      equivalent to the specified arc radius L_km.
    - If `uncert_name` is absent or `use_uncert=False`, only distance weights are used.
    """

    if var_name not in ds:
        raise KeyError(f"{var_name!r} not found in dataset")

    # Extract arrays
    lon = np.asarray(ds[lon_name].values)
    lat = np.asarray(ds[lat_name].values)
    y = np.asarray(ds[var_name].values).astype(float)

    # Optional uncertainty weights (1/variance)
    if use_uncert and (uncert_name in ds):
        sigma = np.asarray(ds[uncert_name].values).astype(float)
        sigma = np.where((sigma > 0) & np.isfinite(sigma), sigma, np.nan)
        w_unc = 1.0 / (sigma**2)
        w_unc = np.where(np.isfinite(w_unc), w_unc, 0.0)
        used_uncert = True
    else:
        w_unc = None
        used_uncert = False

    # KDTree on the sphere
    xyz = _to_xyz(lon, lat)
    tree = cKDTree(xyz)

    # Radius in chord space for neighbor search
    if min_neighbors < 1:
        raise ValueError("min_neighbors must be positive")
    _kernel_weights(np.array([0.0]), L_km, kernel, exponent)
    arc = min((3 if kernel == "gaussian" else 1) * L_km / EARTH_RADIUS_KM, np.pi)
    chord_radius = 2.0 * np.sin(arc / 2.0)

    y_s = np.full_like(y, np.nan, dtype=float)
    y_s_un = np.full_like(y, np.nan, dtype=float)

    neighbors_list = tree.query_ball_point(xyz, r=chord_radius)

    for i, idxs in enumerate(neighbors_list):
        if len(idxs) < min_neighbors:
            continue

        # Distances (km) from point i to its neighbors
        euclid = np.linalg.norm(xyz[idxs] - xyz[i], axis=1)
        d_km = _arc_km_from_chord(euclid)

        # Distance-decay weights
        w = _kernel_weights(d_km, L_km, kernel=kernel, exponent=exponent)

        # Combine with 1/sigma^2 if provided
        if w_unc is not None:
            w = w * w_unc[idxs]

        # Exclude NaNs and zero weights
        valid = np.isfinite(y[idxs]) & (w > 0)
        if np.count_nonzero(valid) < min_neighbors:
            continue

        wv = w[valid]
        yv = y[idxs][valid]

        # Weighted mean
        wsum = np.sum(wv)
        ybar = np.sum(wv * yv) / wsum
        y_s[i] = ybar

        # Weighted standard deviation of neighbors
        # sqrt( sum(w*(x-mu)^2) / sum(w) )
        variance = np.sum(wv * (yv - ybar) ** 2) / wsum
        if w_unc is not None:
            # Independent measurement error of the weighted mean, plus unresolved
            # spatial scatter. An isolated station retains its reported error.
            variance += np.sum((wv / wsum) ** 2 * sigma[idxs][valid] ** 2)
        y_s_un[i] = np.sqrt(variance)

    # Package outputs
    out = ds.copy()
    base_coords = {dim: ds[dim], lon_name: ds[lon_name], lat_name: ds[lat_name]}

    out[f"{var_name}_s"] = xr.DataArray(
        y_s,
        dims=(dim,),
        coords=base_coords,
        attrs={
            "description": f"Radially smoothed {var_name} (weighted mean)",
            "L_km": L_km,
            "kernel": kernel,
            "exponent": exponent,
            "uncert_weighting": used_uncert,
        },
    )
    out[f"{var_name}_s_un"] = xr.DataArray(
        y_s_un,
        dims=(dim,),
        coords=base_coords,
        attrs={
            "description": f"Propagated measurement error plus spatial scatter for {var_name}_s",
            "L_km": L_km,
            "kernel": kernel,
            "exponent": exponent,
            "uncert_weighting": used_uncert,
            "formula": "sqrt(sum(a*(x-mu)^2) + sum(a^2*sigma^2)); a=w/sum(w)",
        },
    )
    return out


def couple_short(
    a: xr.DataArray,
    b: xr.DataArray,
    *,
    dim: str = "x",
    lon_name: str = "lon",
    lat_name: str = "lat",
    id_coord_a: Optional[str] = None,  # ID for a (if coord present; else row index)
    id_coord_b: Optional[str] = None,  # ID for b
    limit_km: float = 100.0,
    return_distance: bool = True,
) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray, bool]:
    """
    Match the *shorter* dataset to the *other* by nearest neighbor on a sphere.

    Returns arrays of length = len(shorter set):
      - ids_short : IDs of the shorter dataset (for reference)
      - ids_other : IDs from the other dataset (NaN if nearest > limit_km)
      - dist_km   : great-circle distance (km) for each row (if return_distance=True)
      - short_is_a: True if `a` was the shorter set; False if `b` was shorter
    """

    # Extract lon/lat as numpy arrays
    lon_a = np.asarray(a[lon_name].values)
    lat_a = np.asarray(a[lat_name].values)
    lon_b = np.asarray(b[lon_name].values)
    lat_b = np.asarray(b[lat_name].values)

    # Helper to pick ID coord or fallback to row index
    def get_ids(da: xr.DataArray, id_coord: Optional[str]) -> np.ndarray:
        if id_coord is not None and id_coord in da.coords:
            return np.asarray(da[id_coord].values)
        return np.arange(da.sizes[dim], dtype=float)  # float so we can set NaN

    ids_a = get_ids(a, id_coord_a)
    ids_b = get_ids(b, id_coord_b)

    # Decide which is shorter → this will be the query set
    if a.sizes[dim] <= b.sizes[dim]:
        lon_short, lat_short, ids_short_np = lon_a, lat_a, ids_a
        lon_other, lat_other, ids_other_all = lon_b, lat_b, ids_b
        short_is_a = True
    else:
        lon_short, lat_short, ids_short_np = lon_b, lat_b, ids_b
        lon_other, lat_other, ids_other_all = lon_a, lat_a, ids_a
        short_is_a = False

    n = len(ids_short_np)

    # KDTree on the *other* set in 3D (unit sphere); query with the *short* set
    tree = cKDTree(_to_xyz(lon_other, lat_other))
    euclid, nn_idx = tree.query(_to_xyz(lon_short, lat_short), k=1)

    # Euclidean chord length -> great-circle arc length -> distance in km
    arc_rad = 2.0 * np.arcsin(np.clip(euclid * 0.5, 0.0, 1.0))
    dist_km = EARTH_RADIUS_KM * arc_rad

    # IDs of the nearest "other" points for each short point; NaN if beyond limit
    ids_other_for_short = ids_other_all[nn_idx].astype(float)
    ids_other_for_short[dist_km > limit_km] = np.nan

    # Build xarray outputs aligned to the shorter set; give lon/lat coords explicit dims
    base_coords = {
        dim: np.arange(n),
        lon_name: (dim, lon_short),
        lat_name: (dim, lat_short),
    }

    ids_short_da = xr.DataArray(
        ids_short_np.astype(float), dims=(dim,), coords=base_coords, name="id_short"
    )
    ids_other_da = xr.DataArray(
        ids_other_for_short, dims=(dim,), coords=base_coords, name="id_other"
    )

    if return_distance:
        dist_da = xr.DataArray(
            dist_km, dims=(dim,), coords=base_coords, name="distance_km"
        )
        return ids_short_da, ids_other_da, dist_da, short_is_a

    return ids_short_da, ids_other_da, None, short_is_a


def make_grw_grid(VLM_REC_MERGED_GRW_large_GRW, dis_sel_2d, ids_2):
    """
    Build a time-by-location dataset for GRW ensemble std on the target (lat, lon, time) grid.

    Parameters
    ----------
    VLM_REC_MERGED_GRW_large_GRW : xr.Dataset
        Dataset containing 'GRW_ensemble_std' with a time dimension.
    dis_sel_2d : xr.Dataset
        Target grid providing 'lat' and 'lon' coordinates (2D).
    ids_2 : np.ndarray
        Indices mapping from flattened target grid to source vector.

    Returns
    -------
    xr.Dataset
        Dataset with variable 'GRW_ensemble_std(lat, lon, time)' on the target grid.
    """
    all_datasets = []

    lat_a = dis_sel_2d["lat"].values
    lon_a = dis_sel_2d["lon"].values

    for time_id in np.arange(
        len(VLM_REC_MERGED_GRW_large_GRW["GRW_ensemble_std"].time)
    ):
        data_set = (
            np.empty([len(dis_sel_2d["VLM_trend_coefficient_mean"].values.flatten())])
            * np.nan
        )
        np.put(
            data_set,
            list(ids_2.astype(int)),
            VLM_REC_MERGED_GRW_large_GRW["GRW_ensemble_std"][:, time_id],
        )
        all_datasets.append(data_set)
    all_datasets = np.asarray(all_datasets).T
    dis_sel_2d_grw = xr.Dataset(
        data_vars=dict(
            GRW_ensemble_std=(
                ["lat", "lon", "time"],
                all_datasets.reshape(
                    len(lat_a),
                    len(lon_a),
                    len(VLM_REC_MERGED_GRW_large_GRW["GRW_ensemble_std"].time),
                ),
            ),
        ),
        coords=dict(
            lat=("lat", lat_a),
            lon=("lon", lon_a),
        ),
    )
    return dis_sel_2d_grw


def regrid_from_flat(data_in, variables, lon_a, lat_a):
    """
    Re-grid 1D scattered values back to a structured (lat, lon) grid.

    Parameters
    ----------
    data_in : xr.Dataset
        Must contain coords 'lon' and 'lat' and variables listed in `variables`.
    variables : list[str]
        Variables in `data_in` to re-grid.
    lon_a, lat_a : array-like
        Target 1D lon/lat coordinates.

    Returns
    -------
    (xr.Dataset, np.ndarray)
        - Dataset with each variable on shape (lat, lon).
        - Index map from source points to target flattened grid.
    """
    xx, yy = np.meshgrid(lat_a, lon_a)
    xx = xx.T.flatten()
    yy = yy.T.flatten()

    data_set = np.zeros([len(yy)]) * 0
    ids_ = np.empty([len(data_in.lon.values)]) * 0
    i = 0

    for lond, latd in zip(data_in.lon.values, data_in.lat.values):
        target = np.argwhere((xx == latd) & (yy == lond))
        ids_[i] = int(target[0, 0])
        i = i + 1

    variable = variables[0]
    data_set = np.zeros([len(yy)]) * np.nan
    np.put(data_set, list(ids_.astype(int)), data_in[variable])

    ds = xr.Dataset(
        data_vars=dict(
            data=(["lat", "lon"], data_set.reshape(len(lat_a), len(lon_a))),
        ),
        coords=dict(
            lat=("lat", lat_a),
            lon=("lon", lon_a),
        ),
    )
    ds = ds.rename({"data": variable})
    if len(variables) > 1:
        for var in variables[1:]:
            data_set = np.zeros([len(yy)]) * np.nan
            np.put(data_set, list(ids_.astype(int)), data_in[var].values)
            ds[var] = (("lat", "lon"), data_set.reshape(len(lat_a), len(lon_a)))

    return ds, ids_


def haversine(coord1, coord2):
    """
    Great-circle distance(s) using the haversine formula.

    Parameters
    ----------
    coord1 : tuple(float, float)
        (lat1, lon1) in degrees.
    coord2 : tuple(float, float) or arrays
        (lat2, lon2) in degrees. Can be arrays to get vectorized distances.

    Returns
    -------
    float or np.ndarray
        Distance(s) in meters.
    """
    R = EARTH_RADIUS_KM * 1000  # Earth radius in meters
    lat1, lon1 = coord1
    lat2, lon2 = coord2

    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)

    a = np.clip(
        np.sin(dphi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2) ** 2,
        0,
        1,
    )
    return 2 * R * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def compute_dist_between_datasets(data1, data2):
    """
    For each (lon, lat) in `data1`, find nearest neighbor in `data2` and return
    index and distance.

    Returns
    -------
    (np.ndarray, np.ndarray)
        - ids_  : indices into `data2` for nearest neighbors
        - dist_ : distances (meters) to those neighbors
    """
    if data2.lon.size == 0:
        raise ValueError("Cannot find neighbors in an empty source dataset")
    tree = cKDTree(_to_xyz(data2.lon.values, data2.lat.values))
    chord, ids_ = tree.query(_to_xyz(data1.lon.values, data1.lat.values))
    return ids_, _arc_km_from_chord(chord) * 1000


def flatten_lon_lat(datain, name="sla", time=True):
    """
    Flatten a 2D (lat, lon) field into a 1D coordinate 'x' with preserved coords.

    Parameters
    ----------
    datain : xr.DataArray
    name : str
        Output variable name.
    time : bool
        If True, keep 'time' dimension in the output.

    Returns
    -------
    xr.Dataset
    """
    lon, lat = np.meshgrid(datain.lon, datain.lat)
    flat_x = datain.transpose(
        *(["time", "lat", "lon"] if time else ["lat", "lon"])
    ).values
    coords = {"lon": (["x"], lon.flatten()), "lat": (["x"], lat.flatten())}
    if time:
        coords["time"] = datain.time.values
        var_dict = {name: (["time", "x"], flat_x.reshape(datain.sizes["time"], -1))}
    else:
        var_dict = {name: (["x"], flat_x.ravel())}

    ds = xr.Dataset(var_dict, coords=coords)
    return ds


def merge_gps_and_vlm_reconstruction(VLM_REC, GPS_VLM):
    """
    Merge GPS trends with VLM reconstruction using uncertainty-aware smoothing.

    Returns
    -------
    (xr.Dataset, xr.Dataset)
        - VLM_REC_MERGED     : concatenated dataset combining GPS trends and the VLM reconstruction
        - VLM_REC_MERGED_GRW : corresponding GRW ensemble std on merged set

    Notes
    -----
    - 'kind' = 0 for VLM reconstruction, 1 for GPS-smoothed contribution.
    - GRW_ensemble_std represents present-day variability (Gaussian random walks).
    - all units are in m/year or m
    """

    VLM_REC = VLM_REC.assign_coords(lon=(VLM_REC.lon + 180) % 360 - 180).copy()
    GPS_VLM = GPS_VLM.assign_coords(lon=(GPS_VLM.lon + 180) % 360 - 180).copy()
    VLM_REC = VLM_REC.drop_vars("x", errors="ignore")
    GPS_VLM = GPS_VLM.drop_vars("x", errors="ignore")
    _, distances = compute_dist_between_datasets(GPS_VLM, VLM_REC)

    # select GPS at locations not covered by the reconstruction
    GPS_VLM_sub = GPS_VLM.isel({"x": distances > 100000})[
        ["trend", "trend_un", "trend_CMR_corrected"]
    ]
    GPS_select = [
        [
            -35.68,
            10.44,
            -20.37,
            42.82,
        ],  # Canary Islands region (off NW Africa in the Atlantic)
        [-179.1, 7.0, -138.6, 40.2],  # Hawaii and central Pacific Ocean
        [
            51.9,
            -56.6,
            92.4,
            -2.2,
        ],  # Western/central Indian Ocean (Madagascar, Mauritius, Seychelles)
        [
            27.0,
            -56.6,
            92.4,
            -35.7,
        ],  # Southwestern Indian Ocean (Mozambique Channel, around Madagascar)
        [
            115.49,
            1.01,
            131.67,
            21.04,
        ],  # Maritime Southeast Asia (Borneo, Philippines, Sulawesi)
        [-177.4, -53.1, -101.6, 6.7],  # South Pacific Ocean (Fiji, Tonga, Samoa region)
    ]

    select_gps_sub = np.asarray([False] * len(GPS_VLM_sub.lon))
    for box in GPS_select:
        select_gps_sub = select_gps_sub | (
            (GPS_VLM_sub.lon < box[2])
            & (GPS_VLM_sub.lat < box[3])
            & (GPS_VLM_sub.lon > box[0])
            & (GPS_VLM_sub.lat > box[1])
        )

    GPS_VLM_sub = GPS_VLM_sub.sel({"x": select_gps_sub})
    GPS_VLM_sub["trend"][:] = GPS_VLM_sub["trend_CMR_corrected"][:].values

    # compute spatial averages based on GPS trends
    smoothed = radial_smoothing(
        GPS_VLM_sub,
        var_name="trend",
        uncert_name="trend_un",  # skip or set use_uncert=False if not available
        dim="x",
        lon_name="lon",
        lat_name="lat",
        L_km=100.0,
        kernel="cos",
        exponent=1.0,
        use_uncert=True,
    )

    # match resolution with the 0.25 degree resolution of the VLM reconstruction
    smoothed["lon"][:] = np.round(smoothed["lon"] * 4) / 4
    smoothed["lat"][:] = np.round(smoothed["lat"] * 4) / 4
    smoothed = (
        smoothed.drop_vars("x", errors="ignore").rename(
            {
                "trend_s": "VLM_trend_coefficient_mean",
                "trend_s_un": "VLM_trend_coefficient_uncertainty",
            }
        )
        / 1000.0
    )
    smoothed["kind"] = copy.deepcopy(smoothed["VLM_trend_coefficient_uncertainty"])
    smoothed["kind"][:] = [1] * len(smoothed["kind"])

    VLM_REC["kind"] = copy.deepcopy(VLM_REC["VLM_trend_coefficient_uncertainty"])
    VLM_REC["kind"][:] = [0] * len(VLM_REC["kind"])
    VLM_REC_MERGED = xr.concat(
        [
            VLM_REC[
                [
                    "VLM_trend_coefficient_mean",
                    "VLM_trend_coefficient_uncertainty",
                    "kind",
                ]
            ],
            smoothed,
        ],
        dim="x",
    )
    smoothed_grw = (
        smoothed[["VLM_trend_coefficient_uncertainty"]]
        .expand_dims(dim={"time": VLM_REC["GRW_ensemble_std"].time}, axis=0)[
            "VLM_trend_coefficient_uncertainty"
        ]
        .rename("GRW_ensemble_std")
        .T.to_dataset()
        * 0
    )

    # GRW_ensemble_std is the ensemble std based on the Gaussian random walks, used to mimic present day variability
    VLM_REC_MERGED_GRW = xr.concat(
        [VLM_REC[["GRW_ensemble_std"]], smoothed_grw], dim="x"
    )

    # lat_merged = VLM_REC_MERGED.lat.values
    lon_merged = VLM_REC_MERGED.lon.values
    select_merged_sub = np.asarray([False] * len(lon_merged))

    # select regions with very little GNSS station coverage where GIA trends will be applied
    unselect_data_boxes = [
        [-119, -39, 51, 91],  # Covers parts of North America
        [-180, -38, 68, 91],  # Covers eastern Siberia to Greenland/Arctic Ocean
        [52, 77, 10, 28],  # Covers Middle East / Arabian Peninsula
        [-20, 50, 69.1, 80],
    ]  # Covers Northern Europe and parts of the Arctic
    for box in unselect_data_boxes:
        select_merged_sub = select_merged_sub | (
            (VLM_REC_MERGED.lon > box[0])
            & (VLM_REC_MERGED.lat < box[3])
            & (VLM_REC_MERGED.lon < box[1])
            & (VLM_REC_MERGED.lat > box[2])
        )
    VLM_REC_MERGED = VLM_REC_MERGED.sel({"x": ~select_merged_sub.values})
    VLM_REC_MERGED_GRW = VLM_REC_MERGED_GRW.sel({"x": ~select_merged_sub.values})
    return VLM_REC_MERGED, VLM_REC_MERGED_GRW


def normalize_grid(ds):
    """Sort a rectilinear grid and normalize its longitude seam without mutation."""
    if "lon" not in ds.coords or "lat" not in ds.coords:
        raise ValueError("Grid requires lon and lat coordinates")
    _to_xyz(ds.lon.values, np.zeros(ds.lon.size))
    if not np.all(np.isfinite(ds.lat)) or np.any(abs(ds.lat) > 90):
        raise ValueError("Invalid latitude coordinates")
    ds = ds.assign_coords(lon=(ds.lon + 180) % 360 - 180).sortby("lon").sortby("lat")
    _, idx = np.unique(ds.lon.values, return_index=True)
    ds = ds.isel(lon=np.sort(idx))
    if np.unique(ds.lat).size != ds.lat.size:
        raise ValueError("Duplicate latitude coordinates")
    return ds


def periodic_interp(ds, lat, lon):
    """Interpolate across the longitude seam; latitude extrapolation only at poles."""
    ds = normalize_grid(ds)
    lon = (np.asarray(lon) + 180) % 360 - 180
    if ds.lon.size < 2 or ds.lat.size < 2:
        raise ValueError("Input grids require at least two latitudes and longitudes")
    # Extend only global grids. Regional inputs must cover their requested domain.
    spacing = float(np.median(np.diff(ds.lon.values)))
    if (
        float(ds.lon.max() - ds.lon.min()) >= 180
        and 360 - float(ds.lon.max() - ds.lon.min()) <= 2 * spacing
    ):
        ds = xr.concat(
            [
                ds.isel(lon=[-1]).assign_coords(lon=ds.lon.values[-1:] - 360),
                ds,
                ds.isel(lon=[0]).assign_coords(lon=ds.lon.values[:1] + 360),
            ],
            dim="lon",
        )
    if np.min(lon) < float(ds.lon.min()) or np.max(lon) > float(ds.lon.max()):
        raise ValueError("Input grid does not cover requested longitudes")
    lat_step = float(np.median(np.diff(ds.lat.values)))
    if (
        np.min(lat) < float(ds.lat.min()) - lat_step
        or np.max(lat) > float(ds.lat.max()) + lat_step
    ):
        raise ValueError("Input grid does not cover requested latitudes")
    return ds.interp(lat=lat, lon=lon, kwargs={"fill_value": "extrapolate"})


def _validate_gia_rates(gia, names):
    """Validate rate units and uncertainty for both blended and GIA-only runs."""
    for name in names:
        declared = gia[name].attrs.get("units", gia[name].attrs.get("unit", "mm/yr"))
        if declared.lower().replace(" ", "") not in {"mm/yr", "mm/year", "mmyr-1"}:
            raise ValueError(f"GIA {name} must be in mm/year")
    if np.any(gia.rsl_sterr < 0):
        raise ValueError("GIA standard errors cannot be negative")


def _interpolate_gia(gia, include_radial=True):
    """Use the same quarter-degree GIA grid for both projection paths."""
    gia = periodic_interp(
        gia, np.arange(-90, 90, 0.25), np.arange(-180, 180, 0.25)
    ).transpose("lat", "lon")
    names = ["rsl_mean", "rsl_sterr"] + (["rad_mean"] if include_radial else [])
    for name in names:
        if not np.all(np.isfinite(gia[name])):
            raise ValueError(f"GIA {name} contains missing values after interpolation")
    if include_radial:
        gia["asl"] = gia.rsl_mean + gia.rad_mean
    return gia


def preprocess_gia_only(input_gia_vlm_file):
    """Load only GIA RSL means/errors, bypassing all observational input data."""
    names = ["rsl_mean", "rsl_sterr"]
    with xr.open_dataset(input_gia_vlm_file) as source:
        absent = set(names) - set(source)
        if absent:
            raise ValueError(
                f"{input_gia_vlm_file}: missing variables {sorted(absent)}"
            )
        gia = source[names].load()
    for name in names:
        if set(gia[name].dims) != {"lat", "lon"}:
            raise ValueError(f"GIA {name} must have lat and lon dimensions")
    _validate_gia_rates(gia, names)
    return {"GIA_VLM_int": _interpolate_gia(gia, include_radial=False)}


def interpolate_and_compute_weights(
    VLM_REC_MERGED, VLM_REC_MERGED_GRW, GIA_VLM, dist2coast
):
    """Prepare spatial fields, retaining GRW only at source sites with actual dates.

    A source_index grid replaces the former global space-by-time GRW cube. Site
    extraction happens before temporal interpolation in postprocessing.
    """
    lon_grid = np.arange(-180, 180, 0.25)
    lat_grid = np.arange(-90, 90, 0.25)
    # Keep the established lightweight coastal-grid decimation, but normalize the
    # seam before interpolation. Source native coordinate resolution may vary.
    coast = dist2coast["z"].transpose("lat", "lon")
    coast = periodic_interp(coast, lat_grid, lon_grid).transpose("lat", "lon")
    selected = np.isfinite(coast.values) & (abs(coast.values) < 150)
    rows, cols = np.nonzero(selected)
    sites = xr.Dataset(
        coords={"lon": ("x", lon_grid[cols]), "lat": ("x", lat_grid[rows])}
    )
    if VLM_REC_MERGED.sizes["x"]:
        ids, distances = compute_dist_between_datasets(sites, VLM_REC_MERGED)
    else:
        # All observations may have been excluded by the established masks.
        rows, cols = np.array([], dtype=int), np.array([], dtype=int)
        ids, distances = np.array([], dtype=int), np.array([])
    coords = {"lat": lat_grid, "lon": lon_grid}
    shape = selected.shape
    spatial = xr.Dataset(coords=coords)
    for name in ["VLM_trend_coefficient_mean", "VLM_trend_coefficient_uncertainty"]:
        values = np.full(shape, np.nan)
        values[rows, cols] = VLM_REC_MERGED[name].values[ids]
        spatial[name] = (("lat", "lon"), values)
    d = np.full(shape, np.nan)
    d[rows, cols] = distances
    spatial["distance"] = (("lat", "lon"), d)
    spatial["z"] = coast
    source_index = np.full(shape, -1, dtype=np.int32)
    source_index[rows, cols] = ids
    spatial["source_index"] = (("lat", "lon"), source_index)
    weights = xr.where(
        np.isfinite(spatial.distance) & (spatial.distance < 100000),
        (np.cos(np.pi * spatial.distance / 100000) + 1) / 2,
        0,
    )
    gia = _interpolate_gia(GIA_VLM)
    grw = VLM_REC_MERGED_GRW["GRW_ensemble_std"].transpose("x", "time")
    return gia, spatial, grw, weights


def load_datasets(
    input_vlm_reconstruction_file: str,
    input_gps_vlm_file: str,
    input_gia_vlm_file: str,
    input_dist2coast_file: str,
    input_example_file: str,
):
    """
    Load required datasets

    Returns
    -------
    (VLM_REC, GPS_VLM, GIA_VLM, dist2coast, example) : tuple of xr.Dataset
    """

    paths = [
        input_vlm_reconstruction_file,
        input_gps_vlm_file,
        input_gia_vlm_file,
        input_dist2coast_file,
        input_example_file,
    ]
    datasets = []
    # Load only necessary variables, closing every input handle before computation.
    required = [
        [
            "VLM_trend_coefficient_mean",
            "VLM_trend_coefficient_uncertainty",
            "GRW_ensemble_std",
        ],
        ["trend", "trend_un", "trend_CMR_corrected"],
        ["rsl_mean", "rsl_sterr", "rad_mean"],
        ["z"],
        [],
    ]
    for path, names in zip(paths, required):
        with xr.open_dataset(path) as source:
            absent = set(names) - set(source)
            if absent:
                raise ValueError(f"{path}: missing variables {sorted(absent)}")
            if "z" in names:
                source = source.rename(
                    {
                        k: v
                        for k, v in {"x": "lon", "y": "lat"}.items()
                        if k in source.dims
                    }
                )
                # Do not read the multi-GB native coastal raster into memory.
                lat_stride = max(
                    1, round(0.1 / float(np.median(abs(np.diff(source.lat)))))
                )
                lon_stride = max(
                    1, round(0.11 / float(np.median(abs(np.diff(source.lon)))))
                )
                source = source.isel(
                    lat=slice(None, None, lat_stride), lon=slice(None, None, lon_stride)
                )
            datasets.append(
                source[names].load() if names else xr.Dataset(attrs=dict(source.attrs))
            )
    VLM_REC, GPS_VLM, GIA_VLM, dist2coast, example = datasets
    for ds, mean, sigma, unit in [
        (
            VLM_REC,
            "VLM_trend_coefficient_mean",
            "VLM_trend_coefficient_uncertainty",
            "m/year",
        ),
        (GPS_VLM, "trend_CMR_corrected", "trend_un", "mm/year"),
    ]:
        if ds[mean].dims != ("x",) or ds[sigma].dims != ("x",):
            raise ValueError(f"{mean} and {sigma} must be one-dimensional on x")
        _to_xyz(ds.lon.values, ds.lat.values)
        if np.any(ds[sigma].values < 0):
            raise ValueError(f"{sigma} cannot be negative")
        valid = np.isfinite(ds[mean]) & np.isfinite(ds[sigma])
        if not bool(valid.all()):
            logger.warning(
                "Dropping %d sites missing %s or its uncertainty",
                int((~valid).sum()),
                mean,
            )
        filtered = ds.isel(x=valid)
        if ds is VLM_REC:
            VLM_REC = filtered
        else:
            GPS_VLM = filtered
        for var in [mean, sigma]:
            declared = ds[var].attrs.get("units", unit).lower().replace(" ", "")
            accepted = (
                {"m/year", "m/yr", "myr-1"}
                if unit == "m/year"
                else {"mm/year", "mm/yr", "mmyr-1"}
            )
            if declared not in accepted:
                raise ValueError(f"{var}: expected {unit}, got {declared}")
    if VLM_REC.sizes["x"] == 0:
        raise ValueError("Reconstruction contains no valid sites")
    if VLM_REC.GRW_ensemble_std.attrs.get("units", "m") != "m":
        raise ValueError("GRW_ensemble_std must contain displacement in metres")
    _validate_gia_rates(GIA_VLM, ["rsl_mean", "rsl_sterr", "rad_mean"])
    return VLM_REC, GPS_VLM, GIA_VLM, dist2coast, example


def global_preprocess_verticallandmotion(
    pipeline_id: str,
    input_vlm_reconstruction_file: str,
    input_gps_vlm_file: str,
    input_gia_vlm_file: str,
    input_dist2coast_file: str,
    input_example_file: str,
):
    """Load, interpolate and combine datasets:

    VLM_REC: VLM reconstruction (Oelsmann et al., 2024)

    GPS: VLM trends from NGL

    GIA_VLM: GIA VLM and effects on RSL

    dist2coast: distance to coast grid

    example: File with the same output format as IPPC AR6 SL projections

    Parameters
    ----------
    pipeline_id : str
        Unique identifier for the pipeline run.

    """

    # Load required datasets
    VLM_REC, GPS_VLM, GIA_VLM, dist2coast, example = load_datasets(
        input_vlm_reconstruction_file=input_vlm_reconstruction_file,
        input_gps_vlm_file=input_gps_vlm_file,
        input_gia_vlm_file=input_gia_vlm_file,
        input_dist2coast_file=input_dist2coast_file,
        input_example_file=input_example_file,
    )
    # Merge GPS and the VLM reconstruction
    VLM_REC_MERGED, VLM_REC_MERGED_GRW = merge_gps_and_vlm_reconstruction(
        VLM_REC, GPS_VLM
    )
    # Interpolate
    GIA_VLM_int, VLM_REC_MERGED, VLM_REC_MERGED_GRW, weights = (
        interpolate_and_compute_weights(
            VLM_REC_MERGED, VLM_REC_MERGED_GRW, GIA_VLM, dist2coast
        )
    )

    # Populate the output dictionary
    outdata = {
        "GIA_VLM_int": GIA_VLM_int,
        "VLM_REC_MERGED": VLM_REC_MERGED,
        "VLM_REC_MERGED_GRW": VLM_REC_MERGED_GRW,
        "weights": weights,
        "example": example,
    }

    return outdata
