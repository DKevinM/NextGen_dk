import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# --------------------------------------------------
# CONFIG
# --------------------------------------------------

# CRS that uses metres (important for distances); adjust to your usual one
TARGET_CRS = "EPSG:4326"  

# File paths (adjust to your repo)
ACA_SHP   = "NextGen_dk/data/ACA_Boundary_2022.shp"
WCAS_SHP  = "NextGen_dk/data/WCAS_2024.shp"

# Remote data sources (GitHub raw)
STATIONS_URL = "https://raw.githubusercontent.com/DKevinM/AB_datapull/main/data/last6h.csv"
PURPLE_URL   = "https://raw.githubusercontent.com/DKevinM/AB_datapull/main/data/AB_PM25_map.json"


# Grid resolutions (metres)
ACA_CELL_M   = 1000   # 2 km
WCAS_CELL_M  = 1000   # 2 km  (change if you want 1 km)
# Later: COUNTY_CELL_M = 500

# Weights
STATION_WEIGHT = 1.0
PURPLE_WEIGHT  = 0.5  # or 0.5

# Max influence distance (metres) – optional
MAX_DIST_M = 50_000  # 50 km


# --------------------------------------------------
# AQHI helpers
# --------------------------------------------------

def est_aqhi_from_pm(pm_val: float) -> float | None:
    """Estimated AQHI from PM2.5: round(PM/10)+1, capped at 10."""
    try:
        v = float(pm_val)
    except (TypeError, ValueError):
        return None
    aq = round(v / 10.0) + 1
    if aq < 0:
        return None
    return min(aq, 10.0)


# --------------------------------------------------
# Data loading
# --------------------------------------------------

def load_airshed(shp_path: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(shp_path)
    if gdf.crs is None:
        raise ValueError(f"{shp_path} has no CRS; please set one before using.")
    return gdf.to_crs(TARGET_CRS)


def load_stations(csv_path: str) -> gpd.GeoDataFrame:
    """
    Expect columns: lon, lat, AQHI (adjust names as needed).
    If you already have a shapefile, switch to gpd.read_file().
    """
    df = pd.read_csv(csv_path)
    # adjust these to your actual column names
    lon_col = "lon"
    lat_col = "lat"
    aqhi_col = "AQHI"

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df[lon_col], df[lat_col]),
        crs="EPSG:4326"
    ).to_crs(TARGET_CRS)

    # numeric AQHI, capped at 10
    gdf["aqhi_val"] = pd.to_numeric(gdf[aqhi_col], errors="coerce")
    gdf["aqhi_val"] = gdf["aqhi_val"].clip(upper=10)

    gdf["weight"] = STATION_WEIGHT
    return gdf


def load_purple(csv_path: str) -> gpd.GeoDataFrame:
    """
    Expect columns: lon, lat, and some pm25 field(s) (pm_corr, pm25, PM2_5, etc.).
    Adjust the field names to match your actual schema.
    """
    df = pd.read_csv(csv_path)

    lon_col_candidates = ["lon", "longitude", "LON"]
    lat_col_candidates = ["lat", "latitude", "LAT"]

    lon_col = next((c for c in lon_col_candidates if c in df.columns), None)
    lat_col = next((c for c in lat_col_candidates if c in df.columns), None)

    if lon_col is None or lat_col is None:
        raise ValueError("Could not find lon/lat columns for PurpleAir data.")

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df[lon_col], df[lat_col]),
        crs="EPSG:4326"
    ).to_crs(TARGET_CRS)

    # PM2.5 candidates; adjust to your actual columns
    pm_candidates = ["pm_corr", "pm25", "PM2_5", "PM2.5", "pm_25"]
    pm_col = next((c for c in pm_candidates if c in gdf.columns), None)

    if pm_col is None:
        raise ValueError("No PM2.5 column found for PurpleAir data.")

    gdf["pm_val"] = pd.to_numeric(gdf[pm_col], errors="coerce")
    gdf["aqhi_val"] = gdf["pm_val"].apply(est_aqhi_from_pm)
    gdf["weight"]   = PURPLE_WEIGHT

    return gdf


# --------------------------------------------------
# Grid creation inside polygon
# --------------------------------------------------

def make_grid_points_within_polygon(poly_gdf: gpd.GeoDataFrame, cellsize_m: float) -> gpd.GeoDataFrame:
    """
    Create a regular grid of points (cell centers) inside the polygon extent,
    then keep only those whose point lies inside the polygon.
    """
    poly_union = poly_gdf.unary_union
    minx, miny, maxx, maxy = poly_union.bounds

    xs = np.arange(minx, maxx + cellsize_m, cellsize_m)
    ys = np.arange(miny, maxy + cellsize_m, cellsize_m)

    points = []
    for x in xs:
        for y in ys:
            p = Point(x, y)
            if poly_union.contains(p):
                points.append(p)

    grid = gpd.GeoDataFrame(geometry=points, crs=poly_gdf.crs)
    grid["x"] = grid.geometry.x
    grid["y"] = grid.geometry.y
    return grid


# --------------------------------------------------
# Blended IDW (AQHI) with source-type weights
# --------------------------------------------------

def blended_idw_aqhi(grid_gdf: gpd.GeoDataFrame,
                     pts_gdf: gpd.GeoDataFrame,
                     power: float = 2.0,
                     max_dist_m: float | None = None) -> np.ndarray:
    """
    For each grid point, compute blended AQHI using:
      AQHI_j = Σ( (w_i * AQHI_i) / d_ij^p ) / Σ( w_i / d_ij^p )
    where w_i is source weight (station vs PurpleAir).
    """
    gx = grid_gdf["x"].to_numpy()
    gy = grid_gdf["y"].to_numpy()

    px = pts_gdf.geometry.x.to_numpy()
    py = pts_gdf.geometry.y.to_numpy()
    val = pts_gdf["aqhi_val"].to_numpy()
    wt  = pts_gdf["weight"].to_numpy()

    n_grid = gx.size
    out = np.full(n_grid, np.nan)

    for j in range(n_grid):
        dx = px - gx[j]
        dy = py - gy[j]
        d  = np.sqrt(dx*dx + dy*dy)

        # ignore invalid or zero distance
        valid = np.isfinite(d) & (d > 0)
        if max_dist_m is not None and np.isfinite(max_dist_m):
            valid &= (d <= max_dist_m)

        if not np.any(valid):
            continue

        d  = d[valid]
        v  = val[valid]
        w0 = wt[valid]

        # drop NaN AQHI
        finite = np.isfinite(v)
        if not np.any(finite):
            continue

        d  = d[finite]
        v  = v[finite]
        w0 = w0[finite]

        w_idw = w0 / (d**power)

        num = np.sum(w_idw * v)
        den = np.sum(w_idw)

        if den > 0:
            out[j] = num / den

    return out


# --------------------------------------------------
# Master function: build ACA & WCAS grids
# --------------------------------------------------

def build_airshed_grids(shp_path: str,
                        cellsize_m: float,
                        stations_gdf: gpd.GeoDataFrame,
                        purple_gdf: gpd.GeoDataFrame,
                        outfile_geojson: str | None = None) -> gpd.GeoDataFrame:
    """
    Build two AQHI IDW grids inside the given airshed:
      - aqhi_stn   : stations-only IDW (weight = 1)
      - aqhi_blend : stations + PurpleAir blended (station=1, PA=PURPLE_WEIGHT)
    """
    # 1) Load airshed polygon and ensure CRS
    poly = load_airshed(shp_path)

    # 2) Prepare stations-only points (aqhi_val + weight)
    stn = stations_gdf.copy()
    stn = stn[stn["aqhi_val"].notna()].copy()
    stn["weight"] = STATION_WEIGHT  # should already be this, but just to be explicit

    # 3) Prepare blended points (stations + PurpleAir)
    pa  = purple_gdf.copy()
    pa  = pa[pa["aqhi_val"].notna()].copy()

    pts_blend = pd.concat([stn, pa], ignore_index=True)
    pts_blend = gpd.GeoDataFrame(pts_blend, geometry="geometry", crs=TARGET_CRS)

    # 4) Make grid points inside the airshed
    grid = make_grid_points_within_polygon(poly, cellsize_m)

    # 5) Stations-only IDW
    aqhi_stn = blended_idw_aqhi(
        grid_gdf=grid,
        pts_gdf=stn,
        power=2.0,
        max_dist_m=MAX_DIST_M
    )

    # 6) Blended IDW (stations + PurpleAir)
    aqhi_blend = blended_idw_aqhi(
        grid_gdf=grid,
        pts_gdf=pts_blend,
        power=2.0,
        max_dist_m=MAX_DIST_M
    )

    # 7) Attach results and cap at 10
    grid["aqhi_stn"]   = np.clip(aqhi_stn, None, 10)
    grid["aqhi_blend"] = np.clip(aqhi_blend, None, 10)

    # Optional display categories
    for col in ["aqhi_stn", "aqhi_blend"]:
        cat_col = col + "_cat"
        vals = grid[col].to_numpy()
        cat = np.full(vals.shape, None, dtype=object)

        mask = np.isfinite(vals)
        cat[mask] = vals[mask].round().astype(int).astype(str)
        cat[(mask) & (vals >= 10)] = "10+"

        grid[cat_col] = cat

    if outfile_geojson is not None:
        grid.to_file(outfile_geojson, driver="GeoJSON")

    return grid


if __name__ == "__main__":
    # Load point data once
    stations_gdf = load_stations(STATIONS_CSV)
    purple_gdf   = load_purple(PURPLE_CSV)

    # ACA grid
    aca_grid = build_airshed_grids(
        shp_path=ACA_SHP,
        cellsize_m=ACA_CELL_M,
        stations_gdf=stations_gdf,
        purple_gdf=purple_gdf,
        outfile_geojson="data/ACA_AQHI_station_vs_blended.geojson"
    )

    # WCAS grid
    wcas_grid = build_airshed_grids(
        shp_path=WCAS_SHP,
        cellsize_m=WCAS_CELL_M,
        stations_gdf=stations_gdf,
        purple_gdf=purple_gdf,
        outfile_geojson="data/WCAS_AQHI_station_vs_blended.geojson"
    )
