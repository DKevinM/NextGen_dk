import json
import urllib.request
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path

# --------------------------------------------------
# CONFIG
# --------------------------------------------------

# CRS that uses metres (important for distances); adjust to your usual one
TARGET_CRS = "EPSG:3400"  

# File paths
BASE_DIR = Path(__file__).resolve().parent
ACA_SHP  = BASE_DIR / "data" / "ACA_Boundary_2022.shp"
WCAS_SHP = BASE_DIR / "data" / "WCAS_2024.shp"

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




def load_airshed(shp_path):
    from pathlib import Path
    p = Path(shp_path)
    print(f"[grid] Loading airshed shapefile from: {p.resolve()}")
    gdf = gpd.read_file(p)
    if gdf.crs is None:
        raise ValueError(f"{p} has no CRS; please set one in the shapefile.")
    return gdf.to_crs(TARGET_CRS)



# --------------------------------------------------
# Data loading
# --------------------------------------------------

def load_stations_from_url() -> gpd.GeoDataFrame:
    """
    Load station locations + *true* AQHI from last6h.csv.

    Logic:
      - Use StationName, Latitude/Longitude from the CSV.
      - Use AQHI directly:
          * if there is a wide 'AQHI'/'aqhi' column, use it, OR
          * if long format, keep only rows where ParameterName is an AQHI variant
            and take 'Value' as AQHI.
      - Take the *latest* AQHI per station based on datetime column.
      - Clip AQHI at 10 for safety, keep NaNs if not numeric.
    """
    df = pd.read_csv(STATIONS_URL)

    # ---- 1) Station name ----
    if "StationName" not in df.columns:
        raise ValueError("No StationName column in last6h.csv")
    stn_col = "StationName"

    # ---- 2) Latitude / longitude ----
    lat_candidates = ["Latitude", "lat", "Lat"]
    lon_candidates = ["Longitude", "lon", "Lon", "lng"]

    lat_col = next((c for c in lat_candidates if c in df.columns), None)
    lon_col = next((c for c in lon_candidates if c in df.columns), None)
    if lat_col is None or lon_col is None:
        raise ValueError("Could not find latitude/longitude columns in last6h.csv")

    # ---- 3) Find AQHI values ----
    # Case A: wide table with AQHI column on each row
    if "AQHI" in df.columns:
        aqhi_df = df.copy()
        aqhi_col = "AQHI"
    elif "aqhi" in df.columns:
        aqhi_df = df.copy()
        aqhi_col = "aqhi"
    # Case B: long table with ParameterName/Value
    elif {"ParameterName", "Value"}.issubset(df.columns):
        # include a few likely variants of the AQHI name;
        # tweak these strings if needed once you inspect the CSV.
        aqhi_names = [
            "AQHI",
            "Air Quality Health Index",
            "Air Quality Health Index (AQHI)"
        ]
        mask = df["ParameterName"].isin(aqhi_names)
        aqhi_df = df[mask].copy()
        aqhi_col = "Value"
        if aqhi_df.empty:
            raise ValueError(
                "Found ParameterName/Value columns but no rows where "
                "ParameterName looks like AQHI; check the exact label."
            )
    else:
        raise ValueError(
            "Could not find AQHI in last6h.csv "
            "(no AQHI column and no ParameterName/Value layout)."
        )

    # ---- 4) Datetime column to pick *latest* AQHI per station ----
    dt_candidates = ["ReadingDate", "DateTime", "date_time", "Timestamp"]
    dt_col = next((c for c in dt_candidates if c in aqhi_df.columns), None)
    if dt_col is None:
        raise ValueError("No datetime column (ReadingDate/DateTime) in last6h.csv")

    aqhi_df[dt_col] = pd.to_datetime(aqhi_df[dt_col], errors="coerce")

    # Keep only rows with valid coords and AQHI candidate
    aqhi_df = aqhi_df.dropna(subset=[lat_col, lon_col, aqhi_col])

    # ---- 5) For each station, keep the last AQHI record ----
    aqhi_df = aqhi_df.sort_values(dt_col).drop_duplicates(stn_col, keep="last")

    # ---- 6) Build GeoDataFrame in WGS84, then to TARGET_CRS ----
    gdf = gpd.GeoDataFrame(
        aqhi_df,
        geometry=gpd.points_from_xy(aqhi_df[lon_col], aqhi_df[lat_col]),
        crs="EPSG:4326"   # lat/lon in the CSV
    ).to_crs(TARGET_CRS)

    # ---- 7) Numeric AQHI for IDW, capped at 10 ----
    gdf["aqhi_val"] = pd.to_numeric(gdf[aqhi_col], errors="coerce").clip(upper=10)
    gdf["weight"]   = STATION_WEIGHT  # 1.0

    return gdf



def load_purple_from_url() -> gpd.GeoDataFrame:
    """
    Load PurpleAir sensors from AB_PM25_map.json (non-OGR JSON),
    build a GeoDataFrame, and compute estimated AQHI.
    """
    # --- 1) Fetch the JSON ---
    with urllib.request.urlopen(PURPLE_URL) as f:
        raw = f.read().decode("utf-8")
    data = json.loads(raw)

    # data could be:
    #  - a list of dicts
    #  - a dict with "data" key
    #  - a GeoJSON-like dict with "features"
    if isinstance(data, dict) and "features" in data:
        # It's GeoJSON-ish
        gdf = gpd.GeoDataFrame.from_features(data["features"], crs="EPSG:4326")
    else:
        # Try to get to a plain list of records
        if isinstance(data, dict) and "data" in data:
            records = data["data"]
        else:
            records = data

        if not isinstance(records, list):
            raise ValueError("Unexpected JSON structure in AB_PM25_map.json")

        df = pd.DataFrame(records)

        # --- 2) Identify lat/lon columns ---
        lat_candidates = ["lat", "Lat", "latitude", "Latitude", "LAT"]
        lon_candidates = ["lon", "Lon", "lng", "Lng", "longitude", "Longitude", "LON"]

        lat_col = next((c for c in lat_candidates if c in df.columns), None)
        lon_col = next((c for c in lon_candidates if c in df.columns), None)

        if lat_col is None or lon_col is None:
            raise ValueError("Could not find lat/lon columns in AB_PM25_map.json")

        # --- 3) Build GeoDataFrame from lat/lon ---
        gdf = gpd.GeoDataFrame(
            df,
            geometry=[Point(xy) for xy in zip(df[lon_col], df[lat_col])],
            crs="EPSG:4326"
        )

    # --- 4) Reproject to target CRS ---
    gdf = gdf.to_crs(TARGET_CRS)

    # --- 5) Find PM2.5 column and compute estimated AQHI ---
    pm_candidates = ["pm_corr", "PM2_5", "pm25", "pm_25", "pm2_5", "pm"]
    pm_col = next((c for c in pm_candidates if c in gdf.columns), None)
    if pm_col is None:
        raise ValueError("No PM2.5 column found in AB_PM25_map.json")

    gdf["pm_val"] = pd.to_numeric(gdf[pm_col], errors="coerce")

    def est_aqhi(pm):
        if pd.isna(pm):
            return np.nan
        aq = round(pm / 10.0) + 1
        return min(max(aq, 0), 10)

    gdf["aqhi_val"] = gdf["pm_val"].apply(est_aqhi)
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

    print(f"[grid] {shp_path} – stations total: {len(stations_gdf)}, non-null AQHI: {len(stn)}")

    # 3) Prepare blended points (stations + PurpleAir)
    pa  = purple_gdf.copy()
    pa  = pa[pa["aqhi_val"].notna()].copy()

    print(f"[grid] {shp_path} – PurpleAir total: {len(purple_gdf)}, non-null AQHI: {len(pa)}")

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
        from pathlib import Path
        outfile_geojson = Path(outfile_geojson)

        # Reproject grid to WGS84 for Leaflet
        grid_out = grid.to_crs("EPSG:4326")

        # Save full (stations+blended) grid
        grid_out.to_file(outfile_geojson, driver="GeoJSON")

        # Also save a stations-only file
        stn_only_path = outfile_geojson.with_name(
            outfile_geojson.stem.replace("_station_vs_blended", "_station_only") + ".geojson"
        )

        stn_only_cols = ["geometry", "aqhi_stn", "aqhi_stn_cat"]
        grid_stn_only_out = grid_out[stn_only_cols].copy()
        grid_stn_only_out.to_file(stn_only_path, driver="GeoJSON")

    return grid
                            

if __name__ == "__main__":
    # Load point data once
    stations_gdf = load_stations_from_url()
    purple_gdf   = load_purple_from_url()


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
