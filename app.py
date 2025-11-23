# app.py

from flask import Flask, request, jsonify, send_from_directory
from datetime import datetime, timedelta
import math
import random
from collections import Counter
import os

app = Flask(__name__, static_folder="static", static_url_path="")

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────

N_PARTICLES_DEFAULT = 15
HOURS_BACK_DEFAULT = 4.0          # total hours back for trajectories
DT_MINUTES = 5                    # timestep in minutes
EARTH_RADIUS_M = 6371000.0

FOOTPRINT_DLAT = 0.05
FOOTPRINT_DLON = 0.05


# ─────────────────────────────────────────────
# UTILITIES
# ─────────────────────────────────────────────

def parse_iso8601(ts_str: str) -> datetime:
    if ts_str.endswith("Z"):
        ts_str = ts_str[:-1]
    return datetime.fromisoformat(ts_str)


def uv_to_latlon_step(lat_deg: float, u_ms: float, v_ms: float, dt_seconds: float):
    lat_rad = math.radians(lat_deg)

    dlat = (v_ms * dt_seconds) / EARTH_RADIUS_M * (180.0 / math.pi)
    dlon = (u_ms * dt_seconds) / (EARTH_RADIUS_M * math.cos(lat_rad)) * (180.0 / math.pi)

    return dlat, dlon


def bearing_deg(lat1, lon1, lat2, lon2):
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dlon = math.radians(lon2 - lon1)

    x = math.sin(dlon) * math.cos(phi2)
    y = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(dlon)
    brng = math.degrees(math.atan2(x, y))
    return (brng + 360.0) % 360.0


def bearing_to_compass_sector(bearing, n_sectors=8):
    labels_8 = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
    sector_size = 360.0 / n_sectors
    idx = int((bearing + sector_size / 2.0) // sector_size) % n_sectors
    return labels_8[idx]


# ─────────────────────────────────────────────
# WIND FIELD ACCESS – PLUG YOUR HRDPS/GFS HERE
# ─────────────────────────────────────────────

# Idea: you adapt your existing script to expose ONE function:
# get_wind_at(lat, lon, when) -> (u_ms, v_ms)
#
# Under the hood, you:
#   - open your HRDPS/GFS GRIB/NetCDF once (global state or cached)
#   - find nearest time slice to `when`
#   - bilinear interpolate u and v at (lat, lon)

def get_wind_at(lat, lon, when: datetime):
    """
    Placeholder version.

    Replace this with your HRDPS/GFS reader using cfgrib/pygrib/xarray:

      1. Open GRIB / NetCDF file (or cache it).
      2. Select nearest time step to `when`.
      3. Bilinear interpolate u, v at (lat, lon).

    For now we just use a synthetic westerly wind with noise.
    """
    base_u = 8.0  # m/s eastward
    base_v = 1.0  # m/s northward
    u = base_u + random.uniform(-1.5, 1.5)
    v = base_v + random.uniform(-1.0, 1.0)
    return u, v


# ─────────────────────────────────────────────
# CORE TRAJECTORY / FOOTPRINT LOGIC
# ─────────────────────────────────────────────

def compute_ensemble_trajectories(
    start_lat,
    start_lon,
    start_time: datetime,
    n_particles=N_PARTICLES_DEFAULT,
    hours_back=HOURS_BACK_DEFAULT,
    dt_minutes=DT_MINUTES,
):
    dt_seconds = dt_minutes * 60.0
    n_steps = int(hours_back * 60.0 / dt_minutes)

    trajectories = []

    for p in range(n_particles):
        lat = start_lat + random.uniform(-0.02, 0.02)
        lon = start_lon + random.uniform(-0.02, 0.02)
        t = start_time

        traj = [(lat, lon, t)]

        for _ in range(n_steps):
            u, v = get_wind_at(lat, lon, t)
            dlat, dlon = uv_to_latlon_step(lat, u, v, dt_seconds)

            # backward trajectory: move opposite wind
            lat -= dlat
            lon -= dlon
            t -= timedelta(seconds=dt_seconds)

            traj.append((lat, lon, t))

        trajectories.append(traj)

    return trajectories


def build_footprint_grid(trajectories, dlat=FOOTPRINT_DLAT, dlon=FOOTPRINT_DLON):
    counts = {}

    for traj in trajectories:
        for (lat, lon, _) in traj:
            iy = int(round(lat / dlat))
            ix = int(round(lon / dlon))
            key = (iy, ix)
            counts[key] = counts.get(key, 0) + 1

    footprint = []
    for (iy, ix), count in counts.items():
        lat_center = iy * dlat
        lon_center = ix * dlon
        footprint.append({"lat": lat_center, "lon": lon_center, "count": count})

    return footprint


def compute_dominant_direction(trajectories):
    bearings = []

    for traj in trajectories:
        if len(traj) < 2:
            continue
        lat_start, lon_start, _ = traj[-1]
        lat_end, lon_end, _ = traj[0]
        brng = bearing_deg(lat_start, lon_start, lat_end, lon_end)
        bearings.append(brng)

    if not bearings:
        return {"sector": None, "mean_bearing": None, "samples": 0}

    mean_bearing = sum(bearings) / len(bearings)
    sectors = [bearing_to_compass_sector(b) for b in bearings]
    sector_counts = Counter(sectors)
    dominant_sector, _ = sector_counts.most_common(1)[0]

    return {
        "sector": dominant_sector,
        "mean_bearing": mean_bearing,
        "samples": len(bearings),
    }


def trajectories_to_geojson_lines(trajectories):
    features = []
    for i, traj in enumerate(trajectories):
        coords = [[lon, lat] for (lat, lon, _) in traj]
        features.append(
            {
                "type": "Feature",
                "properties": {"id": i},
                "geometry": {"type": "LineString", "coordinates": coords},
            }
        )
    return {"type": "FeatureCollection", "features": features}


# ─────────────────────────────────────────────
# API + STATIC
# ─────────────────────────────────────────────

@app.route("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


@app.route("/api/trajectory")
def api_trajectory():
    try:
        lat = float(request.args.get("lat"))
        lon = float(request.args.get("lon"))
    except (TypeError, ValueError):
        return jsonify({"error": "lat and lon are required as floats"}), 400

    time_str = request.args.get("time")
    if time_str:
        start_time = parse_iso8601(time_str)
    else:
        start_time = datetime.utcnow()

    n_particles = int(request.args.get("n", N_PARTICLES_DEFAULT))
    hours_back = float(request.args.get("hours_back", HOURS_BACK_DEFAULT))

    trajectories = compute_ensemble_trajectories(
        start_lat=lat,
        start_lon=lon,
        start_time=start_time,
        n_particles=n_particles,
        hours_back=hours_back,
    )

    footprint = build_footprint_grid(trajectories)
    dom_dir = compute_dominant_direction(trajectories)
    traj_geojson = trajectories_to_geojson_lines(trajectories)

    response = {
        "trajectories": traj_geojson,
        "footprint": footprint,
        "dominant_direction": dom_dir,
        "receptor": {
            "lat": lat,
            "lon": lon,
            "time": start_time.isoformat() + "Z",
        },
    }

    return jsonify(response)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5001, debug=True)



import requests

def get_wind_at(lat, lon, when):
    """
    Pull HRDPS U and V winds from ECCC GeoMet for the given lat/lon/time.
    """

    base = "https://geo.weather.gc.ca/geomet"
    
    # Format bbox for a single point
    bbox = f"{lon},{lat},{lon},{lat}"

    # 1. Get U wind (eastward)
    params_u = {
        "service": "WFS",
        "version": "2.0.0",
        "request": "GetFeature",
        "typeName": "HRDPS.CONTINENTAL_UU",
        "bbox": bbox,
        "outputFormat": "application/json"
        # you *can* also add: "time": when.isoformat()  (GeoMet tries to pick nearest)
    }
    u_json = requests.get(base, params=params_u).json()
    
    try:
        u_val = u_json["features"][0]["properties"]["HRDPS.CONTINENTAL_UU"]
    except:
        u_val = 0.0  # default fallback

    # 2. Get V wind (northward)
    params_v = {
        "service": "WFS",
        "version": "2.0.0",
        "request": "GetFeature",
        "typeName": "HRDPS.CONTINENTAL_VV",
        "bbox": bbox,
        "outputFormat": "application/json"
    }
    v_json = requests.get(base, params=params_v).json()
    
    try:
        v_val = v_json["features"][0]["properties"]["HRDPS.CONTINENTAL_VV"]
    except:
        v_val = 0.0

    return u_val, v_val
