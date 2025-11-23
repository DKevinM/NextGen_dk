from flask import Flask, request, jsonify
from datetime import datetime, timedelta
import math, random, os
from collections import Counter
import requests

app = Flask(__name__)

EARTH_RADIUS_M = 6371000.0
N_PARTICLES_DEFAULT = 15
HOURS_BACK_DEFAULT = 4.0
DT_MINUTES = 5
FOOTPRINT_DLAT = 0.05
FOOTPRINT_DLON = 0.05

def parse_iso8601(ts_str: str) -> datetime:
    if ts_str.endswith("Z"):
        ts_str = ts_str[:-1]
    return datetime.fromisoformat(ts_str)

def uv_to_latlon_step(lat_deg, u_ms, v_ms, dt_seconds):
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

# ── HRDPS wind via GeoMet WFS ─────────────────────────────────────

def get_wind_at(lat, lon, when):
    base = "https://geo.weather.gc.ca/geomet"
    bbox = f"{lon},{lat},{lon},{lat}"

    # U component
    params_u = {
        "service": "WFS",
        "version": "2.0.0",
        "request": "GetFeature",
        "typeName": "HRDPS.CONTINENTAL_UU",
        "bbox": bbox,
        "outputFormat": "application/json",
    }
    u_json = requests.get(base, params=params_u).json()
    u_val = 0.0
    try:
        u_val = u_json["features"][0]["properties"]["HRDPS.CONTINENTAL_UU"]
    except Exception:
        pass

    # V component
    params_v = {
        "service": "WFS",
        "version": "2.0.0",
        "request": "GetFeature",
        "typeName": "HRDPS.CONTINENTAL_VV",
        "bbox": bbox,
        "outputFormat": "application/json",
    }
    v_json = requests.get(base, params=params_v).json()
    v_val = 0.0
    try:
        v_val = v_json["features"][0]["properties"]["HRDPS.CONTINENTAL_VV"]
    except Exception:
        pass

    return u_val, v_val

# ── Trajectory core ───────────────────────────────────────────────

def compute_ensemble_trajectories(start_lat, start_lon, start_time,
                                  n_particles=N_PARTICLES_DEFAULT,
                                  hours_back=HOURS_BACK_DEFAULT,
                                  dt_minutes=DT_MINUTES):
    dt_seconds = dt_minutes * 60.0
    n_steps = int(hours_back * 60.0 / dt_minutes)
    trajectories = []

    for _ in range(n_particles):
        lat = start_lat + random.uniform(-0.02, 0.02)
        lon = start_lon + random.uniform(-0.02, 0.02)
        t = start_time
        traj = [(lat, lon, t)]

        for _ in range(n_steps):
            u, v = get_wind_at(lat, lon, t)
            dlat, dlon = uv_to_latlon_step(lat, u, v, dt_seconds)
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
        bearings.append(bearing_deg(lat_start, lon_start, lat_end, lon_end))

    if not bearings:
        return {"sector": None, "mean_bearing": None, "samples": 0}

    mean_bearing = sum(bearings) / len(bearings)
    sectors = [bearing_to_compass_sector(b) for b in bearings]
    from collections import Counter
    sector_counts = Counter(sectors)
    dom, _ = sector_counts.most_common(1)[0]

    return {"sector": dom, "mean_bearing": mean_bearing, "samples": len(bearings)}

def trajectories_to_geojson_lines(trajectories):
    features = []
    for i, traj in enumerate(trajectories):
        coords = [[lon, lat] for (lat, lon, _) in traj]
        features.append({
            "type": "Feature",
            "properties": {"id": i},
            "geometry": {"type": "LineString", "coordinates": coords},
        })
    return {"type": "FeatureCollection", "features": features}

# ── API route ─────────────────────────────────────────────────────

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

    return jsonify({
        "trajectories": traj_geojson,
        "footprint": footprint,
        "dominant_direction": dom_dir,
        "receptor": {
            "lat": lat,
            "lon": lon,
            "time": start_time.isoformat() + "Z",
        },
    })

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5001))
    app.run(host="0.0.0.0", port=port)
