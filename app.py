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
