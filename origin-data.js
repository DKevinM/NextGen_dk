// origin-data.js
// Minimal station + PurpleAir + NPRI loaders for the origin map

// ----------------- Helpers -----------------
function arrayToFeatureCollection(arr){
  if (!arr || !arr.length) return { type:'FeatureCollection', features: [] };
  const sample = arr[0] || {};
  const latK=['lat','latitude','Latitude','LAT','Lat'].find(k=>k in sample) || 'lat';
  const lonK=['lon','lng','long','longitude','Longitude','LON','Lng'].find(k=>k in sample) || 'lon';
  return {
    type:'FeatureCollection',
    features: arr.map(r => ({
      type:'Feature',
      properties: { ...r },
      geometry: {
        type:'Point',
        coordinates: [ +r[lonK], +r[latK] ]
      }
    })).filter(f =>
      Number.isFinite(f.geometry.coordinates[0]) &&
      Number.isFinite(f.geometry.coordinates[1])
    )
  };
}

function toPointFC(json) {
  const arr = Array.isArray(json)
    ? json
    : (json.data || json.features || Object.values(json || {}));

  return arrayToFeatureCollection(arr || []);
}

// ----------------- 1) AQHI stations from last6h.csv -----------------

const unitsLookup = {
  "AQHI": " ", "Ozone": " ppb", "Total Oxides of Nitrogen": " ppb",
  "Hydrogen Sulphide": " ppb", "Total Reduced Sulphur": " ppb", "Sulphur Dioxide": " ppb",
  "Fine Particulate Matter": " µg/m³", "Total Hydrocarbons": " ppm", "Carbon Monoxide": " ppm",
  "Wind Direction": " degrees", "Relative Humidity": " %", "Outdoor Temperature": " °C",
  "Nitric Oxide": " ppb", "Wind Speed": " km/hr", "Non-methane Hydrocarbons": " ppm",
  "Nitrogen Dioxide": " ppb", "Methane": " ppm"
};

const abbrLookup = {
  "AQHI": "AQHI", "Ozone": "O₃", "Total Oxides of Nitrogen": "NOx",
  "Hydrogen Sulphide": "H₂S", "Total Reduced Sulphur": "TRS", "Sulphur Dioxide": "SO₂",
  "Fine Particulate Matter": "PM2.5", "Total Hydrocarbons": "THC", "Carbon Monoxide": "CO",
  "Wind Direction": "wd", "Relative Humidity": "RH", "Outdoor Temperature": "ET",
  "Nitric Oxide": "NO", "Wind Speed": "ws", "Non-methane Hydrocarbons": "NMHC",
  "Nitrogen Dioxide": "NO₂", "Methane": "CH₄"
};

const shortLookup = {
  "AQHI": "AQHI", "Ozone": "O3", "Total Oxides of Nitrogen": "NOX",
  "Hydrogen Sulphide": "H2S", "Total Reduced Sulphur": "TRS", "Sulphur Dioxide": "SO2",
  "Fine Particulate Matter": "PM2.5", "Total Hydrocarbons": "THC", "Carbon Monoxide": "CO",
  "Wind Direction": "wd", "Relative Humidity": "RH", "Outdoor Temperature": "ET",
  "Nitric Oxide": "NO", "Wind Speed": "ws", "Non-methane Hydrocarbons": "NMHC",
  "Nitrogen Dioxide": "NO2", "Methane": "CH4"
};

let dataByStation = {};
window.dataByStation = dataByStation;

window.dataReady = fetch('https://raw.githubusercontent.com/DKevinM/AB_datapull/main/data/last6h.csv')
  .then(res => res.text())
  .then(text => {
    const rows = text.trim().split('\n');
    const headers = rows.shift().split(',');

    const raw = {};
    rows.forEach(line => {
      const cols = line.split(',');
      const e = Object.fromEntries(headers.map((h,i)=>[h,cols[i]]));

      if (!e.Latitude||!e.Longitude||isNaN(e.Latitude)||isNaN(e.Longitude)) return;

      e.ParameterName = e.ParameterName||"AQHI";
      e.Units = unitsLookup[e.ParameterName]||"";
      e.Abbreviation = abbrLookup[e.ParameterName]||"";
      e.Shortform = shortLookup[e.ParameterName]||"";

      let v = parseFloat(e.Value);
      if (["Ozone","Total Oxides of Nitrogen","Hydrogen Sulphide","Total Reduced Sulphur","Sulphur Dioxide","Nitric Oxide","Nitrogen Dioxide"].includes(e.ParameterName)) {
        v *= 1000;
      }
      if (isNaN(v)) return;
      e.Value = v;

      const utc = new Date(e.ReadingDate);
      e.DisplayDate = utc.toLocaleString("en-CA", {
        timeZone: "America/Edmonton",
        hour12: true
      });
      e.ReadingDate = utc.toISOString();

      raw[e.StationName] = raw[e.StationName] || [];
      raw[e.StationName].push(e);
    });

    Object.entries(raw).forEach(([station, arr]) => {
      arr.sort((a, b) => new Date(b.ReadingDate) - new Date(a.ReadingDate));
      const byParam = {};
      arr.forEach(e => {
        const param = e.ParameterName;
        if (!byParam[param] || new Date(e.ReadingDate) > new Date(byParam[param].ReadingDate)) {
          byParam[param] = e;
        }
      });
      dataByStation[station] = Object.values(byParam);
    });
  });

window.fetchAllStationData = function () {
  const stationNames = Object.keys(dataByStation);
  if (!stationNames.length) return Promise.resolve([]);

  const orderedParams = [
    "AQHI", "Outdoor Temperature", "Relative Humidity", "Wind Speed", "Wind Direction", 
    "Nitrogen Dioxide", "Total Oxides of Nitrogen", "Nitric Oxide", "Ozone",
    "Fine Particulate Matter", "Sulphur Dioxide", "Hydrogen Sulphide", "Total Reduced Sulphur",
    "Carbon Monoxide", "Total Hydrocarbons", "Methane", "Non-methane Hydrocarbons"  
  ];

  const shortformOverride = {
    "Outdoor Temperature": "Temp",
    "Relative Humidity": "Humidity",
    "Wind Speed": "Wind Speed",
    "Wind Direction": "Wind Dir"
  };

  const results = stationNames.map(stationName => {
    const stationData = dataByStation[stationName];
    if (!stationData || !stationData.length) return null;

    const paramLookup = {};
    let latestTime = null;
    for (const r of stationData) {
      paramLookup[r.ParameterName] = r;
      const t = new Date(r.ReadingDate);
      if (!latestTime || t > latestTime) latestTime = t;
    }

    const displayTime = latestTime
      ? latestTime.toLocaleString("en-CA", { timeZone: "America/Edmonton", hour12: true })
      : "Invalid Date";

    const lines = orderedParams
      .filter(p => paramLookup[p] && p !== "AQHI")
      .map(p => {
        const r = paramLookup[p];
        const label = shortformOverride[p] || r.Shortform || p;
        const value = r.Value;
        const unit = r.Units || "";
        return `${label}: ${value}${unit}`;
      });

    const aqhiValue = paramLookup["AQHI"]?.Value || "N/A";
    const lat = stationData[0].Latitude;
    const lon = stationData[0].Longitude;

    return {
      stationName,
      lat: +lat,
      lon: +lon,
      aqhi: aqhiValue,
      html: `
        <div style="font-size:0.9em;">
          <strong>${stationName}</strong><br>
          <small><em>${displayTime}</em></small><br>
          AQHI: ${aqhiValue > 10 ? "10+" : aqhiValue}<br>
          ${lines.join("<br>")}
        </div>
      `
    };
  }).filter(Boolean);

  return Promise.resolve(results);
};


// Global FeatureCollections
window.STATIONS_FC = { type: "FeatureCollection", features: [] };
window.PURPLE_FC   = { type: "FeatureCollection", features: [] };
window.NPRI_FC     = { type: "FeatureCollection", features: [] };

// 1) Stations
window.stationsFCReady = (async () => {
  try {
    await window.dataReady;               // your last6h.csv loader
    const rows = await window.fetchAllStationData();

    window.STATIONS_FC = {
      type: "FeatureCollection",
      features: (rows || [])
        .filter(r => isFinite(+r.lat) && isFinite(+r.lon))
        .map(r => ({
          type: "Feature",
          properties: { ...r },
          geometry: { type: "Point", coordinates: [ +r.lon, +r.lat ] }
        }))
    };
    console.log("[origin] STATIONS_FC features:", window.STATIONS_FC.features.length);
  } catch (e) {
    console.error("[origin] stationsFCReady failed", e);
    window.STATIONS_FC = { type: "FeatureCollection", features: [] };
  }
})();

// 2) PurpleAir
window.purpleFCReady = (async () => {
  try {
    const res  = await fetch("https://raw.githubusercontent.com/DKevinM/AB_datapull/main/data/ACA_PM25_map.json");
    if (!res.ok) throw new Error(`HTTP ${res.status}`);

    const json = await res.json();

    // Use the flexible helper so it works whether the file is an array,
    // an object with .data/.features, and whether it uses lat/Lat/Latitude, etc.
    window.PURPLE_FC = toPointFC(json);

    console.log("[origin] PURPLE_FC features:", window.PURPLE_FC.features.length);
  } catch (e) {
    console.error("[origin] purpleFCReady failed", e);
    window.PURPLE_FC = { type: "FeatureCollection", features: [] };
  }
})();


// 3) NPRI
window.npriFCReady = (async () => {
  try {
    const url =
      "https://maps-cartes.ec.gc.ca/arcgis/rest/services/STB_DGST/NPRI/MapServer/0/query" +
      "?where=1=1&outFields=*&returnGeometry=true&outSR=4326&f=pjson";

    console.log("[origin] NPRI GET", url);
    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const json = await res.json();

    const feats = Array.isArray(json.features) ? json.features : [];

    // Convert Esri JSON { attributes, geometry: {x,y} } → GeoJSON Point
    window.NPRI_FC = {
      type: "FeatureCollection",
      features: feats
        .filter(f => f && f.geometry && typeof f.geometry.x === "number" && typeof f.geometry.y === "number")
        .map(f => ({
          type: "Feature",
          properties: { ...(f.attributes || {}) },
          geometry: {
            type: "Point",
            coordinates: [ f.geometry.x, f.geometry.y ]
          }
        }))
    };

    console.log("[origin] NPRI_FC features:", window.NPRI_FC.features.length);
  } catch (e) {
    console.error("[origin] npriFCReady failed", e);
    window.NPRI_FC = { type: "FeatureCollection", features: [] };
  }
})();





// Once NPRI_FC is ready, plot all NPRI facilities on the map
window.npriFCReady
  ?.then(() => {
    console.log("[origin] plotting NPRI facilities:", NPRI_FC.features.length);
    npriLayerGroup.clearLayers();

    (NPRI_FC.features || []).forEach(f => {
      const ll = getFeatureLatLon(f);
      if (!ll) return;

      const p   = f.properties || {};
      const fac = p.FACILITY_NAME || p.FacilityName || p.facility || "Facility";
      const co  = p.COMPANY_NAME  || p.Company      || p.company  || "";
      const yr  = p.REPORTING_YEAR || p.ReportingYear || p.year || "";
      const label = co ? `${fac} (${co})` : fac;

      const popupHtml = `
        <b>${label}</b><br>
        ${yr ? "Reporting year: " + yr + "<br>" : ""}
        <small>Approx. location: ${ll.lat.toFixed(4)}, ${ll.lon.toFixed(4)}</small>
      `;

      L.circleMarker([ll.lat, ll.lon], {
        radius: 4,
        color: "#800026",
        weight: 1,
        fillOpacity: 0.7
      })
      .bindPopup(popupHtml)
      .addTo(npriLayerGroup);
    });
  })
  .catch(err => {
    console.error("[origin] error plotting NPRI facilities", err);
  });



