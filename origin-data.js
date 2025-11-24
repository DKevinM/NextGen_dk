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

// Build a FeatureCollection for stations
window.stationsFCReady = (async () => {
  try {
    await window.dataReady;
    const rows = await window.fetchAllStationData();
    return {
      type: 'FeatureCollection',
      features: (rows || []).map(r => ({
        type: 'Feature',
        properties: { ...r, sourceType: 'station' },
        geometry: { type: 'Point', coordinates: [r.lon, r.lat] }
      }))
    };
  } catch (e) {
    console.error('[stationsFCReady] failed', e);
    return { type: 'FeatureCollection', features: [] };
  }
})();

// ----------------- 2) PurpleAir -----------------

window.purpleFCReady = (async () => {
  try {
    const res = await fetch('https://raw.githubusercontent.com/DKevinM/AB_datapull/main/data/ACA_PM25_map.json');
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const json = await res.json();
    const fc = toPointFC(json);
    fc.features.forEach(f => {
      f.properties.sourceType = 'purpleair';
    });
    return fc;
  } catch (e) {
    console.error('[purpleFCReady] failed', e);
    return { type:'FeatureCollection', features: [] };
  }
})();

// ----------------- 3) NPRI via Esri Leaflet -----------------

const NPRI_REST_URL = 'https://maps-cartes.ec.gc.ca/arcgis/rest/services/STB_DGST/NPRI/MapServer';

window.npriFCReady = (async () => {
  try {
    // requires esri-leaflet to be loaded in origin.html
    const fc = await new Promise((resolve, reject) => {
      L.esri.query({ url: `${NPRI_REST_URL}/0` })
        .where('1=1')
        .fields(['*'])
        .returnGeometry(true)
        .run((err, geojson) => err ? reject(err) : resolve(geojson));
    });
    fc.features.forEach(f => {
      f.properties.sourceType = 'npri';
    });
    return fc;
  } catch (e) {
    console.error('[npriFCReady] failed', e);
    return { type:'FeatureCollection', features: [] };
  }
})();
