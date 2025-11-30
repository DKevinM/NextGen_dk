# generate_npri_data.py
import requests
import json

def fetch_and_save_npri():
    url = "https://maps-cartes.ec.gc.ca/arcgis/rest/services/STB_DGST/NPRI/MapServer/0/query?where=1=1&outFields=*&returnGeometry=true&f=geojson&outSR=4326"
    
    response = requests.get(url)
    data = response.json()
    
    # Save the GeoJSON file
    with open('data/npri_data.geojson', 'w') as f:
        json.dump(data, f)
    
    print("Data saved to data/npri_data.geojson")

if __name__ == "__main__":
    fetch_and_save_npri()
