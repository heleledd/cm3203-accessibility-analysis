# python .\grid_to_chloropleth.py --input "../data/grid_cells_with_hospital_distances_20260216_100847.geojson"

import geopandas as gpd
import folium
import argparse
import json

def main(geojson_data):
    # Load GeoJSON and set the CRS
    gdf = gpd.GeoDataFrame.from_features(geojson_data['features'])
    gdf.set_crs('EPSG:27700', inplace=True)
    
    # Reproject to WGS84 (lat/lon)
    gdf_wgs84 = gdf.to_crs('EPSG:4326')
    
    # Convert back to GeoJSON dict
    geojson_wgs84 = json.loads(gdf_wgs84.to_json())
    
    # Now use with Folium
    m = folium.Map(
        location=(51.49611, -3.17332),
        zoom_start=13
    )
    
    folium.Choropleth(
        geo_data=geojson_wgs84,
        data=gdf_wgs84,
        columns=['grid_id', 'nearest_hospital'],
        key_on='feature.properties.grid_id',
        fill_color='RdYlGn_r',
        fill_opacity=0.7,
        line_opacity=0.2,
        legend_name='Distance to Nearest Hospital (m)'
    ).add_to(m)
    
    m.save('hospital_chloropleth.html')





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a heatmap of accessibility to supermarkets in Cardiff.')
    parser.add_argument('--input', type=str, default='output/geojson/supermarket_access.geojson', help='Path to the input GeoJSON file')
    
    args = parser.parse_args()
    
    # Load the GeoJSON data
    with open(args.input) as f:
        geojson_data = json.load(f)
    

    # Call main function with loaded data
    main(geojson_data)
