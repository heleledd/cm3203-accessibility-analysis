import geopandas as gpd

# Load files
schools = gpd.read_file("data/script_output/grid_cells_with_school_distances.geojson")
parks = gpd.read_file("data/script_output/grid_cells_with_park_distances.geojson")
hospitals = gpd.read_file("data/script_output/grid_cells_with_hospital_distances.geojson")

# Keep only needed columns (avoid duplicate geometry columns)
parks = parks[["grid_id", "nearest_park"]]
hospitals = hospitals[["grid_id", "nearest_hospital"]]

# Merge by grid_id
combined = schools.merge(parks, on="grid_id", how="left")
combined = combined.merge(hospitals, on="grid_id", how="left")
combined = combined.to_crs(epsg=4326)  # reproject to WGS84

# Save result
combined.to_file("grid_cells_combined_distance.geojson", driver="GeoJSON")
