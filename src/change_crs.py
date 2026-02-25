import geopandas as gpd

edges = gpd.read_file("edges.geojson")
print(edges.crs)  # check what CRS it's currently in

edges = edges.to_crs(epsg=4326)  # reproject to WGS84
edges.to_file("edges_wgs84.geojson", driver="GeoJSON")
