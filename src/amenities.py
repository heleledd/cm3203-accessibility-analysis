import logging

import osmnx as ox
import geopandas as gpd

from config import TARGET_CRS

# !!! need to change this because of new data which has access points!!
def get_park_boundary_nodes(G, park_path):
    """Loads parks from local GeoJSON and maps their boundaries to network nodes."""
    park_gdf = gpd.read_file(park_path).to_crs(TARGET_CRS)
    park_gdf['boundary'] = park_gdf.geometry.boundary

    return park_gdf

def get_osm_features(G, bbox, tags):
    """Downloads point/polygon features from OSM and maps their centroids to network nodes."""
    logging.info(f"Downloading OSM features for tags: {tags}...")
    features_gdf = ox.features_from_bbox(bbox, tags=tags).to_crs(TARGET_CRS)
    features_gdf['centroid'] = features_gdf.geometry.centroid
    features_gdf['nearest_node'] = ox.distance.nearest_nodes(
        G, features_gdf['centroid'].x, features_gdf['centroid'].y
    )
    return features_gdf