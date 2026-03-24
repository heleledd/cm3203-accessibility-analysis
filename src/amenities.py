import logging

import osmnx as ox
import geopandas as gpd

from config import TARGET_CRS


def get_park_data(network, access_points_path, boundaries_path):
    """Loads parks access points and boundaries from local .geojson file"""
    logging.info(f"Loading park data from {access_points_path} and {boundaries_path}...")
    
    # read in the access points
    park_access_points_gdf = gpd.read_file(access_points_path)

    # drop any rows where the geometry is null/missing
    park_access_points_gdf = park_access_points_gdf.dropna(subset=['geometry'])

    # also drop "empty" geometries which can crash the CRS transformer
    park_access_points_gdf = park_access_points_gdf[~park_access_points_gdf.geometry.is_empty]

    # convert to the target CRS
    park_access_points_gdf = park_access_points_gdf.to_crs(TARGET_CRS)

    # calculate nearest node in the street network to each park entrance and store in dataframe
    park_access_points_gdf['nearest_node'] = network.get_node_ids(
        park_access_points_gdf.geometry.x, 
        park_access_points_gdf.geometry.y
    )

    # read in the boundaries and save them in a geopandas dataframe
    park_boundaries_gdf = gpd.read_file(boundaries_path).to_crs(TARGET_CRS)

    return park_access_points_gdf, park_boundaries_gdf


def get_osm_features(network, bbox, tags):
    """Downloads point/polygon features from OSM and maps their centroids to network nodes."""
    logging.info(f"Downloading OSM features for tags: {tags}...")
    
    try:
        features_gdf = ox.features_from_bbox(bbox, tags=tags).to_crs(TARGET_CRS)
        features_gdf['centroid'] = features_gdf.geometry.centroid
        
    
        # snap centroids to the network
        features_gdf['nearest_node'] = network.get_node_ids(
            features_gdf['centroid'].x, 
            features_gdf['centroid'].y
        )

        return features_gdf
    
    except ox._errors.InsufficientResponseError:
        # Catch the specific OSMnx error when no features exist
        logging.warning(f"No features found for tags: {tags}. Skipping...")
        
        # Return an empty GeoDataFrame so the pipeline doesn't break
        return gpd.GeoDataFrame(crs=TARGET_CRS)
        
    except Exception as e:
        # Catch any other random network drops or timeout errors
        logging.error(f"An unexpected error occurred fetching {tags}: {e}")
        return gpd.GeoDataFrame(crs=TARGET_CRS)
    
    