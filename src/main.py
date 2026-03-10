import os
from tqdm import tqdm
from datetime import datetime
import logging
import math
import concurrent.futures
import osmnx as ox
import geopandas as gpd

from config import (
    CARDIFF_BBOX,
    GRID_SIZE_METERS,
    OUTPUT_DATA_DIR,
    PARK_BOUNDARY_PATH,
    PARK_ACCESS_POINTS_PATH,
    NUM_WORKERS
)

from network_and_bbox import get_street_network_graph, reproject_bbox, split_bbox_into_grid
from amenities import get_park_data, get_osm_features


def main(
        bbox=CARDIFF_BBOX,
        grid_size=GRID_SIZE_METERS,
        park_access_points_path=PARK_ACCESS_POINTS_PATH,
        park_boundary_path=PARK_BOUNDARY_PATH,
        num_workers=NUM_WORKERS
    ):
    
    # load street network
    logging.info("Loading Pandana network...")
    network = get_street_network_graph(bbox)
    
    # split the bounding box into a grid
    bbox_reprojected = reproject_bbox(bbox)
    grid_cells_gdf = split_bbox_into_grid(bbox_reprojected, grid_size)

    # load amenity features
    park_access_points_gdf, park_boundaries_gdf = get_park_data(network, park_access_points_path, park_boundary_path)
    gps_gdf = get_osm_features(network, bbox, tags={"amenity": "doctors"})
    schools_gdf = get_osm_features(network, bbox, tags={"amenity": "school"})

    MAX_DIST = 5000 # 5km max search distance
    logging.info("Attaching amenities to the network...")

    network.set_pois('parks', MAX_DIST, 1, park_access_points_gdf.geometry.x, park_access_points_gdf.geometry.y)
    network.set_pois('gps', MAX_DIST, 1, gps_gdf['centroid'].x, gps_gdf['centroid'].y)
    network.set_pois('schools', MAX_DIST, 1, schools_gdf['centroid'].x, schools_gdf['centroid'].y)

    # map grid cells to network nodes
    logging.info("Mapping grid cells to network nodes...")
    centroids = grid_cells_gdf.geometry.centroid
    
    # Snapping grid cells using Pandana's vectorized method
    grid_cells_gdf['nearest_node'] = network.get_node_ids(centroids.x, centroids.y)

    # calculate shortest routes for ALL grid cells instantly
    logging.info("Calculating shortest routes with Pandana...")
    park_dists = network.nearest_pois(MAX_DIST, 'parks', num_pois=1)
    gp_dists = network.nearest_pois(MAX_DIST, 'gps', num_pois=1)
    school_dists = network.nearest_pois(MAX_DIST, 'schools', num_pois=1)

    # Map the distances back to the dataframe (column 1 represents the distance to the 1st nearest POI)
    grid_cells_gdf['nearest_park'] = grid_cells_gdf['nearest_node'].map(park_dists[1])
    grid_cells_gdf['nearest_gp'] = grid_cells_gdf['nearest_node'].map(gp_dists[1])
    grid_cells_gdf['nearest_school'] = grid_cells_gdf['nearest_node'].map(school_dists[1])

    # Check if any centroids are INSIDE a polygon e.g park or school
    logging.info("Applying 0.0m distance for cells inside amenity boundaries...")
    centroids_gdf = gpd.GeoDataFrame(geometry=centroids, crs=grid_cells_gdf.crs)

    # Parks
    parks_intersect = gpd.sjoin(centroids_gdf, park_boundaries_gdf, how='inner', predicate='intersects')
    grid_cells_gdf.loc[parks_intersect.index, 'nearest_park'] = 0.0

    # Schools
    schools_intersect = gpd.sjoin(centroids_gdf, schools_gdf, how='inner', predicate='intersects')
    grid_cells_gdf.loc[schools_intersect.index, 'nearest_school'] = 0.0

    logging.info("Saving results...")
    logging.getLogger().setLevel(logging.DEBUG)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_output_dir = os.path.join(OUTPUT_DATA_DIR, timestamp)
    os.makedirs(run_output_dir, exist_ok=True)

    grid_cells_gdf.to_file(os.path.join(run_output_dir, 'grid_cells_accessibility.geojson'), driver='GeoJSON')
    
    # Save features used with the distance file in the same folder
    park_boundaries_gdf.to_file(os.path.join(run_output_dir, 'park.geojson'), driver='GeoJSON')
    gps_gdf.drop(columns=['centroid'], errors='ignore').to_file(os.path.join(run_output_dir, 'gp.geojson'), driver='GeoJSON')
    schools_gdf.drop(columns=['centroid'], errors='ignore').to_file(os.path.join(run_output_dir, 'school.geojson'), driver='GeoJSON')

    print(f"\nSuccess! Results saved to {run_output_dir}")



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
