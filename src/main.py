import os
from datetime import datetime
import logging
import geopandas as gpd
import osmnx as ox

from config import (
    CARDIFF_BBOX,
    GRID_SIZE_METERS,
    OUTPUT_DATA_DIR,
    PARK_BOUNDARY_PATH,
    PARK_ACCESS_POINTS_PATH,
    TARGET_CRS,
    CITY,
    OSM_AMENITIES_CONFIG 
)

from network_and_bbox import get_street_network_graph, reproject_bbox, split_bbox_into_grid
from amenities import get_park_data, get_osm_features


def main(
        bbox=CARDIFF_BBOX,
        grid_size=GRID_SIZE_METERS,
        park_access_points_path=PARK_ACCESS_POINTS_PATH,
        park_boundary_path=PARK_BOUNDARY_PATH,
    ):
    
    # load street network
    logging.info("Loading Pandana network...")
    network = get_street_network_graph(bbox)
    
    # split the bounding box into a grid
    bbox_reprojected = reproject_bbox(bbox)
    grid_cells_gdf = split_bbox_into_grid(bbox_reprojected, grid_size)

    # LOAD ALL AMENITY DATA FROM THE CONFIG FILE
    logging.info("Loading amenity features...")
    loaded_amenities = {}
    
    # Handle custom park data separately, but add it to our unified dictionary
    park_access_points_gdf, park_boundaries_gdf = get_park_data(network, park_access_points_path, park_boundary_path)
    loaded_amenities['parks'] = {
        'gdf': park_boundaries_gdf,
        'check_intersection': True,
        'filename': 'park.geojson',
        'poi_x': park_access_points_gdf.geometry.x,
        'poi_y': park_access_points_gdf.geometry.y
    }

    # Loop through the config to load all OSM amenities dynamically
    for name, config in OSM_AMENITIES_CONFIG.items():
        logging.info(f"Fetching OSM data for {name}...")
        gdf = get_osm_features(network, bbox, tags=config['tags'])
        loaded_amenities[name] = {
            'gdf': gdf,
            'check_intersection': config['check_intersection'],
            'filename': config['filename'],
            'poi_x': gdf['centroid'].x,
            'poi_y': gdf['centroid'].y
        }

    MAX_DIST = 20000 # 20km max search distance
    
    # ATTACH AMENITIES TO NETWORK
    logging.info("Attaching amenities to the network...")
    for name, data in loaded_amenities.items():
        network.set_pois(name, MAX_DIST, 1, data['poi_x'], data['poi_y'])

    # map grid cells to network nodes
    logging.info("Mapping grid cells to network nodes...")
    centroids = grid_cells_gdf.geometry.centroid
    grid_cells_gdf['nearest_node'] = network.get_node_ids(centroids.x, centroids.y)

    # CALCULATE DISTANCES & MAP TO GRID
    logging.info("Calculating shortest routes with Pandana...")
    for name in loaded_amenities.keys():
        dists = network.nearest_pois(MAX_DIST, name, num_pois=1)
        grid_cells_gdf[f'nearest_{name}'] = grid_cells_gdf['nearest_node'].map(dists[1])

    # APPLY 0.0m DISTANCE FOR CELLS INSIDE POLYGONS
    logging.info("Applying 0.0m distance for cells inside amenity boundaries...")
    centroids_gdf = gpd.GeoDataFrame(geometry=centroids, crs=grid_cells_gdf.crs)

    for name, data in loaded_amenities.items():
        if data['check_intersection'] and data['gdf'] is not None and not data['gdf'].empty:
            intersect = gpd.sjoin(centroids_gdf, data['gdf'], how='inner', predicate='intersects')
            grid_cells_gdf.loc[intersect.index, f'nearest_{name}'] = 0.0

    # Filter water cells
    logging.info("Fetching administrative boundary to filter cells...")
    cardiff_boundary = ox.geocode_to_gdf(CITY).to_crs(TARGET_CRS)
    
    logging.info(f"Grid cells before filtering: {len(grid_cells_gdf)}")
    grid_cells_gdf = gpd.clip(grid_cells_gdf, cardiff_boundary)
    logging.info(f"Grid cells after filtering: {len(grid_cells_gdf)}")

    # Change the crs to lat/long so it can be visualised
    grid_cells_gdf = grid_cells_gdf.to_crs('EPSG:4326')

    # SAVE RESULTS
    logging.info("Saving results...")
    logging.getLogger().setLevel(logging.DEBUG)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_output_dir = os.path.join(OUTPUT_DATA_DIR, f"{CITY.split(',')[0]}_{timestamp}")
    os.makedirs(run_output_dir, exist_ok=True)

    grid_cells_gdf.to_file(os.path.join(run_output_dir, 'grid_cells_accessibility.geojson'), driver='GeoJSON')
    
    # Save feature files dynamically
    for name, data in loaded_amenities.items():
        if data['gdf'] is not None and not data['gdf'].empty:
            clean_gdf = data['gdf'].drop(columns=['centroid'], errors='ignore').to_crs('EPSG:4326')
            clean_gdf.to_file(os.path.join(run_output_dir, data['filename']), driver='GeoJSON')

    print(f"\nSuccess! Results saved to {run_output_dir}")


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
