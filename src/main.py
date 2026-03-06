import os
import geopandas as gpd
from tqdm import tqdm
from datetime import datetime
import multiprocessing
import logging

from config import (
    CARDIFF_BBOX,
    GRID_SIZE_METERS,
    INPUT_DATA_DIR,
    OUTPUT_DATA_DIR,
)

from network_and_bbox import get_street_network_graph, reproject_bbox, split_bbox_into_grid
from amenities import get_park_boundary_nodes, get_osm_features
from shortest_route import get_shortest_park_route, get_shortest_centroid_route

def main(
        bbox=CARDIFF_BBOX,
        grid_size=GRID_SIZE_METERS,
        parks_path=os.path.join(INPUT_DATA_DIR, 'cardiff_parks.geojson'),
        num_workers=None
    ):
    
    # load street network
    G = get_street_network_graph(bbox)
    
    # split the bounding box into a grid
    bbox_reprojected = reproject_bbox(bbox)
    grid_cells_gdf = split_bbox_into_grid(bbox_reprojected, grid_size)

    # load amenity features
    parks_gdf = get_park_boundary_nodes(G, parks_path)
    gps_gdf = get_osm_features(G, bbox, tags={"amenity": "doctors"})
    schools_gdf = get_osm_features(G, bbox, tags={"amenity": "school"})

    # Setup columns in dataframe
    grid_cells_gdf['nearest_park'] = float('nan')
    grid_cells_gdf['nearest_gp'] = float('nan')
    grid_cells_gdf['nearest_school'] = float('nan')

    park_cache, gp_cache, school_cache = {}, {}, {}

    # for cleaner tqdm output
    logging.getLogger().setLevel(logging.WARNING)

    # Main processing loop
    for idx, row in tqdm(grid_cells_gdf.iterrows(), total=len(grid_cells_gdf), desc="Routing from grid cells"):
        centroid = row.geometry.centroid

        grid_cells_gdf.at[idx, 'nearest_park'] = get_shortest_park_route(G, centroid, parks_gdf, park_cache)
        grid_cells_gdf.at[idx, 'nearest_gp'] = get_shortest_centroid_route(G, centroid, gps_gdf, gp_cache)
        grid_cells_gdf.at[idx, 'nearest_school'] = get_shortest_centroid_route(G, centroid, schools_gdf, school_cache)

    logging.getLogger().setLevel(logging.DEBUG)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_output_dir = os.path.join(OUTPUT_DATA_DIR, timestamp)
    os.makedirs(run_output_dir, exist_ok=True)

    grid_cells_gdf.to_file(os.path.join(run_output_dir, 'grid_cells_accessibility.geojson'), driver='GeoJSON')
    
    # Save features used with the distance file in the same folder
    parks_gdf.drop(columns=['boundary'], errors='ignore').to_file(os.path.join(run_output_dir, 'parks.geojson'), driver='GeoJSON')
    gps_gdf.drop(columns=['centroid'], errors='ignore').to_file(os.path.join(run_output_dir, 'gps.geojson'), driver='GeoJSON')
    schools_gdf.drop(columns=['centroid'], errors='ignore').to_file(os.path.join(run_output_dir, 'schools.geojson'), driver='GeoJSON')

    print(f"\nSuccess! Results saved to {run_output_dir}")



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()