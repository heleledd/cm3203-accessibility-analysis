import os
from tqdm import tqdm
from datetime import datetime
import logging
import math
import concurrent.futures
import osmnx as ox

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
from shortest_route import get_shortest_route

def process_chunk(args):
    # unpack from the tuple that was passed in (put together in chunk_args)
    chunk_df, G, park_access, park_bounds, gps, schools = args 

    results = []

    # Each chunk gets its own local cache to avoid IPC locks/overhead
    park_cache, gp_cache, school_cache = {}, {}, {}

    # Pre-calculate nearest nodes for all centroids in the chunk to avoid rebuilding the spatial index repeatedly
    centroids = chunk_df.geometry.centroid
    nearest_nodes = ox.distance.nearest_nodes(G, centroids.x.values, centroids.y.values)

    for (idx, row), start_node in zip(chunk_df.iterrows(), nearest_nodes):
        centroid = row.geometry.centroid
        
        p_dist = get_shortest_route(G, centroid, park_access, park_cache, boundaries_gdf=park_bounds, start_node_id=start_node)
        g_dist = get_shortest_route(G, centroid, gps, gp_cache, start_node_id=start_node)
        s_dist = get_shortest_route(G, centroid, schools, school_cache, boundaries_gdf=schools, start_node_id=start_node)
        
        results.append((idx, p_dist, g_dist, s_dist))
        
    return results

def main(
        bbox=CARDIFF_BBOX,
        grid_size=GRID_SIZE_METERS,
        park_access_points_path=PARK_ACCESS_POINTS_PATH,
        park_boundary_path=PARK_BOUNDARY_PATH,
        num_workers=NUM_WORKERS
    ):
    
    # load street network
    G = get_street_network_graph(bbox)
    
    # split the bounding box into a grid
    bbox_reprojected = reproject_bbox(bbox)
    grid_cells_gdf = split_bbox_into_grid(bbox_reprojected, grid_size)

    # load amenity features
    park_access_points_gdf, park_boundaries_gdf = get_park_data(G, park_access_points_path, park_boundary_path)
    gps_gdf = get_osm_features(G, bbox, tags={"amenity": "doctors"})
    schools_gdf = get_osm_features(G, bbox, tags={"amenity": "school"})

    # Setup columns in dataframe
    grid_cells_gdf['nearest_park'] = float('nan')
    grid_cells_gdf['nearest_gp'] = float('nan')
    grid_cells_gdf['nearest_school'] = float('nan')

    # for cleaner tqdm output
    logging.getLogger().setLevel(logging.WARNING)

    logging.info(f"Starting multiprocessing with {num_workers} cores...")

    # split the grid cells dataframe into chunks 
    num_chunks = num_workers * 4
    chunk_size = math.ceil(len(grid_cells_gdf) / num_chunks)
    
    chunks = [
        grid_cells_gdf.iloc[i : i + chunk_size] 
        for i in range(0, len(grid_cells_gdf), chunk_size)
    ]

    # prepare the inputs for the worker
    chunk_args = [
        (chunk, G, park_access_points_gdf, park_boundaries_gdf, gps_gdf, schools_gdf) 
        for chunk in chunks
    ]

    all_results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_chunk, arg) for arg in chunk_args]

        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing Grid Chunks"):
            all_results.extend(future.result())
    
    # Map the gathered results back to the original dataframe
    for idx, p_dist, g_dist, s_dist in all_results:
        grid_cells_gdf.at[idx, 'nearest_park'] = p_dist
        grid_cells_gdf.at[idx, 'nearest_gp'] = g_dist
        grid_cells_gdf.at[idx, 'nearest_school'] = s_dist


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
