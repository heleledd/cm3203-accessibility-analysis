import osmnx as ox
from shapely.geometry import Point, box
import geopandas as gpd
import networkx as nx
import os
from datetime import datetime
from tqdm import tqdm
import logging
import multiprocessing

# Set up paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

if os.environ.get('INPUT_DATA_DIR'):
    INPUT_DATA_DIR = os.environ.get('INPUT_DATA_DIR')
elif os.environ.get('GITHUB_ACTIONS') == 'true':
    INPUT_DATA_DIR = os.path.join(os.getcwd(), 'data', 'input_data')
else:
    INPUT_DATA_DIR = os.path.join(SCRIPT_DIR, '../data/input_data')

if os.environ.get('OUTPUT_DATA_DIR'):
    OUTPUT_DATA_DIR = os.environ.get('OUTPUT_DATA_DIR')
elif os.environ.get('GITHUB_ACTIONS') == 'true':
    OUTPUT_DATA_DIR = os.path.join(os.getcwd(), 'data')
else:
    OUTPUT_DATA_DIR = os.path.join(SCRIPT_DIR, '../output')

os.makedirs(INPUT_DATA_DIR, exist_ok=True)

TARGET_CRS = 'EPSG:27700'

# ==========================================
# DOWNLOAD STREET NETWORK 
# ==========================================

def get_street_network_graph(bbox):
    """Returns a MultiDiGraph of the walkable network around a bounding box, projected to EPSG:27700."""
    network_path = os.path.join(INPUT_DATA_DIR, 'cardiff_network.graphml')
    if os.path.exists(network_path):
        logging.info("Street network graph already exists. Loading from file...")
        G = ox.io.load_graphml(network_path)
        
        # Ensure it's projected to EPSG:27700
        sample_x = list(G.nodes(data=True))[0][1]['x']
        if sample_x < 1000:
            logging.info("Graph appears to be in WGS84, reprojecting to EPSG:27700...")
            G = ox.projection.project_graph(G, to_crs=TARGET_CRS)
        return G
    else:
        logging.info("Street network graph not found. Downloading from OpenStreetMap...")
        G = ox.graph_from_bbox(bbox, network_type="walk", simplify=False, retain_all=False)
        ox.distance.add_edge_lengths(G)
        G_proj = ox.projection.project_graph(G, to_crs=TARGET_CRS)
        ox.save_graphml(G_proj, filepath=network_path)
        return G_proj

def reproject_bbox(bbox):
    """Reproject a WGS84 bounding box to EPSG:27700."""
    minx, miny, maxx, maxy = bbox
    min_point = Point(minx, miny)
    max_point = Point(maxx, maxy)

    min_point_reprojected = gpd.GeoSeries([min_point], crs='EPSG:4326').to_crs(TARGET_CRS).iloc[0]
    max_point_reprojected = gpd.GeoSeries([max_point], crs='EPSG:4326').to_crs(TARGET_CRS).iloc[0]

    return (
        min_point_reprojected.x,
        min_point_reprojected.y,
        max_point_reprojected.x,
        max_point_reprojected.y
    )

# ==========================================
# DIVIDE THE BOUNDING BOX INTO A GRID
# ==========================================

def split_bbox_into_grid(bbox_reprojected, grid_size):
    """Split a projected bounding box into a grid of smaller cells."""
    grid_path = os.path.join(INPUT_DATA_DIR, f'cardiff_grid_cells_{grid_size}m.geojson')
    if os.path.exists(grid_path):
        logging.info("Cardiff grid cells GeoJSON already exists. Loading from file...")
        gdf = gpd.read_file(grid_path)
        return gdf.to_crs(TARGET_CRS)
    else:
        logging.info("Creating Cardiff Grid Cells GeoJSON...")
        minx, miny, maxx, maxy = bbox_reprojected
        grid_cells, grid_id = [], 0

        x = minx
        while x < maxx:
            y = miny
            while y < maxy:
                cell = box(x, y, min(x + grid_size, maxx), min(y + grid_size, maxy))
                grid_cells.append({'grid_id': grid_id, 'geometry': cell})
                grid_id += 1
                y += grid_size
            x += grid_size

        grid_gdf = gpd.GeoDataFrame(grid_cells, crs=TARGET_CRS)
        grid_gdf.to_file(grid_path, driver="GeoJSON")
        return grid_gdf


# ==========================================
# LOAD AMENITY FEATURES
# ==========================================

def find_nearest_boundary_node(G, boundary):
    num_samples = max(10, int(boundary.length / 20))
    sample_points = [boundary.interpolate(i / num_samples, normalized=True) for i in range(num_samples)]
    xs = [p.x for p in sample_points]
    ys = [p.y for p in sample_points]
    return ox.distance.nearest_nodes(G, xs, ys)

def get_park_boundary_nodes(G, park_path):
    """Loads parks from local GeoJSON and maps their boundaries to network nodes."""
    park_gdf = gpd.read_file(park_path).to_crs(TARGET_CRS)
    park_gdf['boundary'] = park_gdf.geometry.boundary
    logging.info("Mapping park boundary nodes (this may take a moment)...")
    park_gdf['nearest_node'] = park_gdf['boundary'].apply(lambda b: find_nearest_boundary_node(G, b))
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


# ==========================================
# FIND NEAREST AMENITY
# ==========================================

def get_nearest_node_from_array(G, start_node_id, candidate_nodes):
    start_x, start_y = G.nodes[start_node_id]['x'], G.nodes[start_node_id]['y']
    best_node, min_dist = None, float('inf')
    for node_id in candidate_nodes:
        dist = ((start_x - G.nodes[node_id]['x']) ** 2 + (start_y - G.nodes[node_id]['y']) ** 2) ** 0.5
        if dist < min_dist:
            min_dist = dist
            best_node = node_id
    return best_node


def get_shortest_park_route(G, start_node, features_gdf, cache):
    """Boundary-based routing for Parks."""
    if features_gdf.contains(start_node).any():
        return 0.0

    nearest_node_id = ox.distance.nearest_nodes(G, start_node.x, start_node.y)
    
    if nearest_node_id in cache:
        return cache[nearest_node_id]

    start_point = Point(G.nodes[nearest_node_id]['x'], G.nodes[nearest_node_id]['y'])
    features_gdf['_straight_line_dist'] = features_gdf['boundary'].distance(start_point)
    nearest_parks = features_gdf.nsmallest(2, '_straight_line_dist')

    distances = []
    for _, row in nearest_parks.iterrows():
        dest_node = get_nearest_node_from_array(G, nearest_node_id, row['nearest_node'])
        try:
            route = ox.routing.shortest_path(G, nearest_node_id, dest_node, weight='length', cpus=1)
            if route:
                distances.append(round(nx.path_weight(G, route, weight='length'), 2))
        except nx.NetworkXNoPath:
            pass

    closest_dist = min(distances) if distances else None
    cache[nearest_node_id] = closest_dist
    return closest_dist


def get_shortest_centroid_route(G, start_node, features_gdf, cache):
    """Centroid-based routing for GPs and Schools."""
    nearest_node_id = ox.distance.nearest_nodes(G, start_node.x, start_node.y)

    if nearest_node_id in cache:
        return cache[nearest_node_id]

    features_gdf['distance_m'] = features_gdf.geometry.distance(start_node)
    features_ranked = features_gdf.sort_values('distance_m').head(3)

    distances = []
    for _, row in features_ranked.iterrows():
        try:
            route = ox.routing.shortest_path(G, nearest_node_id, row['nearest_node'], weight='length', cpus=1)
            if route:
                distances.append(round(nx.path_weight(G, route, weight='length'), 2))
        except nx.NetworkXNoPath:
            pass

    closest_dist = min(distances) if distances else None
    cache[nearest_node_id] = closest_dist
    return closest_dist


# ==========================================
# MAIN
# ==========================================

def main(bbox, grid_size, parks_path, github_actions):
    # Initialize Network and Grid
    G = get_street_network_graph(bbox)
    bbox_reprojected = reproject_bbox(bbox)
    grid_cells_gdf = split_bbox_into_grid(bbox_reprojected, grid_size)

    # Prepare destination features
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

    # Export results
    if github_actions:
        run_output_dir = OUTPUT_DATA_DIR
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_output_dir = os.path.join(OUTPUT_DATA_DIR, timestamp)
        os.makedirs(run_output_dir, exist_ok=True)

    grid_cells_gdf.to_file(os.path.join(run_output_dir, 'grid_cells_accessibility.geojson'), driver='GeoJSON')
    
    # Save features used with the distance file in the same folder
    parks_gdf.drop(columns=['boundary'], errors='ignore').to_file(os.path.join(run_output_dir, 'parks.geojson'), driver='GeoJSON')
    gps_gdf.drop(columns=['centroid'], errors='ignore').to_file(os.path.join(run_output_dir, 'gps.geojson'), driver='GeoJSON')
    schools_gdf.drop(columns=['centroid'], errors='ignore').to_file(os.path.join(run_output_dir, 'schools.geojson'), driver='GeoJSON')

    print(f"\nSuccess! Results saved to {run_output_dir}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    
    # Cardiff city centre bbox for testing (west, south, east, north)
    # cardiff_bbox = (-3.21, 51.49, -3.16, 51.51)

    cardiff_bbox = (-3.35, 51.37, -3.05, 51.57)
    grid_size_meters = 100

    # Define parks path
    parks_file = 'cardiff_parks.geojson'
    parks_path = os.path.join(INPUT_DATA_DIR, parks_file)

    # let the user know where the script thinks it's running
    is_github_actions = os.environ.get('GITHUB_ACTIONS') == 'true'
    if is_github_actions:
        print("Running in GitHub Actions!")
    else:
        print("Running locally!")

    main(cardiff_bbox, grid_size_meters, parks_path, is_github_actions)