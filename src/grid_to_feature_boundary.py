# imported park data from 'dissolved_parks2.geojson'

import osmnx as ox
import matplotlib.pyplot as plt
from shapely.geometry import Point, box
import geopandas as gpd
import networkx as nx
import os
from datetime import datetime
from tqdm import tqdm
import logging

# Set up paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '../data')

TARGET_CRS = 'EPSG:27700'


""" returns a MultiDigraph of the walkable network around a location point"""
def get_street_network_graph(bbox):
    network_path = os.path.join(DATA_DIR, 'cardiff_network.graphml')
    if os.path.exists(network_path):
        logging.info("Street network graph already exists. Loading from file...")
        G = ox.io.load_graphml(network_path)

        # GraphML may not preserve projection — reproject to EPSG:27700 if needed
        sample_x = list(G.nodes(data=True))[0][1]['x']
        if sample_x < 1000:
            logging.info("Graph appears to be in WGS84, reprojecting to EPSG:27700...")
            G = ox.projection.project_graph(G, to_crs=TARGET_CRS)

        logging.info(f"Street network graph loaded. Sample node x: {list(G.nodes(data=True))[0][1]['x']:.2f}")
        return G
    else:
        logging.info("Street network graph not found. Downloading from OpenStreetMap...")

        G = ox.graph_from_bbox(bbox, network_type="walk", simplify=False, retain_all=False)

        # add edge lengths before projecting
        ox.distance.add_edge_lengths(G)

        G_proj = ox.projection.project_graph(G, to_crs=TARGET_CRS)

        ox.save_graphml(G_proj, filepath=network_path)
        logging.info("Street network graph downloaded and saved to file.")

        return G_proj


"""Split a bounding box into a grid of smaller cells in a GeoDataFrame"""
def split_bbox_into_grid(bbox, grid_size):
    # bbox must already be in EPSG:27700
    grid_path = os.path.join(DATA_DIR, f'cardiff_grid_cells_{grid_size}m.geojson')
    if os.path.exists(grid_path):
        logging.info("Cardiff grid cells GeoJSON file already exists. Loading from file...")
        gdf = gpd.read_file(grid_path)
        # GeoJSON is always stored as WGS84 — reproject back to EPSG:27700
        gdf = gdf.to_crs(TARGET_CRS)
        logging.info(f"Cardiff grid cells loaded and reprojected. Sample centroid: {gdf.geometry.iloc[0].centroid.x:.2f}")
        return gdf
    else:
        logging.info("Cardiff Grid Cells GeoJSON file not found. Creating grid...")
        minx, miny, maxx, maxy = bbox

        grid_cells = []
        grid_id = 0

        x = minx
        while x < maxx:
            y = miny
            while y < maxy:
                cell = box(x, y, min(x + grid_size, maxx), min(y + grid_size, maxy))
                grid_cells.append({
                    'grid_id': grid_id,
                    'geometry': cell
                })
                grid_id += 1
                y += grid_size
            x += grid_size

        grid_gdf = gpd.GeoDataFrame(grid_cells, crs=TARGET_CRS)

        grid_gdf.to_file(grid_path, driver="GeoJSON")
        logging.info("Cardiff grid cells GeoJSON file created and saved.")
        return grid_gdf


def assign_nodes_to_grid_cells(G, grid_gdf):
    """
    Find which nodes from the street network graph fall within each grid cell.

    Args:
        G: NetworkX graph from OSMnx (projected to EPSG:27700)
        grid_gdf: GeoDataFrame of grid cells (in EPSG:27700)

    Returns:
        grid_gdf with a new 'nodes' column containing list of node IDs in each cell
    """
    nodes_data = []
    for node, data in G.nodes(data=True):
        nodes_data.append({
            'node_id': node,
            'geometry': Point(data['x'], data['y'])
        })

    nodes_gdf = gpd.GeoDataFrame(nodes_data, crs=TARGET_CRS)

    nodes_with_grid = gpd.sjoin(nodes_gdf, grid_gdf, how='left', predicate='within')

    grid_gdf['nodes'] = None
    for grid_id in grid_gdf['grid_id']:
        nodes_in_cell = nodes_with_grid[nodes_with_grid['grid_id'] == grid_id]['node_id'].tolist()
        grid_gdf.at[grid_gdf[grid_gdf['grid_id'] == grid_id].index[0], 'nodes'] = nodes_in_cell

    grid_gdf['node_count'] = grid_gdf['nodes'].apply(lambda x: len(x) if x else 0)

    return grid_gdf


def find_nearest_boundary_node(G, boundary):
    """Find the network node nearest to any point on the park boundary.
    Both the graph and boundary must be in EPSG:27700.
    """
    num_samples = max(10, int(boundary.length / 20))  # ~1 sample per 20m

    sample_points = [
        boundary.interpolate(i / num_samples, normalized=True)
        for i in range(num_samples)
    ]

    xs = [p.x for p in sample_points]
    ys = [p.y for p in sample_points]

    candidate_nodes = ox.distance.nearest_nodes(G, xs, ys)

    return candidate_nodes


""" load OS maps park data and find the nodes closest to the boundary """
def get_park_boundary_nodes(G):
    park_path = os.path.join(DATA_DIR, 'dissolved_parks2.geojson')
    park_gdf = gpd.read_file(park_path)

    # GeoJSON is always WGS84 — reproject to match graph
    park_gdf = park_gdf.to_crs(TARGET_CRS)
    logging.info(f"Parks loaded and reprojected. Sample boundary point x: {park_gdf.geometry.iloc[0].centroid.x:.2f}")

    park_gdf['boundary'] = park_gdf.geometry.boundary

    park_gdf['nearest_node'] = park_gdf['boundary'].apply(
        lambda boundary: find_nearest_boundary_node(G, boundary)
    )


    return park_gdf


def get_nearest_node_from_array(G, start_node_id, candidate_nodes):
    """Return the node from candidate_nodes that is geographically closest to start_node_id."""
    start_x = G.nodes[start_node_id]['x']
    start_y = G.nodes[start_node_id]['y']
    
    best_node = None
    min_dist = float('inf')
    
    for node_id in candidate_nodes:
        nx_ = G.nodes[node_id]['x']
        ny_ = G.nodes[node_id]['y']
        dist = ((start_x - nx_) ** 2 + (start_y - ny_) ** 2) ** 0.5
        if dist < min_dist:
            min_dist = dist
            best_node = node_id
    
    return best_node



def get_shortest_park_route(G, start_node, features_gdf, cache=None):
    """start_node must be a shapely Point in EPSG:27700."""
    if cache is None:
        cache = {}

    nearest_node_id, _ = ox.distance.nearest_nodes(
        G, start_node.x, start_node.y, return_dist=True
    )

    start_x = G.nodes[nearest_node_id]['x']
    start_y = G.nodes[nearest_node_id]['y']
    start_point = Point(start_x, start_y)

    if nearest_node_id in cache:
        return cache[nearest_node_id]['closest']

    features_gdf['_straight_line_dist'] = features_gdf['boundary'].distance(start_point)
    nearest_parks = features_gdf.nsmallest(2, '_straight_line_dist')

    distances = {}

    for idx, row in nearest_parks.iterrows():
        # find the geographically closest node in the array
        candidate_nodes = row['nearest_node']
        dest_node = get_nearest_node_from_array(G, nearest_node_id, candidate_nodes)
        try:
            route = ox.routing.shortest_path(
                G, nearest_node_id, dest_node, weight='length', cpus=1
            )
            if route is None:
                continue
            route_length = round(nx.path_weight(G, route, weight='length'), 2)
            tqdm.write(f"Route length from {nearest_node_id} to boundary node {dest_node}: {route_length}")
            distances[dest_node] = route_length
        except nx.NetworkXNoPath:
            logging.warning(f"No path found to feature: {idx}")
            distances[dest_node] = float('inf')

    if distances:
        closest_dist = min(distances.values())
        cache[nearest_node_id] = {'closest': closest_dist}
        return closest_dist
    else:
        logging.warning("No features found.")
        return None


def reproject_bbox(bbox):
    """Reproject a WGS84 bounding box to EPSG:27700."""
    minx, miny, maxx, maxy = bbox
    min_point = Point(minx, miny)
    max_point = Point(maxx, maxy)

    min_point_reprojected = gpd.GeoSeries([min_point], crs='EPSG:4326').to_crs(TARGET_CRS).iloc[0]
    max_point_reprojected = gpd.GeoSeries([max_point], crs='EPSG:4326').to_crs(TARGET_CRS).iloc[0]

    bbox_reprojected = (
        min_point_reprojected.x,
        min_point_reprojected.y,
        max_point_reprojected.x,
        max_point_reprojected.y
    )

    return bbox_reprojected


def main(bbox, grid_size, tags, feature_name):
    G = get_street_network_graph(bbox)

    # Confirm graph is in EPSG:27700
    sample_node = list(G.nodes(data=True))[0]
    logging.info(f"Graph sample node coords: {sample_node[1]['x']:.2f}, {sample_node[1]['y']:.2f}")

    bbox_reprojected = reproject_bbox(bbox)
    logging.info(f"Reprojected bbox: {bbox_reprojected}")

    grid_cells_gdf = split_bbox_into_grid(bbox_reprojected, grid_size)

    # Confirm grid is in EPSG:27700
    sample_centroid = grid_cells_gdf.geometry.iloc[0].centroid
    logging.info(f"Grid sample centroid: {sample_centroid.x:.2f}, {sample_centroid.y:.2f}")

    grid_cells_gdf[feature_name] = float('nan')

    features_gdf = get_park_boundary_nodes(G)

    route_cache = {}

    logging.getLogger().setLevel(logging.WARNING)

    for idx, row in tqdm(grid_cells_gdf.iterrows(), total=len(grid_cells_gdf), desc="Processing grid cells"):
        cell_centroid = row.geometry.centroid
        shortest_distance = get_shortest_park_route(G, cell_centroid, features_gdf, cache=route_cache)
        grid_cells_gdf.at[idx, 'nearest_park'] = shortest_distance

    logging.getLogger().setLevel(logging.DEBUG)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = os.path.join(DATA_DIR, f'grid_cells_with_park_distances_{timestamp}.geojson')
    grid_cells_gdf.to_file(output_path, driver='GeoJSON')
    logging.info(f"Results saved to {output_path}")


def __main__():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(levelname)s:%(asctime)s - %(message)s'
    )

    # (west longitude, south latitude, east longitude, north latitude)
    cardiff_bbox_lat_long = (-3.35, 51.37, -3.05, 51.57)

    feature_tags = {"landuse": ["recreation_ground", "grass"]}
    feature_name = 'nearest_park'

    grid_size = 100

    main(cardiff_bbox_lat_long, grid_size, feature_tags, feature_name)


if __name__ == "__main__":
    __main__()