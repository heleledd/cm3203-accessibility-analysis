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

""" returns a MultiDigraph of the walkable network around a location point"""
def get_street_network_graph(bbox=(-3.35, 51.37, -3.05, 51.57)):
    network_path = os.path.join(DATA_DIR, 'cardiff_network.graphml')
    if os.path.exists(network_path):
        logging.info("Street network graph already exists. Loading from file...")
        G = ox.io.load_graphml(network_path)
        logging.info("Street network graph loaded successfully.")
        return G
    else:
        logging.info("Street network graph not found. Downloading from OpenStreetMap...")

        G =  ox.graph_from_bbox(bbox, network_type="walk", simplify=False, retain_all=False)
        
        # add edge lengths to the graph (instructions say to do this before projecting or simplifying the graph)
        ox.distance.add_edge_lengths(G)
        
        G_proj = ox.projection.project_graph(G, to_crs="EPSG:27700")

        ox.save_graphml(G_proj, filepath=network_path)
        logging.info("Street network graph downloaded and saved to file.")

        return G_proj


"""Split a bounding box into a grid of smaller cells in a GeoDataFrame"""
def split_bbox_into_grid(bbox, grid_size):
    grid_path = os.path.join(DATA_DIR, f'cardiff_grid_cells_{grid_size}m.geojson')
    if os.path.exists(grid_path):
        logging.info("Cardiff grid cells GeoJSON file already exists. Loading from file...")
        gdf = gpd.read_file(grid_path)
        logging.info("Cardiff grid cells GeoJSON file loaded successfully.")
        return gdf
    else:
        logging.info("Cardiff Grid Cells GeoJSON file not found. Processing ONS postcodes CSV file...")
        minx, miny, maxx, maxy = bbox

        grid_cells = []
        grid_id = 0

        # create cells from the bottom-left to the top-right of the bounding box
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

        # Create GeoDataFrame
        grid_gdf = gpd.GeoDataFrame(grid_cells, crs='EPSG:27700')
        
        # # Add centroid for routing calculations
        # grid_gdf['centroid'] = grid_gdf.geometry.centroid

        # save the geodataframe to json file
        grid_gdf.to_file(grid_path, driver="GeoJSON")
        logging.info("Cardiff grid cells GeoJSON file created and saved to file.")
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
    # Extract node coordinates from the graph
    nodes_data = []
    for node, data in G.nodes(data=True):
        nodes_data.append({
            'node_id': node,
            'geometry': Point(data['x'], data['y'])
        })
    
    # Create GeoDataFrame of nodes
    nodes_gdf = gpd.GeoDataFrame(nodes_data, crs='EPSG:27700')
    
    # Spatial join: find which grid cell each node belongs to
    nodes_with_grid = gpd.sjoin(nodes_gdf, grid_gdf, how='left', predicate='within')
    
    # Group nodes by grid_id
    grid_gdf['nodes'] = None
    for grid_id in grid_gdf['grid_id']:
        nodes_in_cell = nodes_with_grid[nodes_with_grid['grid_id'] == grid_id]['node_id'].tolist()
        grid_gdf.at[grid_gdf[grid_gdf['grid_id'] == grid_id].index[0], 'nodes'] = nodes_in_cell
    
    # Add count of nodes per cell
    grid_gdf['node_count'] = grid_gdf['nodes'].apply(lambda x: len(x) if x else 0)
    
    return grid_gdf


""" download features from OpenStreetMap and find the nearest network node to each feature"""
def get_features(G, bbox, tags={"shop": "hairdresser"}):
    # find the nearest features to the start node and get its node id
    features_gdf = ox.features_from_bbox(bbox, tags=tags).to_crs(epsg=27700)

    # get the centroid of the features and find the nearest node to that centroid
    features_gdf['centroid'] = features_gdf.geometry.centroid

    # Find nearest network nodes for all features (vectorized)
    features_gdf['nearest_node'] = ox.distance.nearest_nodes(
        G,
        features_gdf['centroid'].x,
        features_gdf['centroid'].y
    )

    return features_gdf


""" ranks closest features and finds the shortest route between two nodes in the graph and plots it"""
def get_shortest_route(G, start_node, features_gdf, cache=None):
    """Find shortest network distance from start_node to features.

    Uses a simple cache keyed by the nearest graph node id so repeated
    queries from the same graph node reuse previous shortest-path results.
    """
    if cache is None:
        cache = {}

    # find closest node to the start point
    nearest_node_id, distance_to_node = ox.distance.nearest_nodes(G, start_node.x, start_node.y, return_dist=True)

    # If we've already computed shortest-paths from this graph node, reuse them
    if nearest_node_id in cache:
        logging.debug(f"Reusing cached shortest-path for start node {nearest_node_id}")
        return cache[nearest_node_id]['closest']

    # find each feature's distance to the start node
    features_gdf['distance_m'] = features_gdf.geometry.distance(start_node)
    
    # Sort by distance (closest first)
    features_ranked = features_gdf.sort_values('distance_m')

    # print("Closest features ranked by distance:")
    # print(features_ranked[['name', 'distance_m']] if 'name' in features_ranked.columns else features_ranked[['distance_m']])

    # Calculate distance through the network to the top 3 geographically closest results
    distances = {}
    for idx, row in features_ranked.iloc[:3].iterrows():
        try:
            
            # Calculate shortest route to the destination based on length 
            route = ox.routing.shortest_path(
                G, 
                nearest_node_id, 
                row['nearest_node'], 
                weight='length',
                cpus=1
            )

            # calculate the actual distance of the route
            route_length = round(nx.path_weight(G, route, weight='length'), 2)

            tqdm.write(f"Route length from start node {nearest_node_id} to feature node {row['nearest_node']}: {route_length}")

            distances[row['nearest_node']] = route_length
        except nx.NetworkXNoPath:
            logging.warning(f"No path found to feature: {row['name'] if 'name' in row else idx}")
            distances[row['nearest_node']] = float('inf')  # Assign infinity if no path is found
            continue

    
    # discover which one was the shortest route
    if distances:
        closest_feature_node = min(distances, key=distances.get)
        logging.info(f"Closest feature node based on route length: {closest_feature_node} with distance {distances[closest_feature_node]} metres")
        # cache the computed distances (store both mapping and chosen minimum)
        cache[nearest_node_id] = {
            'distances': distances,
            'closest': distances[closest_feature_node]
        }

        # return shortest distance to a feature
        return distances[closest_feature_node]
    else:
        logging.warning("No features found within the specified distance.")
        return None



def main(bbox):
    # get the street network graph from within the bounding box
    G = get_street_network_graph(bbox)

    # change the bbox to the same projection as the graph (EPSG:27700)
    minx, miny, maxx, maxy = bbox
    min_point = Point(minx, miny)
    max_point = Point(maxx, maxy)

    # Reproject the points to the graph's CRS (EPSG:27700)
    min_point_reprojected = gpd.GeoSeries([min_point], crs='EPSG:4326').to_crs('EPSG:27700').iloc[0]
    max_point_reprojected = gpd.GeoSeries([max_point], crs='EPSG:4326').to_crs('EPSG:27700').iloc[0]

    # Update bbox with reprojected coordinates
    bbox_reprojected = (
        min_point_reprojected.x,
        min_point_reprojected.y,
        max_point_reprojected.x,
        max_point_reprojected.y
    )

    # split the bounding box into a grid
    grid_cells_gdf = split_bbox_into_grid(bbox_reprojected, grid_size=100) # 100m grid cells

    # grid_cells_gdf = assign_nodes_to_grid_cells(G, grid_cells_gdf)

    # print(f"Total grid cells: {len(grid_cells_gdf)}")
    # print(f"Cells with nodes: {(grid_cells_gdf['node_count'] > 0).sum()}")
    # print(f"Cells without nodes: {(grid_cells_gdf['node_count'] == 0).sum()}")
    # print(f"Average nodes per cell: {grid_cells_gdf['node_count'].mean():.2f}")

    # prepare a column to store the shortest network distance (metres) to a hospital
    grid_cells_gdf['nearest_hospital'] = float('nan')

    # get the hospital features
    hospital_tags = {"amenity": "hospital"}
    features_gdf = get_features(G, bbox, tags=hospital_tags)

    
    # iterate through the filtered grid cells and find the shortest distance to a feature for each one
    # use a simple cache so repeated queries from the same nearest graph node are reused
    route_cache = {}
    
    # Suppress verbose logging during processing for cleaner output
    logging.getLogger().setLevel(logging.WARNING)
    
    for idx, row in tqdm(grid_cells_gdf.iterrows(), total=len(grid_cells_gdf), desc="Processing grid cells"):
        cell_centroid = row.geometry.centroid
        shortest_distance = get_shortest_route(G, cell_centroid, features_gdf, cache=route_cache)
        # maybe add the distance to the node on top as well??
        
        # store the returned shortest distance (None if not found)
        grid_cells_gdf.at[idx, 'nearest_hospital'] = shortest_distance
    
    # Restore logging level
    logging.getLogger().setLevel(logging.DEBUG)


    # today's date and time for the output file name
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # save the results to a new GeoJSON file
    output_path = os.path.join(DATA_DIR, f'grid_cells_with_hospital_distances_{timestamp}.geojson')
    grid_cells_gdf.to_file(output_path, driver='GeoJSON')


def __main__():
    # set up logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(levelname)s:%(asctime)s - %(message)s'
    )
    
    # (east longitude, south latitude, west longitude, north latitude)
    
    # all of cardiff
    cardiff_bbox_lat_long = (-3.35, 51.37, -3.05, 51.57)


    # cardiff city centre
    smaller_bbox_lat_long = (-3.21, 51.49, -3.16, 51.51)

    main(cardiff_bbox_lat_long)


if __name__ == "__main__":
    __main__()