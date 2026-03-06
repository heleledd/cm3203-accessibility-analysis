import networkx as nx
import osmnx as ox
from shapely.geometry import Point


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