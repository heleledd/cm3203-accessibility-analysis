import networkx as nx
import osmnx as ox

def get_shortest_route(G, start_point, target_gdf, cache, boundaries_gdf=None, start_node_id=None, k_nearest=3):
    """
    Unified routing function for Parks, GPs, and Schools.
    
    Parameters:
    - boundaries_gdf: Optional. Pass park boundaries here to trigger the 0.0m 'inside' check.
    - k_nearest: Number of nearest destinations to attempt routing to.
    """
    
    # if the node is actually inside a park then return distance as 0m
    if boundaries_gdf is not None:
        # sindex.query returns an array of indices. If it has items, it means the point intersects with one or more polygons.
        intersecting_indices = boundaries_gdf.sindex.query(start_point, predicate="intersects")
        if len(intersecting_indices) > 0:
            return 0.0

    # get nearest network node and check cache
    if start_node_id is None:
        nearest_node_id = ox.distance.nearest_nodes(G, start_point.x, start_point.y)
    else:
        nearest_node_id = start_node_id
    
    if nearest_node_id in cache:
        return cache[nearest_node_id]

    # find the nearest target destinations using straight-line distance using an R-Tree
    nearest_indices = target_gdf.sindex.nearest(start_point, return_all=False)[1]

    # Get only the nearest rows, up to k_nearest
    nearest_targets = target_gdf.iloc[nearest_indices].head(k_nearest)

    # attempt routing
    distances = []
    for _, row in nearest_targets.iterrows():
        dest_node = row['nearest_node'] 
        
        try:
            route = ox.routing.shortest_path(G, nearest_node_id, dest_node, weight='length', cpus=1)
            if route:
                distances.append(round(nx.path_weight(G, route, weight='length'), 2))
        except nx.NetworkXNoPath:
            pass

    # determine closest distance from the nearest targets, cache result from start node, and return
    closest_dist = min(distances) if distances else None
    cache[nearest_node_id] = closest_dist
    return closest_dist
