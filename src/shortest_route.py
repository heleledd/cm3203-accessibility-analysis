import networkx as nx
import osmnx as ox

def get_shortest_route(G, start_point, target_gdf, cache, boundaries_gdf=None, k_nearest=3):
    """
    Unified routing function for Parks, GPs, and Schools.
    
    Parameters:
    - boundaries_gdf: Optional. Pass park boundaries here to trigger the 0.0m 'inside' check.
    - k_nearest: Number of nearest destinations to attempt routing to.
    """
    
    # if the node is actually inside a park then return distance as 0m
    if boundaries_gdf is not None:
        # sindex.query returns a 2D array of indices. If the second array has items, it means the point intersects with one or more polygons.
        intersecting_indices = boundaries_gdf.sindex.query(start_point, predicate="intersects")
        if len(intersecting_indices[1]) > 0:
            return 0.0

    # get nearest network node and check cache
    nearest_node_id = ox.distance.nearest_nodes(G, start_point.x, start_point.y)
    
    if nearest_node_id in cache:
        return cache[nearest_node_id]

    # find the nearest target destinations using straight-line distance (target_gdf contains 'Point' geometry type)
    target_gdf['_temp_dist'] = target_gdf.geometry.distance(start_point)
    nearest_targets = target_gdf.sort_values('_temp_dist').head(k_nearest)

    # 4. Attempt Routing
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
