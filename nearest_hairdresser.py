import osmnx as ox

location_point = (51.496103, -3.173560) # 115 mackintosh place
distance = 750  # meters

def get_graph(location_point, distance):
    G =  ox.graph_from_point(location_point, dist=distance, dist_type="network", network_type="walk")
    G_proj = ox.projection.project_graph(G)

    G2 = ox.simplification.consolidate_intersections(
        G_proj,
        rebuild_graph=True,
        tolerance=15,
        dead_ends=True,
    )
    return G2


def __main__():
    print("Start of the code...")