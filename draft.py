import networkx as n
import osmnx as ox
import folium
import random
from shapely.geometry import Point, LineString

def get_network_graph(location_point, distance, building_tags=None):
    G =  ox.graph_from_point(location_point, dist=distance, dist_type="network", network_type="walk")
    G_proj = ox.projection.project_graph(G, to_crs="EPSG:27700")

    G2 = ox.simplification.consolidate_intersections(
        G_proj,
        rebuild_graph=True,
        tolerance=15,
        dead_ends=True,
    )
    ox.distance.add_edge_lengths(G2)
    return G2


def add_residences_as_nodes_to_graph(G, residences):
    # add residences as nodes and connect them to the nearest street edge
    G_copy = G.copy()

    # get nodes and edges as GeoDataFrames to use for finding the nearest edge
    nodes, edges = ox.graph_to_gdfs(G_copy)

    residence_node_ids = []

    # find the nearest edge to the residence
    for idx, residence in residences.iterrows():
        centroid = residence.geometry.centroid

        new_node_id = f"residence_{idx}"

        edges['distance_to_residence'] = edges.geometry.distance(centroid)
        nearest_edge = edges.loc[edges['distance_to_residence'].idxmin()]

        # find the closest point on the edge to the residence
        nearest_point_on_edge = nearest_edge.geometry.interpolate(
            nearest_edge.geometry.project(centroid)
        )

        # Calculate actual distance from residence to street
        distance_to_street = centroid.distance(nearest_point_on_edge)
        
        # Add the residence as a new node
        G_copy.add_node(
            new_node_id,
            y=centroid.y,
            x=centroid.x,
            geometry=centroid,
            node_type="residence"
        )
        
        # Get the start and end nodes of the nearest edge
        u, v, key = nearest_edge
        
        # Create a connection node on the street edge
        street_connection_id = f"street_connection_{idx}"
        G_copy.add_node(
            street_connection_id,
            y=nearest_point_on_edge.y,
            x=nearest_point_on_edge.x,
            geometry=nearest_point_on_edge,
            node_type="street_connection"
        )
        
        # Connect residence to the street connection point
        G_copy.add_edge(
            new_node_id,
            street_connection_id,
            length=distance_to_street,
            geometry=LineString([centroid, nearest_point_on_edge]),
            highway="residence_access"
        )
        
        G_copy.add_edge(
            street_connection_id,
            new_node_id,
            length=distance_to_street,
            geometry=LineString([nearest_point_on_edge, centroid]),
            highway="residence_access"
        )
        
        # Now split the original edge into two parts at the connection point
        # Calculate distances along the original edge
        total_edge_length = nearest_edge['length']
        edge_geom = nearest_edge.geometry
        
        # Distance from u to connection point
        dist_u_to_connection = edge_geom.project(nearest_point_on_edge)
        
        # Distance from connection point to v
        dist_connection_to_v = total_edge_length - dist_u_to_connection
        
        # Remove the original edge
        if G_copy.has_edge(u, v, key):
            G_copy.remove_edge(u, v, key)
        
        # Add two new edges: u -> connection and connection -> v
        G_copy.add_edge(
            u,
            street_connection_id,
            length=dist_u_to_connection,
            geometry=LineString(list(edge_geom.coords[:int(dist_u_to_connection)])) if edge_geom else None,
            highway=nearest_edge.get('highway', 'unclassified')
        )
        
        G_copy.add_edge(
            street_connection_id,
            v,
            length=dist_connection_to_v,
            geometry=LineString(list(edge_geom.coords[int(dist_u_to_connection):])) if edge_geom else None,
            highway=nearest_edge.get('highway', 'unclassified')
        )
        
        # If the original edge was bidirectional, add reverse edges too
        if G_copy.has_edge(v, u):
            G_copy.remove_edge(v, u)
            
            G_copy.add_edge(
                street_connection_id,
                u,
                length=dist_u_to_connection,
                highway=nearest_edge.get('highway', 'unclassified')
            )
            
            G_copy.add_edge(
                v,
                street_connection_id,
                length=dist_connection_to_v,
                highway=nearest_edge.get('highway', 'unclassified')
            )
        G_copy.to_file("geopackages/graph_with_residences.gpkg", driver="GPKG")
        residence_node_ids.append(new_node_id)

    
    return G_copy, residence_node_ids



def get_nearest_hairdresser(G, residential_buildings, hairdressers):
    start_point = (51.496103, -3.173560)
    end_point = (51.496497, -3.167885)

    # Find the nearest nodes in your graph
    orig = ox.nearest_nodes(G, start_point[1], start_point[0])  # note: lon, lat order
    dest = ox.nearest_nodes(G, end_point[1], end_point[0])

    # Now calculate the route
    route = ox.routing.shortest_path(G, orig, dest, weight="length")


def display_with_folium(location_point, G):
    # apparently folium only works with lat/lon

    m = folium.Map(location=location_point, zoom_start=15)

    # Convert graph back to EPSG:4326 for Folium
    G_geo = ox.projection.project_graph(G, to_crs="EPSG:4326")
    nodes, edges = ox.graph_to_gdfs(G_geo)
    
    # Plot street network edges
    for _, edge in edges.iterrows():
        coords = list(edge.geometry.coords)
        folium.PolyLine(
            locations=[(lat, lon) for lon, lat in coords],
            color="gray",
            weight=2,
            opacity=0.7
        ).add_to(m)
    
    return m


def __main__():
    location_point = (51.496103, -3.173560) # 115 mackintosh place
    distance = 750  # meters

    G = get_network_graph(location_point, distance)

    # get residential buildings
    residences_tags = {"building": True, "residential": True}
    residences = ox.features_from_point(location_point, tags=residences_tags, dist=distance).to_crs(epsg=27700) 

    # get hairdressers
    shop_tags = {"shop": "hairdresser"}
    hairdressers = ox.features_from_point(location_point, tags=shop_tags, dist=distance).to_crs(epsg=27700)

    print(f"aquired {len(residences)} residences and {len(hairdressers)} hairdressers")

    add_residences_as_nodes_to_graph(G, residences)

    # folium_map = display_with_folium(location_point, G)

    # # save the folium map to an html file
    # folium_map.save("map.html")
    # print("Map saved to map.html")

if __name__ == "__main__":
    __main__()