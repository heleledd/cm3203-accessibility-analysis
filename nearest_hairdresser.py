import osmnx as ox
import folium
import random

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


def get_nearest_hairdresser(G, residential_buildings):
    print('hellow yellow')

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

    folium_map = display_with_folium(location_point, G)

    # save the folium map to an html file
    folium_map.save("map.html")
    print("Map saved to map.html")

if __name__ == "__main__":
    __main__()