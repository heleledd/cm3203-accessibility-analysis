import osmnx as ox
import folium
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx


def get_network_graph(location_point, distance):
    G =  ox.graph_from_point(location_point, dist=distance, dist_type="network", network_type="walk")
    G_proj = ox.projection.project_graph(G, to_crs="EPSG:27700")

    ox.distance.add_edge_lengths(G)

    # all the data is in this file, just need to extract the relevant columns and make a dictionary
    file_path = "data/ons_postcodes/ONSPD_FEB_2023_UK_CF.csv"

    # extract row 1 (the postcode), row 43 (latitude) and row 44 (longitude) and make a dictionary of postcode: (latitude, longitude)
    postcode_df = pd.read_csv(file_path)

    # Filter to active postcodes only and then only to the columns we need
    active_postcodes = postcode_df[(postcode_df['doterm'].isna()) | (postcode_df['doterm'] == '')]
    new_df = active_postcodes[["pcds", "lat", "long"]]

    # convert the coordinates to the same CRS as the graph


    # Create a network graph with the postcodes as nodes
    for index, row in new_df.iterrows():
        postcode = row["pcds"]
        latitude = row["lat"]
        longitude = row["long"]
        G_proj.add_node(len(G_proj.nodes()) + 1, y=latitude, x=longitude, postcode=True)  
    # maybe project the coordinates to the same CRS as the graph and add them as attributes instead of using y and x?


    # # try adding some nodes with networkx (address to Deborah's shop)
    # G.add_node(max(G.nodes()) + 1, y=51.496490, x=-3.167884)

    ox.io.save_graph_geopackage(G_proj, filepath="geopackages/postcodeilicious.gpkg")
    return G_proj

def get_shortest_route(G, start_node, end_node):
    route = ox.routing.shortest_path(G, start_node, end_node, weight="length")
    fig, ax = ox.plot.plot_graph_route(G, route, route_color="y", route_linewidth=6, node_size=0)
    plt.show()


def __main__():
    location_point = (51.496103, -3.173560) # 115 mackintosh place
    distance = 750  # meters

    G = get_network_graph(location_point, distance)

    # # get residential buildings
    # residences_tags = {"building": True, "residential": True}
    # residences = ox.features_from_point(location_point, tags=residences_tags, dist=distance).to_crs(epsg=27700) 

    # # get hairdressers
    # shop_tags = {"shop": "hairdresser"}
    # hairdressers = ox.features_from_point(location_point, tags=shop_tags, dist=distance).to_crs(epsg=27700)

    # print(f"aquired {len(residences)} residences and {len(hairdressers)} hairdressers")

    # orig = next(iter(G))
    # dest = list(G)[36] # picked a random node
    # get_shortest_route(G, orig, dest)



if __name__ == "__main__":
    __main__()