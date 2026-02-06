import osmnx as ox
import folium
import matplotlib.pyplot as plt
from shapely.geometry import Point
import pandas as pd
import geopandas as gpd
import networkx as nx
import os

""" returns a GeoDataFrame of postcodes and their coordinates"""
def get_postcode_locations(file_path):
    if os.path.exists("data/cardiff_postcodes.geojson"):
        print("Cardiff postcodes GeoJSON file already exists. Loading from file...")
        gdf = gpd.read_file("data/cardiff_postcodes.geojson")
        return gdf
    else:
        print("Cardiff postcodes GeoJSON file not found. Processing ONS postcodes CSV file...")
    
        # read csv file with pandas
        postcode_df = pd.read_csv(file_path)

        # Filter to active postcodes only
        active_postcodes = postcode_df[(postcode_df['doterm'].isna()) | (postcode_df['doterm'] == '')]
        
        # Filter to the postcode's name and its latitude and longitude
        lat_long_df = active_postcodes[["pcds", "lat", "long"]]

        # Create geometry column from lat/long coordinates
        geometry = [Point(xy) for xy in zip(lat_long_df['long'], lat_long_df['lat'])]

        # Create GeoDataFrame with EPSG:4326 (WGS84)
        gdf = gpd.GeoDataFrame(lat_long_df, geometry=geometry, crs='EPSG:4326')

        # Convert to EPSG:27700 (British National Grid)
        gdf = gdf.to_crs('EPSG:27700')

        # save the geodataframe to json file
        gdf.to_file('data/cardiff_postcodes.geojson', driver="GeoJSON")

        return gdf


""" returns a MultiDigraph of the walkable network around a location point"""
def get_street_network_graph(location_point, distance):
    if os.path.exists("data/cardiff_network.graphml"):
        print("Street network graph already exists. Loading from file...")
        G = ox.io.load_graphml("data/cardiff_network.graphml")
        return G
    else:
        print("Street network graph not found. Downloading from OpenStreetMap...")
        G =  ox.graph_from_point(location_point, dist=distance, dist_type="network", network_type="walk")
        G_proj = ox.projection.project_graph(G, to_crs="EPSG:27700")

        ox.distance.add_edge_lengths(G_proj)
        ox.save_graphml(G_proj, filepath="data/cardiff_network.graphml")

        return G_proj


def get_hairdressers(G, location_point, distance):
    # find the nearest hairdresser to the start node and get its node id
    hairdresser_tags = {"shop": "hairdresser"}
    hairdressers_gdf = ox.features_from_point(location_point, tags=hairdresser_tags, dist=distance).to_crs(epsg=27700)

    # get the centroid of the hairdressers and find the nearest node to that centroid
    hairdressers_gdf['centroid'] = hairdressers_gdf.geometry.centroid

    # Find nearest network nodes for all hairdressers (vectorized)
    hairdressers_gdf['nearest_node'] = ox.distance.nearest_nodes(
        G,
        hairdressers_gdf['centroid'].x,
        hairdressers_gdf['centroid'].y
    )

    return hairdressers_gdf


""" finds the shortest route between two nodes in the graph and plots it"""
def get_shortest_route(G, start_node, hairdressers_gdf):
    # find closest nodes to the start point
    nearest_node_id, distance_to_node = ox.distance.nearest_nodes(G, start_node.geometry.iloc[0].x, start_node.geometry.iloc[0].y, return_dist=True)

    # find each hairdresser's distance to the start node
    hairdressers_gdf['distance_m'] = hairdressers_gdf.geometry.distance(start_node.geometry.iloc[0])
    
    # Sort by distance (closest first)
    hairdressers_ranked = hairdressers_gdf.sort_values('distance_m')

    print("Closest hairdressers ranked by distance:")
    print(hairdressers_ranked[['name', 'distance_m']] if 'name' in hairdressers_ranked.columns else hairdressers_ranked[['distance_m']])

    # shortest path through network node for the 10 seemingly closest hairdressers
    route = ox.routing.shortest_path(G, nearest_node_id, hairdressers_ranked.iloc[0]['nearest_node'], weight="length")
    
    # Calculate network distance to each hairdresser
    distances = []
    for idx, row in hairdressers_ranked.iloc[:10].iterrows():
        try:
            # Calculate shortest path length
            distance = ox.routing.shortest_path(
                G, 
                nearest_node_id, 
                row['nearest_node'], 
                weight='length',
                cpus=1
            )
            if distance:
                # Sum up edge lengths along the path
                path_length = sum(ox.utils_graph.get_route_edge_attributes(G, distance, 'length'))
                distances.append(path_length)

                # plot the route
                fig, ax = ox.plot.plot_graph_route(G, distance, route_color="y", route_linewidth=6, node_size=0)
                plt.show()
            else:
                distances.append(float('inf'))  # No path found
        except:
            distances.append(float('inf'))  # Error in routing

    # plot the route
    

def main(location_point, distance, postcode_file_path):
    # get the street network graph and the postcode locations
    G = get_street_network_graph(location_point, distance)
    
    # not using postcode right now
    postcode_gdf = get_postcode_locations(postcode_file_path)

    # CF24 4RT - (Mackintosh Place)
    start_node = postcode_gdf.loc[postcode_gdf['pcds'] == 'CF24 4RT']

    haidressers_gdf = get_hairdressers(G, location_point, distance)

    get_shortest_route(G, start_node, haidressers_gdf)

def __main__():
    location_point = (51.496103, -3.173560) # 115 mackintosh place
    distance = 750  # meters

    postcode_file_path = "data/ons_postcodes/ONSPD_FEB_2023_UK_CF.csv"

    main(location_point, distance, postcode_file_path)


if __name__ == "__main__":
    __main__()