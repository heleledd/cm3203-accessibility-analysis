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
    if os.path.exists("../data/cardiff_postcodes.geojson"):
        print("Cardiff postcodes GeoJSON file already exists. Loading from file...")
        gdf = gpd.read_file("../data/cardiff_postcodes.geojson")
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
        gdf.to_file('../data/cardiff_postcodes.geojson', driver="GeoJSON")

        return gdf


""" returns a MultiDigraph of the walkable network around a location point"""
def get_street_network_graph(location_point, distance):
    if os.path.exists("../data/cardiff_network.graphml"):
        print("Street network graph already exists. Loading from file...")
        G = ox.io.load_graphml("../data/cardiff_network.graphml")
        return G
    else:
        print("Street network graph not found. Downloading from OpenStreetMap...")
        G =  ox.graph_from_point(location_point, dist=distance, dist_type="network", network_type="walk")
        
        # add edge lengths to the graph (instructions say to do this before projecting or simplifying the graph)
        ox.distance.add_edge_lengths(G)
        
        G_proj = ox.projection.project_graph(G, to_crs="EPSG:27700")

        ox.save_graphml(G_proj, filepath="../data/cardiff_network.graphml")

        return G_proj


""" download features from OpenStreetMap and find the nearest network node to each feature"""
def get_features(G, location_point, distance, tags={"shop": "hairdresser"}):
    # find the nearest features to the start node and get its node id
    features_gdf = ox.features_from_point(location_point, tags=tags, dist=distance).to_crs(epsg=27700)

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
def get_shortest_route(G, start_node, features_gdf):
    # find closest nodes to the start point
    nearest_node_id, distance_to_node = ox.distance.nearest_nodes(G, start_node.geometry.iloc[0].x, start_node.geometry.iloc[0].y, return_dist=True)

    # find each feature's distance to the start node
    features_gdf['distance_m'] = features_gdf.geometry.distance(start_node.geometry.iloc[0])
    
    # Sort by distance (closest first)
    features_ranked = features_gdf.sort_values('distance_m')

    print("Closest features ranked by distance:")
    # print(features_ranked[['name', 'distance_m']] if 'name' in features_ranked.columns else features_ranked[['distance_m']])

    # Calculate distance through the network to the top 5 geographically closest results
    distances = {}
    for idx, row in features_ranked.iloc[:5].iterrows():
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

            print(f"Route length from start node {nearest_node_id} to feature node {row['nearest_node']}: {route_length}")

            distances[row['nearest_node']] = route_length
        except nx.NetworkXNoPath:
            print(f"No path found to feature: {row['name'] if 'name' in row else idx}")
            distances[row['nearest_node']] = float('inf')  # Assign infinity if no path is found
            continue

    
    # discover which one was the shortest route
    if distances:
        closest_feature_node = min(distances, key=distances.get)
        print(f"Closest feature node based on route length: {closest_feature_node} with distance {distances[closest_feature_node]} metres")
        
        # return shortest distance to a feature
        return distances[closest_feature_node]
    else:
        print("No features found within the specified distance.")
        return None


""" Filters the postcode GeoDataFrame to only include postcodes within a certain distance of a location point"""
def get_postcodes_within_distance(postcode_gdf, location_point, distance):
    # Convert location point to a Shapely Point and project to the same CRS as the postcode GeoDataFrame
    location_point_geom = Point(location_point[1], location_point[0]) 
    location_point_proj = gpd.GeoSeries([location_point_geom], crs='EPSG:4326').to_crs(postcode_gdf.crs).iloc[0]
    
    # Calculate distance from location to all postcodes
    postcode_gdf['distance_to_location'] = postcode_gdf.geometry.distance(location_point_proj)
    
    # Filter to postcodes within distance
    postcodes_within = postcode_gdf[postcode_gdf['distance_to_location'] <= distance]
    
    return postcodes_within


def main(location_point, distance, postcode_file_path):
    # get the street network graph and the postcode locations
    G = get_street_network_graph(location_point, distance)
    
    # not using postcode right now
    postcode_gdf = get_postcode_locations(postcode_file_path)

    # filter the postcodes to only those within the distance of the location point (in keeping with the network we're getting from OSM)
    filtered_postcodes = get_postcodes_within_distance(postcode_gdf, location_point, distance)

    # prepare a column to store the shortest network distance (metres) to a supermarket
    filtered_postcodes['shortest_supermarket_m'] = float('nan')

    # get the supermarket features
    supermarket_tags = {"shop": "supermarket"}
    features_gdf = get_features(G, location_point, distance, tags=supermarket_tags)

    
    # iterate through the filtered postcodes and find the shortest distance to a feature for each one
    for idx, postcode in filtered_postcodes.iterrows():
        print(f"\nFinding closest feature for postcode: {postcode['pcds']} at location {postcode.geometry}")
        shortest_distance = get_shortest_route(G, filtered_postcodes.loc[[idx]], features_gdf)
        
        # store the returned shortest distance (None if not found)
        filtered_postcodes.at[idx, 'shortest_supermarket_m'] = shortest_distance

    # save the results to a new GeoJSON file
    filtered_postcodes.to_file('../data/filtered_postcodes_with_distances.geojson', driver='GeoJSON')


def __main__():
    location_point = (51.496103, -3.173560) # 115 mackintosh place
    distance = 2000  # metres

    postcode_file_path = "../data/ons_postcodes/ONSPD_FEB_2023_UK_CF.csv"

    main(location_point, distance, postcode_file_path)


if __name__ == "__main__":
    __main__()
