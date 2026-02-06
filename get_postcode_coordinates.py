# make a dictionary of postcodes and their coordinates
import pandas as pd
import networkx as nx
import osmnx as ox

# all the data is in this file, just need to extract the relevant columns and make a dictionary
file_path = "data/ons_postcodes/ONSPD_FEB_2023_UK_CF.csv"

# extract row 1 (the postcode), row 43 (latitude) and row 44 (longitude) and make a dictionary of postcode: (latitude, longitude)
postcode_df = pd.read_csv(file_path)

# Filter to active postcodes only
active_postcodes = postcode_df[(postcode_df['doterm'].isna()) | (postcode_df['doterm'] == '')]

new_df = active_postcodes[["pcds", "lat", "long"]]

# Create a network graph with the postcodes as nodes
G = nx.Graph()
for index, row in new_df.iterrows():
    postcode = row["pcds"]
    latitude = row["lat"]
    longitude = row["long"]
    G.add_node(postcode, pos=(latitude, longitude))

ox.io.save_graph_geopackage(G, filepath="geopackages/postcodilicious.gpkg")
