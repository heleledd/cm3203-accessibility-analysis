import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt


location_point = (51.496103, -3.173560)
distance = 750  # meters

G =  ox.graph_from_point(location_point, dist=distance, dist_type="network", network_type="walk")

G_proj = ox.projection.project_graph(G)


G2 = ox.simplification.consolidate_intersections(
    G_proj,
    rebuild_graph=True,
    tolerance=15,
    dead_ends=True,
)

ox.io.save_graph_geopackage(G2, filepath="mack_less_nodes.gpkg")

# get residential buildings
residential_buildings_tags = {"building": True, "residential": True}
shop_tags = {"shop": True}

residential_buildings = ox.features_from_point(location_point, tags=residential_buildings_tags, dist=distance)
shops = ox.features_from_point(location_point, tags=shop_tags, dist=distance)

#--- Project both to match the graph CRS ---
residential_buildings_proj = residential_buildings.to_crs(epsg=32630)  # or match G_proj's CRS
shops_proj = shops.to_crs(epsg=32630)

residential_buildings_proj.to_file("residential_buildings.gpkg", driver="GPKG")
shops_proj.to_file("shops.gpkg", driver="GPKG")
