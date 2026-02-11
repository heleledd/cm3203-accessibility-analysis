import folium
import osmnx as ox

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


def folium_chloropleth_tutorial():
    ...

def main():
    
    folium_chloropleth_tutorial()


def __main__():
    location_point = (51.496103, -3.173560) # 115 mackintosh place
    distance = 2000  # metres

    main()


if __name__ == "__main__":
    __main__()



# folium_map = display_with_folium(location_point, G)

# # save the folium map to an html file
# folium_map.save("map.html")
# print("Map saved to map.html")