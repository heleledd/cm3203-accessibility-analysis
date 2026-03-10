import os
import logging
import osmnx as ox
from shapely.geometry import Point, box
import geopandas as gpd

from config import TARGET_CRS, INPUT_DATA_DIR

def get_street_network_graph(bbox):
    """Returns a MultiDiGraph of the walkable network around a bounding box, projected to EPSG:27700."""
    network_path = os.path.join(INPUT_DATA_DIR, 'cardiff_network.graphml')
    if os.path.exists(network_path):
        logging.info("Street network graph already exists. Loading from file...")
        G = ox.io.load_graphml(network_path)
        
        # Ensure it's projected to EPSG:27700
        sample_x = list(G.nodes(data=True))[0][1]['x']
        if sample_x < 1000:
            logging.info("Graph appears to be in WGS84, reprojecting to EPSG:27700...")
            G = ox.projection.project_graph(G, to_crs=TARGET_CRS)
        return G
    else:
        logging.info("Street network graph not found. Downloading from OpenStreetMap...")
        G = ox.graph_from_bbox(bbox, network_type="walk", simplify=True, retain_all=False)
        ox.distance.add_edge_lengths(G)
        G_proj = ox.projection.project_graph(G, to_crs=TARGET_CRS)
        ox.save_graphml(G_proj, filepath=network_path)
        return G_proj


def reproject_bbox(bbox):
    """Reproject a WGS84 bounding box to EPSG:27700."""
    minx, miny, maxx, maxy = bbox
    min_point = Point(minx, miny)
    max_point = Point(maxx, maxy)

    min_point_reprojected = gpd.GeoSeries([min_point], crs='EPSG:4326').to_crs(TARGET_CRS).iloc[0]
    max_point_reprojected = gpd.GeoSeries([max_point], crs='EPSG:4326').to_crs(TARGET_CRS).iloc[0]

    return (
        min_point_reprojected.x,
        min_point_reprojected.y,
        max_point_reprojected.x,
        max_point_reprojected.y
    )


def split_bbox_into_grid(bbox_reprojected, grid_size):
    """Split a projected bounding box into a grid of smaller cells."""
    grid_path = os.path.join(INPUT_DATA_DIR, f'cardiff_grid_cells_{grid_size}m.geojson')
    if os.path.exists(grid_path):
        logging.info("Cardiff grid cells GeoJSON already exists. Loading from file...")
        gdf = gpd.read_file(grid_path)
        return gdf.to_crs(TARGET_CRS)
    else:
        logging.info("Creating Cardiff Grid Cells GeoJSON...")
        minx, miny, maxx, maxy = bbox_reprojected
        grid_cells, grid_id = [], 0

        x = minx
        while x < maxx:
            y = miny
            while y < maxy:
                cell = box(x, y, min(x + grid_size, maxx), min(y + grid_size, maxy))
                grid_cells.append({'grid_id': grid_id, 'geometry': cell})
                grid_id += 1
                y += grid_size
            x += grid_size

        grid_gdf = gpd.GeoDataFrame(grid_cells, crs=TARGET_CRS)
        grid_gdf.to_file(grid_path, driver="GeoJSON")
        return grid_gdf
