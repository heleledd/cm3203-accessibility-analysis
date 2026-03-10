import os
import logging
import pandas as pd
import numpy as np
import osmnx as ox
import pandana
import geopandas as gpd
from shapely.geometry import Point, box

from config import TARGET_CRS, INPUT_DATA_DIR

import warnings
warnings.filterwarnings('ignore')

pd.options.display.float_format = '{:.2f}'.format

def get_street_network_graph(bbox):
    """Returns a projected Pandana network (EPSG:27700) and caches it to CSV."""
    node_path = os.path.join(INPUT_DATA_DIR, 'street_network_nodes.csv')
    edge_path = os.path.join(INPUT_DATA_DIR, 'street_network_edges.csv')
    
    # --- AUTO-PURGE OLD CACHE ---
    # We are changing how the CSV is saved, so we must delete the old ones.
    if os.path.exists(node_path):
        check_nodes = pd.read_csv(node_path, nrows=5)
        # If the old 'osmid' index is still in the file, delete it
        if 'osmid' in check_nodes.columns or check_nodes.index.max() > 2000000:
            logging.warning("Incompatible cache detected. Auto-deleting and rebuilding...")
            os.remove(node_path)
            if os.path.exists(edge_path):
                os.remove(edge_path)
    
    # --- LOAD OR DOWNLOAD DATA ---
    if os.path.exists(node_path) and os.path.exists(edge_path):
        logging.info("Loading pre-projected street network from CSV...")
        # Notice we are no longer using index_col=0!
        nodes = pd.read_csv(node_path)
        edges = pd.read_csv(edge_path)
        
    else:
        logging.info("Downloading network via OSMnx and projecting to EPSG:27700...")
        G = ox.graph_from_bbox(bbox, network_type="walk", simplify=True)
        G_proj = ox.projection.project_graph(G, to_crs=TARGET_CRS)
        nodes, edges = ox.graph_to_gdfs(G_proj)
        
        # Map huge 64-bit OSM IDs to safe integers starting from 0
        node_mapping = {osmid: np.int32(i) for i, osmid in enumerate(nodes.index)}
        
        # Replace u and v with mapped values
        edges = edges.reset_index()
        edges['u'] = edges['u'].map(node_mapping)
        edges['v'] = edges['v'].map(node_mapping)
        
        # Drop the old 64-bit OSM ID index completely
        nodes = nodes.reset_index(drop=True)
        
        # Save without an index!
        nodes.to_csv(node_path, index=False)
        edges.to_csv(edge_path, index=False)
        
    # --- THE FINAL FIX: PURE DATA ---
    # We enforce strict 32-bit ints and 64-bit floats on raw Series
    
    x = nodes['x'].astype(np.float64)
    y = nodes['y'].astype(np.float64)
    
    u = edges['u'].astype(np.int32)
    v = edges['v'].astype(np.int32)
    length = edges[['length']].astype(np.float64)

    # Explicitly lock the simple row-number indices to 32-bit to satisfy Cython
    x.index = x.index.astype(np.int32)
    y.index = y.index.astype(np.int32)
    length.index = length.index.astype(np.int32)
    
    # Initialize network
    network = pandana.Network(x, y, u, v, length)
    return network

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