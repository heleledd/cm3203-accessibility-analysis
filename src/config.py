import os

TARGET_CRS = 'EPSG:27700'

# Cardiff bounding box (west, south, east, north) in WGS84
CARDIFF_BBOX = (-3.35, 51.37, -3.05, 51.57)

GRID_SIZE_METERS = 100

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_DATA_DIR = os.path.join(SCRIPT_DIR, '../data/input_data')
OUTPUT_DATA_DIR = os.path.join(SCRIPT_DIR, '../output')
