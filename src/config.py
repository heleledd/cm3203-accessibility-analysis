import os

TARGET_CRS = 'EPSG:27700'

# Cardiff bounding box (west, south, east, north) in WGS84
# CARDIFF_BBOX = (-3.35, 51.37, -3.05, 51.57) # coordinates for all of Cardiff

CARDIFF_BBOX = (-3.180, 51.494, -3.196, 51.500) # smaller box for testing

GRID_SIZE_METERS = 100

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_DATA_DIR = os.path.join(SCRIPT_DIR, '../data/input_data')
OUTPUT_DATA_DIR = os.path.join(SCRIPT_DIR, '../output')

PARK_ACCESS_POINTS_PATH = os.path.join(INPUT_DATA_DIR, 'cardiff_greenspace_access_points.geojson')
PARK_BOUNDARY_PATH = os.path.join(INPUT_DATA_DIR, 'cardiff_greenspace_sites.geojson')

# TODO: potentially add argparse here to customise the script