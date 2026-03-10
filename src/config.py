import os

TARGET_CRS = 'EPSG:27700'

# Cardiff bounding box (west, south, east, north) in WGS84
# CARDIFF_BBOX = (-3.278, 51.447, -3.066, 51.554) # coordinates for all of Cardiff
# CARDIFF_BBOX = (-3.196, 51.494, -3.180, 51.500) # smaller box for testing

bbox_string = os.environ.get('BBOX', '-3.196, 51.494, -3.180, 51.500')
CARDIFF_BBOX = tuple(float(coord.strip()) for coord in bbox_string.split(','))

GRID_SIZE_METERS = 100

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_DATA_DIR = os.path.join(SCRIPT_DIR, '../data/input_data')
OUTPUT_DATA_DIR = os.path.join(SCRIPT_DIR, '../data/output')

access_points_filename = os.environ.get('ACCESS_POINTS_FILE', 'cardiff_greenspace_access_points.geojson')
boundary_filename = os.environ.get('BOUNDARY_FILE', 'cardiff_greenspace_sites.geojson')

PARK_ACCESS_POINTS_PATH = os.path.join(INPUT_DATA_DIR, access_points_filename)
PARK_BOUNDARY_PATH = os.path.join(INPUT_DATA_DIR, boundary_filename)
