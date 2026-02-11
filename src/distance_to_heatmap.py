# AI generated script to create an interactive heatmap showing distance to nearest supermarket for Cardiff postcodes using Folium
# uses data generated from the residence_to_feature.py script, which calculates the shortest network distance from each postcode to the nearest supermarket and saves it in a GeoJSON file

import json
import folium
from folium.plugins import HeatMap
import webbrowser
import os
import sys

def load_geojson(filepath):
    """Load GeoJSON data from a file"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            geojson_data = json.load(f)
        print(f"✓ Successfully loaded GeoJSON file: {filepath}")
        print(f"✓ Found {len(geojson_data['features'])} postcodes")
        return geojson_data
    except FileNotFoundError:
        print(f"✗ Error: File not found - {filepath}")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"✗ Error: Invalid JSON format in {filepath}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ Error loading file: {e}")
        sys.exit(1)

def create_heatmap(geojson_data):
    """Create an interactive heatmap showing distance to nearest supermarket for Cardiff postcodes"""
    
    # Extract data for heatmap
    # Format: [latitude, longitude, intensity_value]
    heat_data = []
    
    for feature in geojson_data['features']:
        lat = feature['properties']['lat']
        lon = feature['properties']['long']
        distance = feature['properties']['distance_to_location']
        
        # Higher distance = worse access = higher heat intensity
        # Normalize the distance value for better visualization
        intensity = distance / 100  # Divide by 100 to scale appropriately
        
        heat_data.append([lat, lon, intensity])
    
    # Calculate center point for the map
    avg_lat = sum([f['properties']['lat'] for f in geojson_data['features']]) / len(geojson_data['features'])
    avg_lon = sum([f['properties']['long'] for f in geojson_data['features']]) / len(geojson_data['features'])
    
    # Create base map centered on Cardiff
    m = folium.Map(
        location=[avg_lat, avg_lon],
        zoom_start=14,
        tiles='OpenStreetMap'
    )
    
    # Add the heatmap layer
    # Red areas = far from supermarkets, Blue/Green areas = close to supermarkets
    HeatMap(
        heat_data,
        min_opacity=0.3,
        max_zoom=18,
        radius=25,
        blur=15,
        gradient={
            0.0: 'green',
            0.3: 'yellow',
            0.6: 'orange',
            1.0: 'red'
        }
    ).add_to(m)
    
    # Add markers for each postcode (optional - can be toggled on/off)
    marker_group = folium.FeatureGroup(name='Postcode Markers', show=False)
    
    for feature in geojson_data['features']:
        lat = feature['properties']['lat']
        lon = feature['properties']['long']
        postcode = feature['properties']['pcds']
        distance = feature['properties']['distance_to_location']
        
        # Color code markers based on distance
        if distance < 500:
            color = 'green'
        elif distance < 1000:
            color = 'lightgreen'
        elif distance < 1500:
            color = 'orange'
        else:
            color = 'red'
        
        folium.CircleMarker(
            location=[lat, lon],
            radius=5,
            popup=f"<b>{postcode}</b><br>Distance to supermarket: {distance:.0f}m",
            tooltip=postcode,
            color=color,
            fill=True,
            fillColor=color,
            fillOpacity=0.7
        ).add_to(marker_group)
    
    marker_group.add_to(m)
    
    # Add layer control
    folium.LayerControl().add_to(m)
    
    # Add a custom legend
    legend_html = '''
    <div style="position: fixed; 
                bottom: 50px; right: 50px; width: 220px; height: 160px; 
                background-color: white; border:2px solid grey; z-index:9999; 
                font-size:14px; padding: 10px">
    <p style="margin-bottom: 5px;"><b>Distance to Nearest Supermarket</b></p>
    <p style="margin: 3px;"><i class="fa fa-circle" style="color:green"></i> &lt; 500m (Excellent)</p>
    <p style="margin: 3px;"><i class="fa fa-circle" style="color:yellow"></i> 500-1000m (Good)</p>
    <p style="margin: 3px;"><i class="fa fa-circle" style="color:orange"></i> 1000-1500m (Fair)</p>
    <p style="margin: 3px;"><i class="fa fa-circle" style="color:red"></i> &gt; 1500m (Poor)</p>
    <p style="margin-top: 10px; font-size: 11px;">Toggle markers using layer control ↖</p>
    </div>
    '''
    m.get_root().html.add_child(folium.Element(legend_html))
    
    # Save the map
    output_file = 'cardiff_supermarket_heatmap.html'
    m.save(output_file)
    print(f"✓ Heatmap created successfully: {output_file}")
    
    # Get absolute path
    abs_path = os.path.abspath(output_file)
    print(f"✓ File location: {abs_path}")
    
    # Open in browser
    webbrowser.open('file://' + abs_path)
    print("✓ Opening in your default browser...")
    
    return abs_path

if __name__ == "__main__":
    print("Creating Cardiff Supermarket Access Heatmap...")
    print("-" * 50)
    
    # Check if filename provided as command line argument
    if len(sys.argv) > 1:
        geojson_file = sys.argv[1]
    else:
        # Default filename - change this to your file path
        geojson_file = 'data/filtered_postcodes_with_distances.geojson'
        print(f"No file specified, using default: {geojson_file}")
        print(f"Usage: python {sys.argv[0]} <path_to_geojson_file>")
        print("-" * 50)
    
    # Load the GeoJSON data
    geojson_data = load_geojson(geojson_file)
    
    # Create the heatmap
    create_heatmap(geojson_data)
    
    print("-" * 50)
    print("Done! The map should open in your browser.")
    print("\nHeatmap Legend:")
    print("  🟢 Green areas: Close to supermarkets (<500m)")
    print("  🟡 Yellow areas: Moderate distance (500-1000m)")
    print("  🟠 Orange areas: Further away (1000-1500m)")
    print("  🔴 Red areas: Far from supermarkets (>1500m)")