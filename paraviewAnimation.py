import glob
import os
import re
from paraview.simple import *

# Configuration
input_pattern = "Custom_new_summation_Re100_450_200_1rf_TriangleEquilateral_1prc_*.pvtp"  # Adjust pattern as needed
property_name = "velocity"  # Property to visualize (replace with the actual property name)
output_directory = "./frames"  # Directory for saving frames
os.makedirs(output_directory, exist_ok=True)  # Create output directory if it doesn't exist

# Image resolution settings
image_width = 3840  # Desired width for high-res output
image_height = 2160  # Desired height for high-res output

# Function to extract numeric part of the filename for sorting
def extract_numeric(filename):
    match = re.search(r"_1prc_(\d+)\.pvtp$", filename)
    return int(match.group(1)) if match else float('inf')  # Use inf for unmatched files

# Get a sorted list of files matching the pattern
file_paths = sorted(glob.glob(input_pattern), key=extract_numeric)
if not file_paths:
    print(f"No files found matching pattern: {input_pattern}")
    exit(1)

# Create a new render view
render_view = CreateView("RenderView")
render_view.ViewSize = [image_width, image_height]  # Match resolution with desired output
render_view.Background = [1.0, 1.0, 1.0]  # White background

# Place the render view on the main monitor (0,0 typically aligns with the primary display)
# render_view.Position = [0, 0]

# Remove the axis indicator
render_view.OrientationAxesVisibility = 0

# Loop through the files
for i, file_path in enumerate(file_paths):
    print(f"Processing file: {file_path}")
    
    # Open the data file
    data = OpenDataFile(file_path)
    
    # Ensure the property to visualize is available
    if property_name not in data.PointData.keys():
        print(f"Property '{property_name}' not found in file: {file_path}")
        Delete(data)
        continue
    
    # Set the active scalar/vector property
    data.PointArrayStatus = [property_name]
    representation = Show(data, render_view)
    
    # Set the representation properties
    representation.Representation = "Points"
    representation.PointSize = 10
    representation.RenderPointsAsSpheres = 1  # Enable rendering points as spheres
    representation.ColorArrayName = [None, property_name]  # Color by the selected property
    representation.LookupTable = GetLookupTableForArray(property_name, 1)
    
    # Use ResetCamera for initial adjustment
    ResetCamera(render_view)
    # GetActiveView().ResetCamera(True, 0.7)
    
    # Tighten the camera view by modifying the parallel scale slightly
    render_view.CameraParallelScale *= 0.55  # Adjust scale to remove excess margins
    
    # Export the current frame with custom resolution
    frame_filename = f"{output_directory}/frame_{i:03d}.png"
    SaveScreenshot(frame_filename, render_view, ImageResolution=[image_width, image_height])
    print(f"Saved frame: {frame_filename}")
    
    # Hide the current data before loading the next
    Hide(data, render_view)
    Delete(data)

# Clean up and finish
Delete(render_view)
print("Animation frames exported.")
