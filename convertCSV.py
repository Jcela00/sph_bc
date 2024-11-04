from paraview.simple import *
import sys
import os

def convert_vtk_to_csv(vtk_file, csv_file):
    # Load the VTK file
    reader = OpenDataFile(vtk_file)
    
    # Retrieve and set up the animation scene for time management
    animationScene = GetAnimationScene()
    animationScene.UpdateAnimationUsingDataTimeSteps()

    # Save data as CSV with time information
    writer = CSVWriter(FileName=csv_file)
    writer.Input = reader
    writer.FieldAssociation = "Point Data"  # or "Cells"
    writer.AddTime = 1  # Enable time information

    # Write each time step individually if time steps are present
    for time in reader.TimestepValues:
        animationScene.TimeKeeper.Time = time  # Set the animation scene to the current time step
        writer.UpdatePipeline(time)  # Write the current time step

# Get filename and range from command line
filename = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

# Create output directory if it doesn't exist
output_dir = "CSV_Data/" + filename
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Loop through the specified range of files
for i in range(start, end + 1):
    pvtp_file = f"{filename}_{i}.pvtp"
    csv_file = f"{output_dir}/file_{i}.csv"
    convert_vtk_to_csv(pvtp_file, csv_file)
