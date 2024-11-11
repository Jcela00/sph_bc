from paraview.simple import *
import sys
import os

def convert_vtk_to_csv_old(vtk_file, csv_file):
    # Load the VTK file
    reader = OpenDataFile(vtk_file)

    # Save data as CSV
    writer = CreateWriter(csv_file, reader)
    writer.FieldAssociation = "Point Data"  # or "Cells"
    writer.UpdatePipeline()

    # Cleanup
    Delete(reader)
    Delete(writer)

def convert_vtk_to_csv_new(vtk_file, csv_file):
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
    if hasattr(reader, 'TimestepValues'):  # Check if the file has time steps
        for time in reader.TimestepValues:
            animationScene.TimeKeeper.Time = time  # Set the animation scene to the current time step
            writer.UpdatePipeline(time)
    else:
        # If no time steps are available, write once
        writer.UpdatePipeline()

    # Cleanup
    Delete(reader)
    Delete(writer)

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
    
    if "probes" in pvtp_file:
        convert_vtk_to_csv_old(pvtp_file, csv_file)
    else:
        convert_vtk_to_csv_new(pvtp_file, csv_file)
