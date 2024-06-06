from paraview.simple import *
import sys

def convert_vtk_to_csv(vtk_file, csv_file):
    # Load the VTK file
    reader = OpenDataFile(vtk_file)

    # Save data as CSV
    writer = CreateWriter(csv_file, reader)
    writer.FieldAssociation = "Point Data" # or "Cells"
    writer.UpdatePipeline()
    # writer.Write()

# get filename from command line
filename = sys.argv[1]
count = int(sys.argv[2])

# make dir if not exist
import os
if not os.path.exists("CSV_Data/" + filename):
    os.makedirs("CSV_Data/" +filename)

for i in range(count+1):
    pvtp_file = filename+ "_" + str(i) + ".pvtp"
    csv_file = "CSV_Data/" + filename + "/file_"+ str(i) + ".csv"
    convert_vtk_to_csv(pvtp_file, csv_file)
