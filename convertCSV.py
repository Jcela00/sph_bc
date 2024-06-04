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
pvtp_file_base = sys.argv[1]
count = int(sys.argv[2])

# make dir if not exist
import os
if not os.path.exists("CSV_Data/" + pvtp_file_base):
    os.makedirs("CSV_Data/" +pvtp_file_base)

for i in range(count+1):
    pvtp_file = pvtp_file_base + "_" + str(i) + ".pvtp"
    csv_file = "CSV_Data/" + pvtp_file_base + "/" + pvtp_file_base + "_"+ str(i) + ".csv"
    convert_vtk_to_csv(pvtp_file, csv_file)
