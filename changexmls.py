import os
import xml.etree.ElementTree as ET

def update_umax_in_file(file_path, new_umax):
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Find the "physics" element and update the "umax" attribute
    for physics in root.findall('simulation'):
        physics.set('time', str(new_umax))

    # Write the modified XML back to the file
    tree.write(file_path, encoding="utf-8", xml_declaration=True)

def update_umax_in_folder(folder_path, new_umax):
    for filename in os.listdir(folder_path):
        if filename.endswith(".xml"):
            file_path = os.path.join(folder_path, filename)
            update_umax_in_file(file_path, new_umax)
            print(f"Updated {filename}")

# Example usage:
folder_path = "XMLs/XML_TC/"  # Replace with the actual folder path
new_umax = 20
update_umax_in_folder(folder_path, new_umax)