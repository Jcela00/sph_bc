import os
import xml.etree.ElementTree as ET

def update_value_on_file(file_path, category, field, new_value):
    # Read the original XML content to preserve its formatting
    with open(file_path, "r", encoding="utf-8") as f:
        xml_content = f.read()

    # Parse XML content without reformatting
    tree = ET.ElementTree(ET.fromstring(xml_content))
    root = tree.getroot()

    # Update the specified field in the given category
    for elem in root.findall(category):
        old_field_value = f'{field}="{elem.get(field)}"'
        new_field_value = f'{field}="{new_value}"'
        xml_content = xml_content.replace(old_field_value, new_field_value)

    # Write back the updated content
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(xml_content)

def update_value_on_folder(folder_path, category, field, new_value):
    for filename in os.listdir(folder_path):
        if filename.endswith(".xml"):
            file_path = os.path.join(folder_path, filename)
            update_value_on_file(file_path, category, field, new_value)
            print(f"Updated {filename}")

# Example usage:
folder_path = "XML/poiseuilleDiff/"  # Replace with the actual folder path
new_val = 1.0
update_value_on_folder(folder_path, 'physics', 'Bfactor', new_val)
