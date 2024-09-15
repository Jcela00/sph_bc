import os
from pathlib import Path
import re

# Paths
input_folder = 'XML_TC'  # Folder where the original XML files are located
output_folder = 'XML_TCNew'  # Folder where modified XML files will be saved

# Create output folder if it doesn't exist
Path(output_folder).mkdir(parents=True, exist_ok=True)

# Regular expression to find the bcType attribute
bc_type_pattern = re.compile(r'\s*bcType="[^"]*"')

# Loop over all XML files in the input folder
for filename in os.listdir(input_folder):
    if filename.endswith(".xml"):
        # Full path to input XML file
        input_file = os.path.join(input_folder, filename)

        # Read the content of the file
        with open(input_file, 'r', encoding='utf-8') as file:
            content = file.readlines()

        # Find the line that contains 'bcType="..."' and insert new attributes after it
        for i, line in enumerate(content):
            if bc_type_pattern.search(line):
                # Insert the new attributes after the line with 'bcType'
                content.insert(i + 1, '        bcX="non_periodic" bcY="non_periodic"\n')
                break

        # Write the modified content to a new file in the output folder
        output_file = os.path.join(output_folder, filename)
        with open(output_file, 'w', encoding='utf-8') as file:
            file.writelines(content)

print(f"Modified XML files have been saved to '{output_folder}'.")
