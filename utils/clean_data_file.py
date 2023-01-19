import os
import re

# This script preprocesses files from the Human Mortality Database so they can be handled by the pandas library
if __name__ == "__main__":

    # Change to the location of your data files
    data_folder = "/home/spsalmon/aging/mortality_data/new_population_all_countries/"

    # Get the list of all files in the data folder
    files = sorted([os.path.join(data_folder, x)
                   for x in os.listdir(data_folder)])

    # Rename those files according to the country name specified in the first line of each file
    for file in files:
        with open(file, 'r') as f:
            content = f.read()
            country_name = content.split(',', 1)[0].replace(' ', '_')
            fname = os.path.join(data_folder, country_name +
                                 "_population_mortality.txt")
            os.rename(file, fname)

    # Iterates through all the files, replacing every multiple spaces and tabs by single spaces
    for file in files:
        # Read in the file
        with open(file, 'r') as f:
            lines = f.readlines()

        lines = lines[2:]
        new_lines = []
        for line in lines:
            line = re.sub('\s+', ' ', line)+"\n"
            new_lines.append(line)
        # Write the file out again
        with open(file, 'w') as f:
            f.writelines(new_lines)
