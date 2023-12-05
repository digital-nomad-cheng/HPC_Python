import os
import re
import math
import json

# Function to extract smallest resmax and total time taken from a file
def extract_data(file_path):
    smallest_resmax = float('inf')
    total_time = None

    with open(file_path, 'r') as file:
        lines = file.readlines()

        for line in lines:
            if line.startswith('resmax:'):
                resmax_value = float(line.split(': ')[1])
                if not math.isnan(resmax_value):  # Check if resmax is not NaN
                    smallest_resmax = min(smallest_resmax, resmax_value)
                else:
                    smallest_resmax = float('nan')
            elif re.match(r'\s+\d+ function calls \(.+\) in \d+\.\d+ seconds', line):
                total_time = float(line.split()[-2])
                break

    return smallest_resmax, total_time

# Directory containing your .txt files
folder_path = './'

# Get a list of all .txt files in the directory
txt_files = [file for file in os.listdir(folder_path) if file.endswith('.txt')]

data_dict = {}

for txt_file in txt_files:
    file_path = os.path.join(folder_path, txt_file)
    
    # Extracting variables from the filename
    filename_parts = txt_file.split('_')
    cupy_cpu = filename_parts[0]
    solver = filename_parts[1]
    n = re.findall(r'\d+', filename_parts[3])[0]  # Extracting digits from 'mesh_<n>.txt'
    
    smallest_resmax, total_time = extract_data(file_path)
    
    # Creating the dictionary structure
    if solver not in data_dict:
        data_dict[solver] = {}
    if cupy_cpu not in data_dict[solver]:
        data_dict[solver][cupy_cpu] = {}
    data_dict[solver][cupy_cpu][n] = {
        'time': total_time,
        'smallest_resmax': smallest_resmax
    }

# Sort the dictionary by 'n'
for solver in data_dict:
    for cupy_cpu in data_dict[solver]:
        data_dict[solver][cupy_cpu] = dict(sorted(data_dict[solver][cupy_cpu].items(), key=lambda x: int(x[0])))

# Save the sorted dictionary to a JSON file
output_file = 'parsed_data.json'
with open(output_file, 'w') as json_file:
    json.dump(data_dict, json_file, indent=4)

print(f"Data sorted and saved to {output_file}")
