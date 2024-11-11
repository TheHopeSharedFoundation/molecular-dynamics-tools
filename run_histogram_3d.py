############################################################################
# Copyright 2024 The Hope Shared Foundation                                #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
# http://www.apache.org/licenses/LICENSE-2.0                               #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
#                                                                          #
# The Hope Shared Foundation, September 2024                               #
# https://hopesharedfoundation.org/                                        #
# Towards Alleviating Suffering                                            #
############################################################################

# Import necessary libraries
import numpy as np
from scipy import stats
import re
import os
import glob 
import argparse
from datetime import datetime
from pathlib import Path

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()

parser = argparse.ArgumentParser(prog='3D-bin alpha sphere and atom positions',description='This program 3D-bins alpha sphere and atom positions',epilog='Copyright 2024 The Hope Shared Foundation')

parser.add_argument('--filepath',help="This is the directory that contains the alpha-sphere PQR files or atom TXT files",required=True)
parser.add_argument('--type',default="alpha",help="This is the type of data being binned.  alpha for alpha-spheres and atom for atom coordinates",required=True)
parser.add_argument('--basename',default=current_time,help="This is the basename of the output file which will be created.",required=True)
parser.add_argument('--cutoff',default=current_time,help="B-factor cutoff value",required=True)

args = parser.parse_args()

################################################################# Set parameters

x_bin_min = -100
x_bin_max = 100
y_bin_min = -100
y_bin_max = 100
z_bin_min = -100
z_bin_max = 100

bin_range = [[x_bin_min,x_bin_max],[y_bin_min,y_bin_max],[z_bin_min,z_bin_max]]
x_number_bins = 200
y_number_bins = 200
z_number_bins = 200
x_bin_increment = abs(x_bin_max - x_bin_min) / x_number_bins
y_bin_increment = abs(y_bin_max - y_bin_min) / y_number_bins
z_bin_increment = abs(z_bin_max - z_bin_min) / z_number_bins
bin_number = [x_number_bins,y_number_bins,z_number_bins]

################################################################# Extract alpha sphere coordinates and van der Waals radii

if args.type == "alpha":
    
    coordinates = []
    radii = []
    atom_index = 0

    with open(args.filepath,'r') as f:
        for file_row in f:
            row = re.search(r"^ATOM\s+\d+\s+(\w+)\s+\w+\s+\d+\s+([-+]?\d+\.\d*e?[-+]?\d*)\s+([-+]?\d+\.\d*e?[-+]?\d*)\s+([-+]?\d+\.\d*e?[-+]?\d*)\s+[-+]?\d+\.\d*e?[-+]?\d*\s+([-+]?\d+\.\d*e?[-+]?\d*)", file_row)
            if row:
                atom_index += 1
                elem=row.group(1)
                x=row.group(2)
                y=row.group(3)
                z=row.group(4)
                vdw=row.group(5)
                coordinates.append([float(x), float(y), float(z)])
                radii.append(float(vdw))

elif args.type == "atom":

    coordinates = []
    radii = []

    with open(args.filepath,'r') as f:
        for file_row in f:
            row = re.search(r"([-+]?\d+\.\d*e?[-+]?\d*)\s+([-+]?\d+\.\d*e?[-+]?\d*)\s+([-+]?\d+\.\d*e?[-+]?\d*)", file_row)
            if row:
                x=row.group(1)
                y=row.group(2)
                z=row.group(3)
                vdw=1
                print(float(x),float(y), float(z), float(vdw))
                coordinates.append([float(x), float(y), float(z)])
                radii.append(float(vdw))
                    
################################################################## Bin alpha sphere coordinates and calculate the mean van der Waals radius for each bin

# Bin alpha sphere coordinates, calculate the mean van der Waals radius for each bin, and output a single PDB file with normalized bin counts in the occupancy column and mean van der Waals in the b-factor column.
bin_size = 1
x_start = 0
y_start = 0
z_start = 0
atom_type = "ATOM"
atom_index = 0
atom_name = "Ne"
alt_loc = ""
res_name = "ALP"
chain = "A"
res_index = 1
insertion_code = ""
occupancy = 0.00
b_factor = 0.00
element_sym = "Ne"

coordinates_arranged = np.array(coordinates)
coordinates_arranged.shape = (len(coordinates),3)

res = stats.binned_statistic_dd(coordinates_arranged, radii, bins=bin_number,range=bin_range, statistic = 'count')
res2 = stats.binned_statistic_dd(coordinates_arranged,radii,bins=bin_number,range=bin_range,statistic ='mean',binned_statistic_result=res)

maximum = np.max(res.statistic)
print(f"MAXIMUM: {maximum}")

output_pdb_file = args.basename + "_3dBinned.pdb"

count = 0

x_old = 0
y_old = 0
z_old = 0
x_new = 0
y_new = 0
z_new = 0

pdb_output = []

with open(output_pdb_file,'w') as pdb_output_handle:
    for x in range(0,x_number_bins):
        x_coord = x_bin_min + (x * x_bin_increment)
        for y in range(0,y_number_bins):
            y_coord = y_bin_min + (y * y_bin_increment)
            for z in range(0,z_number_bins):
                z_coord = z_bin_min + (z * z_bin_increment)
                b_factor = res.statistic[x][y][z] / maximum
                occupancy = res2.statistic[x][y][z]
                if b_factor >= float(args.cutoff):
                    print(b_factor,occupancy)
                    atom_index += 1
                    atom_code = format(atom_index,'#07x')
                    pdb_output_handle.write(f"{atom_type:6s}{atom_code[2:]:5s} {atom_name:^4s}{alt_loc:1s}{res_name:3s} {chain:1s}{res_index:4d}{insertion_code:1s}   {x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}{occupancy:6.2f}{b_factor:6.2f}          {element_sym:>2s}" + "\n")


