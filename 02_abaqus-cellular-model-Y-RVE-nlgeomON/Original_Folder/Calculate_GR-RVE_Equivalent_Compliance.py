### IMPORTS ###
#from __future__ import unicode_literals
from math import *
import imp
import sys
import os
import csv

################# SET PATHS #################
mainpath = sys.argv[10]
os.chdir(mainpath)
modulespath = os.path.join(mainpath, 'Modules')
sys.path.append(modulespath)
umatpath = os.path.join(modulespath, 'TissuesMatModel.f') # umat path
geompath = os.path.join(modulespath, 'growth_ring_data.csv') # geometry path

################# IMPORT CLASS #################
import GR_RVE_Class
imp.reload(GR_RVE_Class)
from GR_RVE_Class import GR_RVE_Class
RVE = GR_RVE_Class()

################# FUNCTIONS #################
def get_geometry_properties(file_path):
    # Load CSV file and extract 'layer_thick [mm]' for EW, TW, and LW
    layer_thickness = {}
    with open(file_path, "r") as f:
        reader = csv.reader(f, delimiter=",")
        headers = [h.strip() for h in next(reader)]
        col_indices = {name: i for i, name in enumerate(headers)}       
        if "layer_thick [mm]" not in col_indices:
            raise ValueError("'layer_thick [mm]' column not found in the CSV file.")        
        for row in reader:
            tissue = row[0].strip()  # Assuming tissue type is in the first column
            if tissue in ["EW", "TW", "LW"]:
                layer_thickness[tissue] = float(row[col_indices["layer_thick [mm]"]].strip())
    return [layer_thickness.get(tissue, None) for tissue in ["EW", "TW", "LW"]]
# end: def get_geometry_properties

def write_result_to_csv(result, filename):
    # List of pair identifiers (adjust as necessary)
    pairs = ['11', '22', '33', '44', '55', '66', '12', '13', '23']
    
    # Build header list: ['time11', 'D11', 'time22', 'D22', ...]
    headers = []
    for pair in pairs:
        headers.append('time' + pair)
        headers.append('D' + pair)
    
    # Create a dictionary to hold all columns as lists, and determine max number of rows.
    data_columns = {}
    max_rows = 0
    for key in headers:
        # Retrieve the value from result; default to empty list if not present.
        val = result.get(key, [])
        # If the value isn't a list, convert it to a single-item list.
        if not isinstance(val, list):
            val = [val]
        data_columns[key] = val
        # Update the maximum number of rows.
        if len(val) > max_rows:
            max_rows = len(val)
    
    # Write the CSV file.
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        # Loop over rows index and build a dict for each row.
        for i in range(max_rows):
            row = {}
            for key in headers:
                # If there is no value at index i, fill with empty string.
                row[key] = data_columns[key][i] if i < len(data_columns[key]) else ''
            writer.writerow(row)
    return
# end: write_result_to_csv

################# SET MOISTURE #################
w = 0.07 # moisture content (from 0 to 0.3)
################# SETUP GEOMETRY #################
# Geometry inputs
b = 1.0 # longitudinal thickness [microm]
t = get_geometry_properties(geompath)
# Create part
modelName = RVE.setup_growth_ring_part(t, b)

################# SETUP MATERIALS #################
RVE.setup_tissue_material()

################# SETUP MESH #################
# mesh inputs
meshSize = min(t)/3. # global mesh size
minSizeFactor = 1. # smallest allowable fraction of global size
# Create mesh
RVE.setup_mesh(meshSize=meshSize, minSizeFactor=minSizeFactor)

################# CALCULATE COMPLIANCE #################
load = [1.5e-2, 1.5e-2, 5e-1, 2e-2, 2.3e-2, 1.5e-2] # load for each elementary case
step_time = 150 # step time
# Set creep=True to run visco-elastic simulations, creep=False to run pure elastic simulations
creep = True
# Set sim=True if odb files already exist, sim=False to run the simulations
run = True # true to run simulations, false to read existing odb files
compliance_coeff = RVE.compute_3D_compliance_matrix(load, w, run, umatpath, step_time, creep)

################# SAVE RESULTS #################
tissue = 'GR'
savepath = os.path.join(mainpath, tissue)
if not os.path.exists(savepath):
    os.makedirs(savepath)
# Example of how to call the function:
savepath = os.path.join(mainpath, tissue)
if not os.path.exists(savepath):
    os.makedirs(savepath)
write_result_to_csv(compliance_coeff, os.path.join(savepath, tissue + "_creep_compliance_coeffs.csv"))

print("Done!")
