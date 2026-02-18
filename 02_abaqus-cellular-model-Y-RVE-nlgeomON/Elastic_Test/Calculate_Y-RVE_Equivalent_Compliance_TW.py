################# IMPORTS #################
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
umatpath = os.path.join(modulespath, 'LayersMatModel.f') # umat path
geompath = os.path.join(modulespath, 'growth_ring_data.csv') # geometry path

################# IMPORT CLASS #################
import Y_RVE_Class
imp.reload(Y_RVE_Class)
from Y_RVE_Class import Y_RVE_Class
# Create class instance
RVE = Y_RVE_Class()

################# FUNCTIONS #################
def get_geometry_properties(file_path, material):
    # Load csv file and extract geometry properties
    with open(file_path, "r") as f:
        reader = csv.reader(f, delimiter=",")  
        headers = [h.strip() for h in next(reader)]
        col_indices = {name: i for i, name in enumerate(headers)}
        for row in reader:
            if row[0].strip() == material:  # Match material
                return tuple(float(row[col_indices[col]].strip()) for col in ["t [microm]", "Wr [microm]", "Wt1 [microm]", "Wt2 [microm]"])
    raise ValueError("Material '{}' not found.".format(material))
# end: def get_geometry_properties

def write_result_to_csv(result, filename):
    # List of pair identifiers (adjust as necessary)
    pairs = ['11', '22', '33', '44', '55', '66', '12', '13', '23']
    
    # Build header list: ['time11', 'D11', 'time22', 'D22', ...]
    headers = []
    for pair in pairs:
        #headers.append('time' + pair)
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
w = 0.12
################# SET TISSUE TYPE #################
tissue = 'TW'
################# SETUP GEOMETRY #################
t, Wr, Wt1, Wt2 = get_geometry_properties(geompath, tissue)
b = 1.0 # longitudinal thickness [microm]
h = Wt1
if (Wt2-Wt1)!=0:
    theta = atan((Wt2-Wt1)/Wr)
    l = (Wt2-Wt1)/(2.*sin(theta))
else:
    theta = 0.
    l = Wr/2.  
layers_thick = [0.175, 0.175, 0.125, 0.125, 0., 0.035] # thickness of layer ML, P, S1+, S1-, S2, S3 [microm]
tS2 = t - sum(layers_thick) # S2 thickness [microm]
layers_thick[4] = tS2
layers_MFA = [0, 0, 60, -60, 15., -75] # MFA of layer M, P, S1+, S1-, S2, S3 [deg]
# Create RVE
modelName = RVE.setup_multilayers_geometry(l=l, h=h, theta=theta, lt=layers_thick, b=b)

################# SETUP MATERIALS #################
RVE.create_layers_material(layers_MFA)

################# SETUP MESH #################
# mesh inputs
meshSize = 0.5 # global mesh size
minSizeFactor = 0.5 # smallest allowable fraction of global size
# Create mesh
RVE.setup_mesh(meshSize=meshSize, minSizeFactor=minSizeFactor)

################# CALCULATE COMPLIANCE #################
load = [4, 2, 110, 25, 15, 10] # load for each elementary case
step_time = 1 # step time
# Set creep=True to run visco-elastic simulations, creep=False to run pure elastic simulations
creep = False
# Set sim=True if odb files already exist, sim=False to run the simulations
run = True # true to run simulations, false to read existing odb files
compliance_coeff = RVE.compute_3D_compliance_matrix(load, w, run, umatpath, step_time, creep)

################# SAVE RESULTS #################
# Save in Elastic_Test folder
savepath = os.path.join(mainpath, tissue)
if not os.path.exists(savepath): os.makedirs(savepath)
write_result_to_csv(compliance_coeff, os.path.join(savepath, tissue + "_elastic_compliance_coeffs.csv")) # save csv file

# Save in Original_Folder for creep simulations
savepath = os.path.join(os.path.dirname(mainpath), "Original_Folder", tissue)
if not os.path.exists(savepath): os.makedirs(savepath)
write_result_to_csv(compliance_coeff, os.path.join(savepath, tissue + "_elastic_compliance_coeffs.csv")) # save csv file

print("Done!")