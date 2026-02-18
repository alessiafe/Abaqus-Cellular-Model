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
umatpath = os.path.join(modulespath, 'TracheidsMatModel.f') # umat path
structpath = os.path.join(modulespath, 'Structure_mat') # mat folder path

################# IMPORT CLASS #################
import Tissue_Class
imp.reload(Tissue_Class)
from Tissue_Class import Tissue_Class
# Create class instance
tiss = Tissue_Class()

################# FUNCTIONS #################
def write_result_to_csv(result, filename):
    # Define mapping from TRL keys to RTL keys
    trl_to_rtl = {
        '11': '22',
        '22': '11',
        '44': '55',
        '55': '44',
        '13': '23',
        '23': '13',
        '33': '33',
        '12': '12',
        '66': '66',
    }

    # Target order of pairs in RTL
    rtl_pairs = ['11', '22', '33', '44', '55', '66', '12', '13', '23']
    
    # Build headers for time and D values in RTL
    headers = []
    for pair in rtl_pairs:
        headers.append('time' + pair)
        headers.append('D' + pair)

    # Prepare remapped result dictionary
    remapped_result = {}
    for pair in rtl_pairs:
        # Determine the corresponding TRL pair
        trl_pair = trl_to_rtl[pair]
        
        # Handle 'Dxx' values
        d_key_rtl = 'D' + pair
        d_key_trl = 'D' + trl_pair
        val_d = result.get(d_key_trl, [])
        remapped_result[d_key_rtl] = val_d if isinstance(val_d, list) else [val_d]

        # Handle 'timexx' values
        t_key_rtl = 'time' + pair
        t_key_trl = 'time' + trl_pair
        val_t = result.get(t_key_trl, [])
        remapped_result[t_key_rtl] = val_t if isinstance(val_t, list) else [val_t]

    # Determine max number of rows
    max_rows = max(len(v) for v in remapped_result.values())

    # Write to CSV
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for i in range(max_rows):
            row = {}
            for key in headers:
                row[key] = remapped_result[key][i] if i < len(remapped_result[key]) else ''
            writer.writerow(row)
    return
# end: write_result_to_csv

################# SET MOISTURE #################
w = 0.12
################# SET TISSUE TYPE #################
tissue = 'TW'
geompath = os.path.join(structpath, 'testspruce_'+tissue+'3_mstruct') # mat file path
caepath = os.path.join(structpath, 'testspruce_'+tissue+'3.cae') # cae file path
################# SETUP GEOMETRY #################
b = 1.0 # longitudinal thickness [microm]
# Create tissue
#modelName = tiss.setup_geometry(geompath, b)
#tiss.save_model(caepath)
modelName = tiss.open_model(caepath, geompath)

################# SETUP MATERIALS #################
tiss.setup_tissue_material(tissue)

################# SETUP MESH #################
# Create surface partition
tiss.create_surface_partition(0.5, sides=True)
# Mesh inputs
meshSizeFiber = 2.5 # global size of fibers
minSizeFiber = 0.9 # min size factor of fibers
meshSizeMatrix = 1. # global size of matrix
minSizeMatrix = 1. # min size factor of matrix
# Create mesh
tiss.create_mesh(meshSizeFiber, meshSizeMatrix, minSizeFiber, minSizeMatrix)

################# CALCULATE COMPLIANCE #################
load = [1.5, 1.3, 13500, 1400, 2500, 1.2] # load for each elementary case [mN]
step_time = 150 # step time
# Set creep=True to run visco-elastic simulations, creep=False to run pure elastic simulations
creep = True
# Set sim=True if odb files already exist, sim=False to run the simulations
run = True # true to run simulations, false to read existing odb files
compliance_coeff = tiss.compute_3D_compliance_matrix(load, w, run, umatpath, step_time, creep)
print(compliance_coeff)

################# SAVE RESULTS #################
savepath = os.path.join(mainpath, tissue)
if not os.path.exists(savepath):
    os.makedirs(savepath)
# Example of how to call the function:
savepath = os.path.join(mainpath, tissue)
if not os.path.exists(savepath):
    os.makedirs(savepath)
write_result_to_csv(compliance_coeff, os.path.join(savepath, tissue + "_creep_compliance_coeffs.csv"))

print("Done!")
