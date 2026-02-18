from subprocess import run, PIPE
from pathlib import Path
import os
import sys

# --------------------------------------
# Set working directory
# --------------------------------------
"""script_path = Path(__file__).resolve() # get path of current script
cwd = script_path.parent # get directory"""
cwd = Path(sys.argv[1])

# --------------------------------------
# Include ifort
# --------------------------------------
env = os.environ.copy()
ifort_dir = "/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64"
env["PATH"] = ifort_dir + ":" + env["PATH"]

# --------------------------------------
# Run Y-RVE simulations
# --------------------------------------
tissues = ['EW', 'TW', 'LW']
# Run Y-RVE creep simulations for each tissue
for tissue in tissues:
    print("Running {} ...".format(tissue))
    script_path = cwd / ('Calculate_Y-RVE_Equivalent_Compliance_{}.py'.format(tissue)) # path to script file
    output = run("abaqus cae noGUI=" + str(script_path) + " -- " + str(cwd), stdout=PIPE, stderr=PIPE,
                 universal_newlines=True, shell=True, env=env) # run script
print('Done!')
