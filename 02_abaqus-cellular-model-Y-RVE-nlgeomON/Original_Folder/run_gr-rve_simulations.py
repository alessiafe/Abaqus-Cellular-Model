from subprocess import run, PIPE
from pathlib import Path
import os
import sys

# --------------------------------------
# Set working directory
# --------------------------------------
#script_path = Path(__file__).resolve() # get path of current script
cwd = Path(sys.argv[1]) #script_path.parent # get directory

# --------------------------------------
# Include ifort
# --------------------------------------
env = os.environ.copy()
ifort_dir = "/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64"
env["PATH"] = ifort_dir + ":" + env["PATH"]

# --------------------------------------
# Run GR-RVE simulations
# --------------------------------------
print("Running GR ...")
script_path = cwd / ('Calculate_GR-RVE_Equivalent_Compliance.py') # path to script file
output = run("abaqus cae noGUI=" + str(script_path) + " -- " + str(cwd), stdout=PIPE, stderr=PIPE,
                universal_newlines=True, shell=True, env=env) # run script
print('Done!')