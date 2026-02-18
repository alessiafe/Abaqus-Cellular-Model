# --------------------------------------
# Imports
# --------------------------------------
from subprocess import run
from pathlib import Path

# --------------------------------------
# Set working directory
# --------------------------------------
script_path = Path(__file__).resolve() # get path of current script
cwd = script_path.parent # get working directory

###################### STEP 1 ######################
# --------------------------------------
# Run Y-RVE creep simulations for each tissue
# --------------------------------------
script1_path = cwd / 'run_tissue_simulations.py' # path to script file
run(["python3", script1_path, cwd], check=True)  # run script

###################### STEP 2 ######################
# --------------------------------------
# Run GR-RVE creep simulation
# --------------------------------------
script2_path = cwd / 'run_gr-rve_simulations.py' # path to script file
run(["python3", script2_path, cwd], check=True) # run script
