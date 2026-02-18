# --------------------------------------
# Imports
# --------------------------------------
from subprocess import run, PIPE
from pathlib import Path
import os
import sys
import csv

# --------------------------------------
# Set working directory
# --------------------------------------
#script_path = Path(__file__).resolve() # get path of current script
#cwd = script_path.parent # get working directory
cwd = Path(sys.argv[1]) # get working directory from input arguments

# --------------------------------------
# Include ifort
# --------------------------------------
env = os.environ.copy()
ifort_dir = "/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64"
env["PATH"] = ifort_dir + ":" + env["PATH"]

# --------------------------------------
# Functions
# --------------------------------------
def save_eng_const_to_inc(savepath, tissues):
    #savepath = os.path.abspath(savepath)
    module_path = os.path.join(savepath, "Modules")
    if not os.path.exists(module_path):
        os.makedirs(module_path)

    inc_file = os.path.join(module_path, "eng_const.inc")

    def read_constants(tissue):
        csv_file = os.path.join(savepath, tissue, f"{tissue}_elastic_compliance_coeffs.csv")
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            row = next(reader)

        S11 = float(row["D11"])
        S22 = float(row["D22"])
        S33 = float(row["D33"])
        S44 = float(row["D44"])
        S55 = float(row["D55"])
        S66 = float(row["D66"])
        S12 = float(row["D12"])
        S13 = float(row["D13"])
        S23 = float(row["D23"])

        E1 = 1.0 / S11
        E2 = 1.0 / S22
        E3 = 1.0 / S33
        G23 = 1.0 / S44
        G13 = 1.0 / S55
        G12 = 1.0 / S66
        nu12 = -S12 / S11
        nu13 = -S13 / S11
        nu23 = -S23 / S22

        return {
            "KW_E1": E1, "KW_E2": E2, "KW_E3": E3,
            "KW_nu23": nu23, "KW_nu13": nu13, "KW_nu12": nu12,
            "KW_G23": G23, "KW_G13": G13, "KW_G12": G12
        }

    def format_block(constants):
        lines = []
        for key, value in constants.items():
            lines.append("      {0} = ({1:.3f}D0)\n".format(key, value))
        return lines

    with open(inc_file, 'w') as f:
        f.write("C Engineering constants per material\n")

        for i, tissue in enumerate(tissues):
            clause = "IF" if i == 0 else "ELSE IF"
            f.write("      {0} (MATNAME == '{1}') THEN\n".format(clause, tissue))
            f.writelines(format_block(read_constants(tissue)))

        # Hardcoded TEST case
        f.write("      ELSE IF (MATNAME == 'TEST') THEN\n")
        f.write("""\
      KW_E1 = (1.0D0)
      KW_E2 = (1.0D0)
      KW_E3 = (1.0D0)
      KW_nu23 = (0.3D0)
      KW_nu13 = (0.3D0)
      KW_nu12 = (0.3D0)
      KW_G23 = (0.3846D0)
      KW_G13 = (0.3846D0)
      KW_G12 = (0.3846D0)\n""")

        # Fallback for unknown MATNAME
        f.write("      ELSE\n")
        f.write("      PRINT *, 'Error: Material \"', MATNAME, '\" not supported.'\n")
        f.write("      RETURN\n")
        f.write("      END IF\n")# end: def save_eng_const_to_inc

# --------------------------------------
# Run Y-RVE simulations
# --------------------------------------
tissues = ['EW', 'TW', 'LW'] # list of tissue types
# Run creep simulations for each tissue
for tissue in tissues:
    print("Running {} ...".format(tissue))
    script_path = cwd / ('Calculate_Y-RVE_Equivalent_Compliance_{}.py'.format(tissue)) # path to script file
    output = run("abaqus cae noGUI=" + str(script_path) + " -- " + str(cwd), stdout=PIPE, stderr=PIPE,
                 universal_newlines=True, shell=True, env=env) # run script

# Save engineering constants
save_eng_const_to_inc(cwd, tissues)
savepath = os.path.join(os.path.dirname(cwd), "Original_Folder")
save_eng_const_to_inc(savepath, tissues)

print('Done!')
