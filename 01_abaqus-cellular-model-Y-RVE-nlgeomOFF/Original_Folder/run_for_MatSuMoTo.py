# --------------------------------------
# Imports
# --------------------------------------
import csv
from subprocess import run
from pathlib import Path
import sys

# --------------------------------------
# Set working directory
# --------------------------------------
script_path = Path(__file__).resolve() # get path of current script
cwd = script_path.parent # get directory
#cwd = Path(sys.argv[1])

# --------------------------------------
# Define functions
# --------------------------------------
def get_gamma_layers(filename):
    # Reads the last data row from an existing CSV file.
    # Args:
    # - filename (str or Path): Path to the CSV file.   
    # Returns:
    # - dict: A dictionary containing the attempt number and gamma values.

    filename = Path(filename)
    if not filename.exists():
        raise FileNotFoundError(f"CSV file '{filename}' does not exist.")
    
    with open(filename, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile)
        try:
            _ = next(reader)  # Read header row (if any)
        except StopIteration:
            raise ValueError("CSV file is empty.")
        
        last_row = None
        for row in reader:
            if row:  # Ignore empty rows.
                last_row = row
        if last_row is None:
            raise ValueError("CSV file does not contain any data rows.")
        
        try:
            gamma1_val = float(last_row[0])
            gamma2_val = float(last_row[1])
            gamma3_val = float(last_row[2])
            gamma4_val = float(last_row[3])
        except Exception as e:
            raise ValueError("Error converting CSV values: " + str(e))
        
        return {
            "gamma1": gamma1_val,
            "gamma2": gamma2_val,
            "gamma3": gamma3_val,
            "gamma4": gamma4_val
        }

def generate_include_files_for_tissues(main_path, include_path, tissues=["EW", "TW", "LW"]):
    # Generates a single combined include file (gamma_values.inc) with conditional logic
    # based on MATNAME. Replaces old per-tissue include file generation.

    # Args:
    # - main_path (str or Path): Directory containing subfolders for each tissue.
    # - include_path (str or Path): Directory where gamma_values.inc will be written.
    # - tissues (list): List of tissue names (e.g., ['EW', 'TW', 'LW'])
    
    main_path = Path(main_path)
    include_path = Path(include_path)
    output_file = include_path / "tissues_gamma_values.inc"

    def get_gamma_lines(csv_path):
        values = get_gamma_layers(csv_path)
        return [
            f"KW_Gamma_i(1) = {values['gamma1']:.5f}D0\n",
            f"KW_Gamma_i(2) = {values['gamma2']:.5f}D0\n",
            f"KW_Gamma_i(3) = {values['gamma3']:.5f}D0\n",
            f"KW_Gamma_i(4) = {values['gamma4']:.5f}D0\n"
        ]

    with open(output_file, 'w') as f:
        f.write("C Combined gamma values from all tissues\n")
        for i, tissue in enumerate(tissues):
            csv_path = main_path / tissue / f"{tissue}_gamma.csv"
            condition = "IF" if i == 0 else "ELSE IF"
            f.write(f"        {condition} (MATNAME == '{tissue}') THEN\n")
            lines = get_gamma_lines(csv_path)
            f.writelines(["        " + line for line in lines])
        f.write("        END IF\n")
    return


###################### STEP 1 ######################
# Run Y-RVE creep simulations for each tissue
script1_path = cwd / 'run_y-rve_simulations.py'
run(["python3", script1_path, cwd], check=True)
# Fit and store gamma_i for each tissue
fit1_path = cwd / "fit_prony_y-rve.py"
run(["python3", fit1_path, cwd], check=True)

###################### STEP 2 ######################
# Generate inc file with gamma_i for each tissue for umat routine
include2_path = cwd / "Modules"
generate_include_files_for_tissues(cwd, include2_path)
# Run GR-RVE creep simulation
script2_path = cwd / 'run_gr-rve_simulations.py'
run(["python3", script2_path, cwd], check=True)
# Fit and store gamma_i for each tissue
fit2_path = cwd / "fit_prony_gr-rve.py"
run(["python3", fit2_path, cwd], check=True)

###################### STEP 3 ######################
# Delete odb files
for odb in cwd.glob('*.odb'):
    try:
        odb.unlink()
    except Exception as e:
        print(f"Failed to delete {odb.name}: {e}")