# Calculation of elastic properties of Y-RVE tissues and growth ring

This folder contains all scripts and files used for the **calculation of the elastic properties (compliances and engineering constants)** of the **Y-RVE tissues** (EW, TW, LW), as well as the **growth ring (GR)** composed of them. All simulations are conducted with non-linear geometry deactivated (**nlgeom = OFF**).

The global coordinate system follows the anatomical directions of wood tissues:  
**x** → radial, **y** → tangential, **z** → longitudinal.

## Folder structure and workflow overview

- **`run_elastic_test.py`**  
  Main script to be **run from the terminal**. It sequentially launches:

  1. **`run_y-rve_simulations.py`**  
     Computes the elastic properties of the individual Y-RVE tissues (EW, TW, LW).  
     It calls the corresponding scripts for each tissue:

     - **`Calculate_Y-RVE_Equivalent_Compliance_EW.py`**
     - **`Calculate_Y-RVE_Equivalent_Compliance_TW.py`**
     - **`Calculate_Y-RVE_Equivalent_Compliance_LW.py`**

     Each script performs Abaqus simulations under **nine elementary load-control cases**  
     (3 axial, 3 shear, 3 biaxial). From these, the equivalent compliance tensor and engineering constants are computed. Results are stored in the current folder and in the `Original_Folder` (for later creep simulations).

  2. **`run_gr-rve_simulations.py`**  
     Computes the equivalent elastic properties of the growth ring based on the previously calculated Y-RVE properties. It launches **`Calculate_GR-RVE_Equivalent_Compliance.py`**, which performs the Abaqus simulations and saves results as above.

  All simulations are performed at **constant moisture content**. Therefore, for each script `Calculate_...py`, the moisture content must be set as user input, and **the same value must be used for all tissues and for the GR-RVE** (currently 0.12).

- **`Collect_engineering_constants.ipynb`**  
  Collects and organizes the engineering constants from all simulations into a summary table stored as **`collect_engineering_constants.csv`**.

- **`Fit_tissues_engineering_constants_for_UMAT.ipynb`**  
  Fits the engineering constants over multiple moisture contents (if available) and prepares a Fortran-formatted text file for UMAT implementation. *(This workflow is no longer used)*

- **`Modules`**  
  Contains all classes, Python libraries, Fortran source files, and include files required to run the Y-RVE and GR-RVE simulations:

  - **`growth_ring_data.csv`**  
  File containing geometric data required to generate the Y-RVE models and the growth ring.

  - **`Y_RVE_Class.py`**  
  Python class defining the 3D Y-RVE geometry (multi-layer cell wall model), boundary conditions, loading cases, and data extraction routines.

  - **`GR_RVE_Class.py`**  
  Python class defining the 3D growth-ring geometry, boundary conditions, loading cases, and data extraction routines.

  - **`PBCfunction.py`**  
  Library containing the periodic boundary condition (PBC) functions.

  - **`TissuesMatModel.f`**  
  Fortran material model implementing the constitutive behaviour of the EW, TW, and LW tissues for the growth ring model.

  - **`LayersMatModel.f`**  
  Fortran material model implementing the constitutive behaviour of the multi-layer cell wall model used for the Y-RVE tissues.

  - **`eng_const.inc`**  
  Include file storing the equivalent engineering constants of the tissues for the growth ring simulations.

  - **`activate_ve_creep.inc`**  
  Include file controlling the activation of viscoelastic/creep behaviour in the Fortran material models. *Must remain disabled for elastic tests.*

  - **`tissues_gamma_values.inc`**  
  Include file storing the gamma values for the tissues (same set for all tissues) used in the growth ring simulations. Not used in the elastic simulations, but required for the Fortran material routine to run (the same routine used for creep simulations).

  - **`layers_gamma_values.inc`**  
  Include file storing the gamma values for the cell-wall layers (same set for all layers) used in the Y-RVE simulations. Not used in the elastic simulations, but required for the Fortran material routine to run (the same routine used for creep simulations).

## Note:
In Abaqus, the following units are adopted and used in the data stored in the `.csv` files:  
- Length, area, and volume: µm, µm², µm³  
- Time: s  
- Pressure and stress: GPa  
- Force: mN  
- Mass: 10³ kg = ton  
- Energy: nJ

The `.odb` files generated during the simulations follow a sequential naming scheme  
**Job-1, Job-2, …, Job-9**, corresponding to the following elementary load cases:

1. **Job-1** → x axial load  
2. **Job-2** → y axial load  
3. **Job-3** → z axial load  
4. **Job-4** → yz shear load  
5. **Job-5** → xz shear load  
6. **Job-6** → xy shear load  
7. **Job-7** → x + y biaxial load  
8. **Job-8** → x + z biaxial load  
9. **Job-9** → y + z biaxial load  

These loads are applied in the global coordinate system (**x** → radial, **y** → tangential, **z** → longitudinal). Since the `.odb` filenames are identical for each tissue and for the growth ring, the files are **overwritten every time a script is run**. Therefore, after the full workflow, the only remaining `.odb` files in the folder will be those of the **growth ring** simulations.
