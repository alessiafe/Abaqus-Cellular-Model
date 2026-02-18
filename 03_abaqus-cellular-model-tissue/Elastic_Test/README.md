# Calculation of elastic properties of disordered tissues and growth ring

This folder contains all scripts and files used for the **calculation of the elastic properties (compliances and engineering constants)** of the **disordered tissues** (EW, TW, LW), as well as the **growth ring (GR)** composed of them. All simulations are conducted with non-linear geometry activated (**nlgeom = ON**).

The global coordinate system follows the anatomical directions of wood tissues:  
**x** → tangential, **y** → radial, **z** → longitudinal for the tissues
**x** → radial, **y** → tangential, **z** → longitudinal for the growth ring.

## Folder structure and workflow overview

- **`Set_S-layer_elastic_properties`**
  Contains all scripts and files used for the **micromaterial calculations of the smeared cell wall model** where the secondary wall layers (S1, S2, S3) are merged into one Equivalent Single Layer (ESL), while the primary wall (P) and middle lamella (ML) form an isotropic compound middle lamella (CML) that embeds all S-tracheids:

  - **Micromaterial_Calculation.py**  
  Library of functions to compute the hygroelastic properties of the wood cell wall layers, as well as the equivalent properties of planar and 3D laminate composites.
  
  - **basic_components.npy**  
  NumPy array storing the stiffness matrices and hygroexpansion coefficients of the basic wood cell wall components (cellulose, hemicellulose, and lignin).

  - **materials_composite.npy**  
  NumPy array storing data required to compute the composite materials  
  (rule, fiber, fiber volume fraction, matrix, matrix volume fraction, and geometry coefficients)  
  from the basic components (rule, engineering constants, hygroexpansion coefficients).  

  - **Fit_cell_wall_engineering_constants_for_tissues.ipynb**
   Fits the engineering constants of CML (matrix) and S-layer (fiber) over multiple moisture contents (if available) and prepares a Fortran inc file (**`tracheids_eng_const.inc`**) for UMAT implementation which is stored in the current folder and in the `Original_Folder` (for later creep simulations).

  - **`Fit_tissues_engineering_constants_for_GR.ipynb`**  
  Fits the engineering constants of tissues over multiple moisture contents (if available) and prepares a Fortran-formatted text file for UMAT implementation. *(This workflow is no longer used)*

- **`run_elastic_test.py`**  
  Main script to be **run from the terminal**. It sequentially launches:

  1. **`run_tissue_simulations.py`**  
     Computes the elastic properties of the individual disordered tissues (EW, TW, LW).  
     It calls the corresponding scripts for each tissue:

     - **`Calculate_Tissue_Equivalent_Compliance_EW.py`**
     - **`Calculate_Tissue_Equivalent_Compliance_TW.py`**
     - **`Calculate_Tissue_Equivalent_Compliance_LW.py`**

     Each script performs Abaqus simulations under **nine elementary load-control cases**  
     (3 axial, 3 shear, 3 biaxial). From these, the equivalent compliance tensor and engineering constants are computed. Results are stored in the current folder and in the `Original_Folder` (for later creep simulations).

  2. **`run_gr-rve_simulations.py`**  
     Computes the equivalent elastic properties of the growth ring based on the previously calculated tissue properties. It launches **`Calculate_GR-RVE_Equivalent_Compliance.py`**, which performs the Abaqus simulations and saves results as above.

  All simulations are performed at **constant moisture content**. Therefore, for each script `Calculate_...py`, the moisture content must be set as user input, and **the same value must be used for all tissues and for the GR-RVE** (currently 0.12).

- **`Collect_engineering_constants.ipynb`**  
  Collects and organizes the engineering constants [GPa] from all simulations into a summary table stored as **`collect_engineering_constants.csv`**.

- **`Modules`**  
  Contains all classes, Python libraries, Fortran source files, and include files required to run the tissue and GR-RVE simulations:

  - **Structure_mat**
  Contains:  
  - all `CAE models` (and corresponding `.jnl` files) **generated from the `.mat` files** also included in this folder, and 
  - the `.mat` files defining the coordinates of the points for each fiber and the outer contour of the matrix (previously **produced via MATLAB**).

  - **`growth_ring_data.csv`**  
  File containing geometric data required to generate the Y-RVE models and the growth ring.

  - **`Tissue_Class.py`**  
  Python class defining the 3D tissue geometry (smeared cell wall model), boundary conditions, loading cases, and data extraction routines.

  - **`GR_RVE_Class.py`**  
  Python class defining the 3D growth-ring geometry, boundary conditions, loading cases, and data extraction routines.

  - **`PBCfunction.py`**  
  Library containing the periodic boundary condition (PBC) functions.

  - **`TissuesMatModel.f`**  
  Fortran material model implementing the constitutive behaviour of the EW, TW, and LW tissues for the growth ring model.

  - **`TracheidsMatModel.f`**  
  Fortran material model implementing the constitutive behaviour of the smeared cell wall model used for the disordered tissues.

  - **`tracheids_eng_const.inc`**  
  Include file storing the equivalent engineering constants of the CML (matrix) and S-layer (fiber) for the tissue simulations.

  - **`tissues_eng_const.inc`**  
  Include file storing the equivalent engineering constants of the tissues for the growth ring simulations.

  - **`activate_ve_creep.inc`**  
  Include file controlling the activation of viscoelastic/creep behaviour in the Fortran material models. *Must remain disabled for elastic tests.*

  - **`tracheids_gamma_values.inc`**  
  Include file storing the gamma values for the CML (matrix) and S-layer (fiber) (same set for all materials) used in the tissue simulations. Not used in the elastic simulations, but required for the Fortran material routine to run (the same routine used for creep simulations).

  - **`tissues_gamma_values.inc`**  
  Include file storing the gamma values for the tissues used in the growth ring simulations. Not used in the elastic simulations, but required for the Fortran material routine to run (the same routine used for creep simulations).

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

These loads are applied in the global coordinate system. Since the `.odb` filenames are identical for each tissue and for the growth ring, the files are **overwritten every time a script is run**. Therefore, after the full workflow, the only remaining `.odb` files in the folder will be those of the **growth ring** simulations.
