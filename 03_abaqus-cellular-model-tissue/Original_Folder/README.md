# Simulation of viscoelastic creep of disordered tissues and growth ring

This folder contains all scripts and files used for the **simulation of viscoelastic creep** of the **disordered tissues** (EW, TW, LW), as well as the **growth ring (GR)** composed of them. All simulations are conducted with non-linear geometry activated (**nlgeom = ON**). The elastic properties must be known (calculated in `Elastic_Test`).

The global coordinate system follows the anatomical directions of wood tissues:  
**x** → tangential, **y** → radial, **z** → longitudinal for the tissues
**x** → radial, **y** → tangential, **z** → longitudinal for the growth ring.

## Folder structure and workflow overview

- **`run_for_MatSuMoTo.py`**  
  Main script to be **run from the terminal**. It sequentially launches:

  1. **`run_tissue_simulations.py`**  
     Simulates the viscoelastic creep of the individual Y-RVE tissues (EW, TW, LW).  
     It calls the corresponding scripts for each tissue:

     - **`Calculate_Tissue_Equivalent_Compliance_EW.py`**
     - **`Calculate_Tissue_Equivalent_Compliance_TW.py`**
     - **`Calculate_Tissue_Equivalent_Compliance_LW.py`**

     Each script performs Abaqus simulations under **nine elementary load-control cases**  
     (3 axial, 3 shear, 3 biaxial). From these, the equivalent creep compliance tensor is computed and stored in `.csv` files.

  2. **`fit_prony_tissue.py`**  
     Fits the tissue creep compliances using a Kelvin–Voigt (Prony series) model with a fixed number of elements (4) and fixed characteristic times (0.1, 1, 10, 100). Then, the **gamma ratios** are computed as the ratio between the elastic compliance and the compliance of each KV element for each direction; these are then averaged (as the values are nearly identical across directions). The resulting gamma values are stored in a `.csv` file (new values are appended; the last row is always used) and written into the include file **`tissues_gamma_values.inc`** for the growth ring Fortran material model.

  3. **`run_gr-rve_simulations.py`**  
     Simulates the viscoelastic creep of the growth ring based on the previously calculated Y-RVE creep properties. It launches **`Calculate_GR-RVE_Equivalent_Compliance.py`**, which performs the Abaqus simulations and stores the resulting creep compliance in `.csv` format.

  4. **`fit_prony_gr-rve.py`**  
     Fits the GR-RVE creep compliances using the same Kelvin–Voigt model used for the tissues. Gamma ratios are computed and stored in a `.csv` file (new values are appended).

  All simulations are performed at **constant moisture content**. Therefore, for each script `Calculate_...py`, the moisture content must be set as user input, and **the same value must be used for all tissues and for the GR-RVE** (currently 0.12). Moreover, the creep time must be set (currently 150).

- **`Fit_tissue_creep_compliances_from_Abaqus.ipynb`**  
  Jupyter notebook for inspecting creep compliance results, recomputing the Prony fits, and generating plots.

- **`Modules`**  
  Contains all classes, Python libraries, Fortran source files, and include files required to run the tissues and GR-RVE simulations:

  - **Structure_mat**
  Contains:  
  - all `CAE models` (and corresponding `.jnl` files) **generated from the `.mat` files** also included in this folder, and 
  - the `.mat` files defining the coordinates of the points for each fiber and the outer contour of the matrix (previously **produced via MATLAB**).

  - **`growth_ring_data.csv`**  
  File containing the geometric data required to generate the Y-RVE models and the growth ring.

  - **`Tissue_Class.py`**  
  Python class defining the 3D tissue geometry (smeared cell wall model), boundary conditions, loading cases, and data extraction routines.

  - **`GR_RVE_Class.py`**  
  Python class defining the 3D growth-ring geometry, boundary conditions, loading cases, and data extraction routines.

  - **`PBCfunction.py`**  
  Library implementing periodic boundary condition (PBC) functions.

  - **`TissuesMatModel.f`**  
  Fortran material model implementing the constitutive behaviour of the EW, TW, and LW tissues for the growth ring model.

  - **`TracheidsMatModel.f`**  
  Fortran material model implementing the constitutive behaviour of the smeared cell wall model used for the disordered tissues.

  - **`tracheids_eng_const.inc`**  
  Include file storing the equivalent engineering constants of the CML (matrix) and S-layer (fiber) for the tissue simulations.

  - **`tissues_eng_const.inc`**  
  Include file storing the equivalent engineering constants of the tissues for the growth ring simulations.

  - **`activate_ve_creep.inc`**  
  Include file controlling the activation of viscoelastic creep behaviour in the Fortran material models.

  - **`tracheids_gamma_values.inc`**  
  Include file storing the gamma values for the CML (matrix) and S-layer (fiber) (same set for all materials) used in the tissue simulations.

  - **`tissues_gamma_values.inc`**  
  Include file storing the gamma values for the tissues used in the growth ring simulations.

  - **`tracheids_gamma_initial_values.csv`**  
  Stores initial gamma values of the matrix and fiber. Used when the surrogate model is running to optimize these values.

## Note:
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
