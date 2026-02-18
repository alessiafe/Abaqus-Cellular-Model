# Viscoelastic creep simulations of hierarchical model

This folder contains the scripts and files required to compute and optimize the **viscoelastic creep** behavior of **wood hierarchical model** spanning three scales, from the material inputs at the cell wall scale to representative tissues (EW, TW, LW) that contribute to the growth ring model ([Ferrara and Wittel, 2025](https://doi.org/10.1007/s00707-025-04568-3)).

The workflow proceeds sequentially across the hierarchy, calculating the elastic compliance tensor of the cell wall layer, then the one of the cell wall, the one of the different tissue types and finally with those tensors the one of the growth ring RVE that result in the calculated macroscopic response. The same procedure is applied to the **viscoelastic creep**. Viscoelastic behavior at the cell-wall layer scale is represented using a Kelvin–Voigt (Prony) series with 4 elements and fixed characteristic times (0.1, 1, 10, 100). A set of **gamma** factors is introduced and each of them is defined as the ratio between the elastic compliance and the corresponding KV-element compliance. By computing the RVEs in sequence (cell wall → tissue → growth ring), the full macroscopic creep response is obtained.

At the cell-wall scale, elastic parameters are taken from the literature, and model topologies are fixed according to microscopic pictures of spruce tissues. The unknown gamma values at the cell-wall scale are identified through an **inverse surrogate-based optimization algorithm** that compares the creep response of the growth ring model with the experimental macroscopic measurements.

All simulations are conducted with non-linear geometry activated (**nlgeom = ON**).

## Folder structure and workflow overview

- **`Elastic_Test`**  
  Contains the scripts for computing the **elastic properties** (compliances and engineering constants) of the Y-RVE tissues and the growth ring (see dedicated README inside).

-  **`Original_Folder`**  
  Contains the scripts for computing the **viscoelastic creep** (KV-element compliances and gamma ratios) of the Y-RVE tissues and the growth ring. It relies on the elastic properties pre-computed in **`Elastic_Test`** and can be executed independently of the surrogate optimization routine (see dedicated README inside).
  
- **`MATSuMoTo-master`**  
  Contains the MATSuMoTo package for running inverse surrogate-based optimization.

- **`cellular_model_main.m`**  
  Main MATLAB script that runs the MATSuMoTo surrogate-based optimization algorithm. It relies on the following functions:

  - **`cellular_model_datainput.m`**  
  Defines the inputs for the surrogate optimization (number of unknown parameters, parameter bounds, objective function).

  - **`cellular_model_run_simulations.m`**  
  Copies the contents of **`Original_Folder`** into a new attempt folder and launches the corresponding simulation round via `run_for_MatSuMoTo.py`. It embeds the objective function, which compares the simulated creep response of the growth-ring model with the experimental reference data stored in **`macro_master_gamma.csv`**.

  The results of the algorithm are automatically stored in `Results.mat`, while the optimized **gamma** ratios (one set for each tissue and one for the growth ring) are saved in `best_gamma.csv`.

- **`macro_master_gamma.csv`**  
  Contains the gamma values for each component of the compliance tensor derived from the experimental macroscopic measurements.

- **`Plot_MatSuMoTo_vs_master_curves.ipynb`**  
  Jupyter notebook used to plot and compare the results of the surrogate algorithm with the experimentally derived **master curves**.
