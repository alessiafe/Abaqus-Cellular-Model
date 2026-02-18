# Viscoelastic creep simulations of a hierarchical model

This folder contains the scripts and files required to compute and optimize the **viscoelastic creep** behavior of a **hierarchical wood model** spanning three scales, from material inputs at the cell wall scale to representative tissues (EW, TW, LW) contributing to the growth ring model ([Ferrara and Wittel, 2025](https://doi.org/10.1007/s00707-025-04568-3)). The elastic and viscoelastic properties are sequentially calculated across the hierarchy.

Three alternative versions of the model are available in the following folders, each accompanied by a dedicated README file:

- **`01_abaqus-cellular-model-Y-RVE-nlgeomOFF`**  
  This model is based on the discrete layered cell wall approach. Material properties are defined using composite mixing rules at the scale of the individual cell wall layers (ML, P, S1, S2, S3). The cellular structure is represented as a hexagonal grid, and its behavior is evaluated on an irreducible Y-shaped representative volume element (Y-RVE), formed by segments of three neighboring multi-layered tracheids. A growth ring RVE is then defined and discretized into three different materials: EW, TW, and LW. All tissue simulations are performed with non-linear geometry deactivated (**nlgeom = OFF**).

- **`02_abaqus-cellular-model-Y-RVE-nlgeomON`**  
  This model is identical to **`01_abaqus-cellular-model-Y-RVE-nlgeomOFF`**, except that tissue simulations are conducted with non-linear geometry activated (**nlgeom = ON**).

- **`03_abaqus-cellular-model-tissue`**  
  This model is based on the smeared cell wall approach, which represents an extension of the discrete layered cell wall model. In this case, multiple layers are merged into a single equivalent layer (ESL) based on laminate theory. Specifically, the secondary wall layers (S1, S2, S3) are combined into one ESL, while the primary wall (P) and middle lamella (ML) form an isotropic compound middle lamella (CML) embedding all S-tracheids. The cellular structure is derived from microscopy images of RT cross-sections to build representative tissue models characterized by geometric irregularity and disorder. The cell walls are modeled as an ESL embedded in the CML. All tissue simulations are performed with non-linear geometry activated (**nlgeom = ON**).
