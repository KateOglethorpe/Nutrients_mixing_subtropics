# **Mesoscale Eddy Stirring & Microscale Turbulence in Subtropical Gyres**

This repository contains **output data** and **MATLAB scripts** from the paper: **"The Role of Mesoscale Eddy Stirring and Microscale Turbulence in Sustaining Biological Production in Subtropical Gyres."**
Co-authors: K. Oglethorpe (ko389@cam.ac.uk), Dr C. Spingys, Dr B. Fernández-Castro, Prof. A. C. Naveira Garabato, Prof. R. G. Williams. 
Paper submitted to *Global Biogeochemical Cycles*.
DOI: 10.5281/zenodo.14771059

## **Description**
The six output files (.mat) contain estimates of **nitrate supply to the winter mixed layer** (`WML`) **and thermocline** by:
- **Mesoscale eddy stirring**, derived from the convergences of isopycnal diffusive nitrate fluxes (`F_iso_diff`).
- **Microscale turbulence**, derived from the sum of convergences of diapycnal diffusive nitrate fluxes (`F_dia_diff`) and diapycnal advective nitrate fluxes (`F_dia_adv`).

Units: mol N m⁻² yr⁻¹

| File Name               | Description |
|-------------------------|-------------|
| `wml_F_iso_diff.mat`    | Isopycnal diffusive nitrate supply to the WML |
| `wml_F_dia_diff.mat`    | Diapycnal diffusive nitrate supply to the WML |
| `wml_F_dia_adv.mat`     | Diapycnal advective nitrate supply to the WML |
| `therm_F_iso_diff.mat`  | Isopycnal diffusive nitrate supply to the Thermocline |
| `therm_F_dia_diff.mat`  | Diapycnal diffusive nitrate supply to the Thermocline |
| `therm_F_dia_adv.mat`   | Diapycnal advective nitrate supply to the Thermocline |

The three MATLAB scripts (.m) files contain code to calculate these rates of nitrate supply from global climatologies and estimates (`Calculations.m`) and plot the paper figures (`North_Atlantic_figures.m` and `Global_figures.m`).

## **Data Sources**
This research utilizes publicly available datasets:

- **World Ocean Atlas 2018 (WOA2018) Data**  
  Access via the **National Oceanic and Atmospheric Administration (NOAA)**:  
  [World Ocean Atlas 2018](https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/)  
  Garcia et al. (2018), Locarnini et al. (2018), Zweng et al. (2019)

- **Global Estimates of Isopycnal Diffusivity**  
  Groeskamp et al. (2020), available via **Figshare**:  
  [DOI: 10.6084/m9.figshare.12554555.v2](https://doi.org/10.6084/m9.figshare.12554555.v2)

- **Gibbs Seawater (GSW) Oceanographic Toolbox**  
  Available at **TEOS-10 website**:  
  [GSW Toolbox](https://www.teos-10.org/software.htm)  
  Jackett & McDougall (1997)
