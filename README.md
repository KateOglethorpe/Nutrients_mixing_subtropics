# **Mesoscale Eddy Stirring & Microscale Turbulence in Subtropical Gyres**

This repository contains **six output files** and **MATLAB code** from the research project:  
**"The Role of Mesoscale Eddy Stirring and Microscale Turbulence in Sustaining Biological Production in Subtropical Gyres."**
Co-authors: K. Oglethorpe (ko389@cam.ac.uk), Dr C. Spingys, Dr B. Fernández-Castro, Prof. A. C. Naveira Garabato, Prof. R. G. Williams.

## **Description**
The six output files (.mat) contain estimates of **nitrate supply to the Winter Mixed Layer (WML) and Thermocline (THERM)** by:
- **Mesoscale eddy stirring**, derived from the convergences of isopycnal diffusive nitrate fluxes (`F_iso_diff`).
- **Microscale turbulence**, derived from the sum of convergences of diapycnal diffusive nitrate fluxes (`F_dia_diff`) and diapycnal advective nitrate fluxes (`F_dia_adv`).

## **Files**
| File Name               | Description |
|-------------------------|-------------|
| `wml_F_iso_diff.mat`    | Isopycnal diffusive nitrate flux convergence in WML |
| `wml_F_dia_diff.mat`    | Diapycnal diffusive nitrate flux convergence in WML |
| `wml_F_dia_adv.mat`     | Diapycnal advective nitrate flux convergence in WML |
| `therm_F_iso_diff.mat`  | Isopycnal diffusive nitrate flux convergence in Thermocline |
| `therm_F_dia_diff.mat`  | Diapycnal diffusive nitrate flux convergence in Thermocline |
| `therm_F_dia_adv.mat`   | Diapycnal advective nitrate flux convergence in Thermocline |

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
