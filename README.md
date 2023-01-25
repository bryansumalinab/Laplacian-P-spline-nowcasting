# Bayesian nowcasting with Laplacian-P-splines
Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes

## About this repository
This repository contains R codes used to generate results from the paper "Bayesian nowcasting with Laplacian-P-splines" by Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes.

## R codes for simulation:
The R codes used in the simulation study can be found in the **Simulation** folder. There are two subfolders for the simulation codes, the **LPS** and **VanDeKassteele** folders:
1. **LPS** - This folder contains the codes for the simulation study using our proposed Laplacian-P-splines (LPS) method. Inside this folder are the three R scripts **LPS_f11.R**, **LPS_f12.R**, and **LPS_f2.R**, which generate results using the three different functions considered in the simulation, namely: functions $f_{11}(t)$, $f_{12}(t)$, and $f_2(t)$, respectively. To run these codes, the script **Nowcasting_sim.R** needs to be loaded, which is also contained in the LPS folder.
2. **VanDeKassteele** - This folder contains R codes for the simulation using van de Kassteele et al. (2019) method. Similar to the LPS folder, it contains three R scripts: **van_f11.R**, **van_f12.R**, and **van_f2.R**, corresponding to functions $f_{11}(t)$, $f_{12}(t)$, and $f_2(t)$, respectively. Inside the vanDeKassteele folder is the folder named **van_functions.R** containing R codes that need to be loaded before running the three simulation codes.

## R codes for real data application:
The script **mortality_nowcasting.R** is used for nowcasting of mortality data for one nowcast date. The script **mortality_diffnd.R** contains R codes that loops through all the different nowcast dates considered in the analysis of mortality data. The codes for plotting the nowcast and (estimated) delay density are contained in these scripts.

## Source functions
There are two source functions for the LPS method:
1. **Nowcasting.R** - this include the day of the week effects in the model that is used for the analysis of 2021 Belgian COVID-19 mortality data (**mort2021.xlsx**).
2. **Nowcasting_sim.R** - source function to perform the simulation study (without day of the week effects) using LPS method.
