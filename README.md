# Bayesian nowcasting with Laplacian-P-splines
Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes

## About this repository
This repository contains R codes used to generate results from the paper "Bayesian nowcasting with Laplacian-P-splines" by Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes.

## R codes for simulation:
The R codes used in the simulation study can be found in the **Simulation** folder. There are two subfolders for the simulation codes, the **LPS** and **VanDeKassteele** folders:
1. **LPS** - This folder contains the codes for the simulation study using our proposed Laplacian-P-splines (LPS) method. Inside this folder are the three R scripts **LPS_f11.R**, **LPS_f12.R**, and **LPS_f2.R**, which generate results using the three different functions considered in the simulation, namely: functions $f_{11}(t)$, $f_{12}(t)$, and $f_2(t)$, respectively. To run these codes, the script **Nowcasting_sim.R** needs to be loaded, which is also contained in the LPS folder.
2. **VanDeKassteele** - This folder contains R codes for the simulation using van de Kassteele et al. (2019) method. Similar to the LPS folder, it contains three R scripts: **van_f11.R**, **van_f12.R**, and **van_f2.R**, corresponding to functions $f_{11}(t)$, $f_{12}(t)$, and $f_2(t)$, respectively. Inside the vanDeKassteele folder is the folder named **van_functions.R** containing R codes that need to be loaded before running the three simulation codes.

Each of the above codes will produce simulation results such as prediction interval coverage, mean absolute percentage error (MAPE), and symmetric mean absolute percentage error (SMAPE) for the specified number of iterations (in our case 1000). The results are only for a specified nowcast date, and the scripts must be re-run to obtain results for the other nowcast dates. This can be set at the beginning of the R script with the object name *date.now*.

## Data
The data used in the paper is the 2021 COVID-19 mortality data in Belgium which can be found in the **Data** folder named **mort2021.xlsx**.

## R codes for real data application:
The R codes used for the analysis of mortality data are contained in the Real Data Application folder which contains three R scripts: **mortality_nowcasting.R**, **mortality_diffnd.R** and **Nowcasting.R**.

1. **mortality_nowcasting.R** - The script **mortality_nowcasting.R** is used for the nowcasting of mortality data for only one nowcast date. This code will produce the nowcast plot and the delay density plot for a specified nowcast date.
2. **mortality_diffnd.R** - The script **mortality_diffnd.R** on the other hand, contains R codes that loop through all the different nowcast dates considered in the analysis of mortality data. This will produce the nowcast plot (Figure 5 in the paper) and estimated delay density plot (Figure 6 in the paper) for eight nowcast dates.
3. **Nowcasting.R** - This source function needs to be loaded before running the codes **mortality_nowcasting.R** and **mortality_diffnd.R**.
