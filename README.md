# Bayesian nowcasting with Laplacian-P-splines
Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes

## About this repository
This repository contains R codes used to generate results from the paper "Bayesian nowcasting with Laplacian-P-splines" by Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes.

## R codes for simulation:
The R codes used in the simulation study can be found in the **Simulation** folder. There are three subfolders for the simulation codes, the (1) **LPSNB** (2) **LPSPoisson** and (3) **VanDeKassteele** folders:
1. **LPSNB** - This folder contains the codes for the simulation study using our proposed Laplacian-P-splines (LPS) method assuming a *negative binomial distribution* on the number of cases. Inside this folder are the four R scripts **LPSNB_f11.R**, **LPSNB_f12.R**, and **LPSNB_f21.R**, **LPSNB_f22.R**, which generate results using the four different functions considered in the simulation, namely: functions $f_{11}(t)$, $f_{12}(t)$, $f_{21}(t)$, and $f_{22}(t)$, respectively. To run these codes, the scripts in the folder **LPS_functions** need to be loaded, which is also contained in the LPSNB folder.
2. **VanDeKassteele** - This folder contains R codes for the simulation using van de Kassteele et al. (2019) method. Similar to the LPS folder, it contains four R scripts: **van_f11.R**, **van_f12.R**, **van_f21.R** and **van_f22.R**, corresponding to functions $f_{11}(t)$, $f_{12}(t)$, $f_{21}(t)$ and $f_{22}(t)$, respectively. Inside the vanDeKassteele folder is the folder named **van_functions.R** containing functions that need to be loaded before running the three simulation codes.
3. **LPSPoisson** - This folder contains R codes for the LPS method, assuming a *Poisson* distribution on the number of cases.

Each of the above codes will produce simulation results such as prediction interval coverage, mean absolute percentage error (MAPE), and symmetric mean absolute percentage error (SMAPE) for the specified number of iterations. The results are only for a specified nowcast date, and the scripts must be re-run to obtain results for the other nowcast dates. This can be set at the beginning of the R script with the object name *date.now*.

## Data
The data used in the paper is the 2021/2022 COVID-19 mortality/incidence data in Belgium which can be found in the **Data** folder named **mort2021.xlsx** and **incidence2022.xlsx**.

## R codes for real data application:
The script used for the analysis of mortality data is contained in the **Real Data Application** folder.
