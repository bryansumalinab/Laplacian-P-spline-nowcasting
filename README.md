# Bayesian nowcasting with Laplacian-P-splines
Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes

## About this repository
This repository contains R codes used to generate results from the paper "Bayesian nowcasting with Laplacian-P-splines" by Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes.

## R codes for simulation:
The R codes for performing the simulation study are contained in the **LPS** and **VanDeKassteele** folders, corresponding to the Laplacian-P-splines (LPS) and Van de Kassteele et al. (2019) method, respectively.

## R codes for real data application:
The script **mortality_nowcasting.R** is used for nowcasting of mortality data for one nowcast date. The script **mortality_diffnd.R** contains R codes that loopsthrough all the different nowcast dates considered in the analysis of mortality data. The codes for plotting the nowcast and (estimated) delay density are contained in these scripts.

## Source functions
There are two source functions for the LPS method:
1. **Nowcasting.R** - this include the day of the week effects in the model that is used for the analysis of 2021 Belgian COVID-19 mortality data (**mort2021.xlsx**).
2. **Nowcasting_sim.R** - source function to perform the simulation study (without day of the week effects) using LPS method.
