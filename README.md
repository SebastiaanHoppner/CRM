# CRM
This repository contains the R package and the corresponding papers on Cellwise Robust M-regression (2019) by P. Filzmoser, S. Höppner, I. Ortner, S. Serneels and T. Verdonck (https://doi.org/10.1016/j.csda.2020.106944) and Outlyingness: Which variables contribute most? (2019) (https://link.springer.com/article/10.1007/s11222-018-9831-5).

The R package `crmReg` contains the implementation of the Cellwise Robust M-regression (CRM) algorithm and the SPArse DIrections of Maximal Outlyingness (SPADIMO) algorithm. The `crmReg` package also includes the predict function for fitted CRM regression models, the data preprocessing function used by CRM as well as the function that creates the cellwise heatmaps in the paper. Furthermore, this repository contains the R scripts, as used in the paper, for the simuation studies and the application of CRM on the Nutrients data set.

The R package can be installed from CRAN: https://CRAN.R-project.org/package=crmReg  
or through the `devtools` package: `devtools::install_github("SebastiaanHoppner/CRM/crmReg")`

