# CRM
This repository contains the R package `crmReg` and the corresponding articles:

Debruyne, M., Höppner, S., Serneels, S., and Verdonck, T. (2019). Outlyingness: Which variables contribute most? Statistics and Computing, 29 (4), 707–723. DOI: 10.1007/s11222-018-9831-5  
https://link.springer.com/article/10.1007/s11222-018-9831-5  
Filzmoser, P., Höppner, S., Ortner, I., Serneels, S., and Verdonck, T. (2020). Cellwise Robust M Regression. Computational Statistics and Data Analysis, 147:106944. DOI: 10.1016/j.csda.2020.106944  
https://doi.org/10.1016/j.csda.2020.106944

The R package `crmReg` contains the implementation of the Cellwise Robust M-regression (CRM) algorithm and the SPArse DIrections of Maximal Outlyingness (SPADIMO) algorithm. The `crmReg` package also includes the predict function for fitted CRM regression models, the data preprocessing function used by CRM as well as the function that creates the cellwise heatmaps. Furthermore, this repository contains the R scripts, as used in Filzmoser et al. (2020), for the simuation studies and the application of CRM on the Nutrients data set.

The R package can be installed from CRAN: https://CRAN.R-project.org/package=crmReg  
or through the `devtools` package: `devtools::install_github("SebastiaanHoppner/CRM/crmReg")`

