This folder contains R scripts and data files needed to run the summer and winter multi-species occupancy models and goodness-of-fit tests described in Griggs et al. (20XX).

MSOM scripts are adapted from code from Dyck et al. (2022; https://github.com/marissadyck/Multi_occ_LWW/tree/main). Both MSOM scripts fit single-season, multi-species occupancy models [following the approach from Rota et al. (2016)] using the unmarked package in R. Goodness of fit is assessed using parametric bootstrap tests with SSE, Pearson chi-square, and Freeman-Tukey discrepancy statistics following the approach used in Andrade-Ponce et al. (2022; https://zslpublications.onlinelibrary.wiley.com/doi/full/10.1111/jzo.13002).

------------------------------------------------------------
General Information
------------------------------------------------------------

Author Information: 
Authors: Johanna Griggs, Chris Darimont (University of Victoria), Cole Burton (University of British Columbia), Persia Khan (University of Victoria), Tyler Jessen (University of Victoria), and William Housty (Heiltsuk Integrated Resource Management Department)
Institution: Applied Conservation Science Lab, University of Victoria
Address: Victoria, BC, Canada
Email: johannagriggs@mac.com

Date of data collection: 2023-2025 

Geographic location of data collection: Bella Bella, BC, Canada and surrounding areas 

Information about funding sources that supported the collection of the data: This work was conducted in partnership with the Heiltsuk Integrated Resource Management Department (HIRMD). Funding for Johanna Griggs and the Heiltsuk Camera Project was provided by the University of Victoria, Raincoast Conservation Foundation, the UVic Edward Bassett Family Bursary, the Together for Wildlife Scholarship (managed by the Habitat Conservation Trust Foundation), and the Natural Science and Engineering Council of Canada Graduate Research Scholarship (Master’s).

------------------------------------------------------------
Folder structure and description of contents 
------------------------------------------------------------
Please keep the folder structure exactly as shown below. 
The scripts use relative paths based on the main project folder; the working directory should be set to the main folder, Griggs_et_al_MSOMs.


Griggs_et_al_MSOMs/
├── Griggs_et_al_MSOMs_ReadMe
├── scripts/
│   ├── covariates_correlation_matrices.R
│   ├── MSOM_summer_biweekly.R
│   └── MSOM_winter_biweekly.R
└── inputs (all files below are in same folder)
│   └── om_predict_fix.R (source code provided from Ken Kellner via Dyck et al. to fix package bugs related to the predict function)
│   └── Summer detection history dataframes in .csv file format (deer_summer, wolf_summer)
│   └── Winter detection history dataframes in .csv file format (deer_winter, wolf_winter)
│   └── Detection covariate dataframes in .csv file format used for both seasons (det_cov_trail_strength, det_cov_year)
│   └── Summer-specific detection (observation) covariate dataframes in .csv file format
│		 └── summer_biweekly_cam_effort
│		 └── det_cov_deer_summer_binary
│		 └── det_cov_deer_summer_prop_det_biweekly
│		 └── det_cov_wolf_summer_binary
│		 └── det_cov_wolf_summer_prop_det_biweekly
│	└── Winter-specific detection (observation) covariate dataframes in .csv file format
│		 └── winter_biweekly_cam_effort
│		 └── det_cov_deer_winter_binary
│		 └── det_cov_deer_winter_prop_det_biweekly
│		 └── det_cov_wolf_winter_binary
│		 └── det_cov_wolf_winter_prop_det_biweekly
│   └── Occupancy (site) covariate dataframes in .csv file format used for both seasons
│		 └── occu_site_covs_biweekly_MSOM (standardized)
│		 └── orig_site_covs_biweekly_MSOM (raw, for back-transforming plot axes)
│   └── List of all site-year combinations (site_year_list_biweekly)
├── outputs
│   ├── summer_mod/
│       └── [CSV files produced from MSOM_summer_biweekly.R script]
│       └── [GOF test results (rds. files) produced from MSOM_summer_biweekly.R script]
│   └── winter_mod/
│       └── [CSV files produced from MSOM_winter_biweekly.R script]
│       └── [GOF test results (rds. files) produced from MSOM_winter_biweekly.R script]

------------------------------------------------------------
How to run the code
------------------------------------------------------------

Open the project folder in RStudio and make sure the working directory is the main folder: Griggs_et_al_MSOMs

Then you can run the script for either season:

scripts/MSOM_summer_biweekly.R
or:
scripts/MSOM_winter_biweekly.R

You can also run the covariates_correlation_matrices.R script to: 
- check for collinearity between occupancy and detection covariates.
- produce correlation matrices and GVIF table

------------------------------------------------------------
Notes on Goodness-of-Fit
------------------------------------------------------------

Goodness-of-fit is assessed twice, once for the global model to estimate c-hat which is subsequently used in model selection and then again for the top performing model.
It should be noted that the parametric bootstrap goodness-of-fit tests can take a long time to run because they use 1000 simulations.

The saved parboot objects are included in the outputs/ folder, the scripts can load those results directly instead of rerunning the bootstrap.

The Pearson chi-square statistic is included because it is commonly used for occupancy-model goodness-of-fit testing. However, for sparse binary detection data, the Pearson statistic can be sensitive to very small fitted expected values because each contribution is divided by the expected value. For this reason, the scripts also calculate SSE and Freeman-Tukey discrepancy statistics. These additional diagnostics are used to evaluate whether any Pearson chi-square results are being strongly influenced by rare, extreme bootstrap values.

Note on Winter MSOM GOF results: 
The winter global model produced a highly right-skewed Pearson chi-square bootstrap distribution because a small number of simulated and refitted datasets assigned near-zero expected probabilities to detected wolves, resulting in extremely large cell-level Pearson contributions. These rare values disproportionately inflated the mean simulated Pearson statistic, making the corresponding Pearson-based overdispersion estimate (c-hat) uninformative. SSE and Freeman–Tukey bootstrap distributions remained stable and did not show comparable extreme values. Winter overdispersion was therefore estimated using the SSE statistic (c-hat = 1.12), while the Pearson and Freeman–Tukey results were retained as supplementary goodness-of-fit diagnostics.