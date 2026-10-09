# High-Dimensional Matrix Estimation through Elliptical Factor Models

This repository contains code and data for the paper **“High-Dimensional Matrix Estimation through Elliptical Factor Models”**.

## Files

- `code/factor_number_simulation.R`: factor number estimation
- `code/scatter_matrix_vH.R`: simulation code for scatter matrix estimation
- `code/POET_SS_H0_vH.R`: simulation code for POET_SS estimator and its target
- `code/covariance_matrix_vH.R`: simulation code for covariance matrix estimation
- `code/precision_matrix_vH.R`: simulation code for precision matrix estimation
- `code/SP500_vH.R`: empirical analysis code for S&P 500 data
- `sp500_anymember_excess_return_wide.rar`: S&P 500 return data
- `sp500_excess_return_wide.rar`: S&P 500 return data

## Requirements

The code is written in R. Required packages include:

`matrixStats`, `parallel`, `MASS`, `SPCAvRP`, `mvtnorm`, `mnormt`, `ICSNP`, `ICS`, `SpatialNP`, `MNM`, `flare`, `glassoFast`

## Usage

Run the scripts in the `code/` folder for simulations and empirical analysis.  
Please extract the `.rar` data files before running the S&P 500 analysis.

## License

MIT License.
