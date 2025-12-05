# SIR Model: Epidemic Curve Fitting

This project implements and calibrates a Susceptible–Infectious–Recovered (SIR) model to real influenza data. It highlights skills in mathematical modeling, numerical methods, nonlinear optimization, and data analysis using Python.  

## Summary

The script:

- Defines the SIR system of differential equations 
- Fits transmission (`β`) and recovery (`γ`) parameters using least-squares minimization
- Imports real influenza case data from Excel  
- Compares the fitted model to observed cases and analyzes residuals  

## Outputs

The script produces:

- Observed vs. fitted infection curves  
- Residual diagnostics  
- Infection and recovery trajectories from the calibrated SIR model