# SIR Model: Epidemic Curve Fitting

This project implements and calibrates a Susceptible–Infectious–Recovered (SIR) model to real influenza data. It highlights skills in mathematical modeling, numerical methods, nonlinear optimization, and data analysis using Python.  

## Summary

The script:

- Defines the SIR system of differential equations 
- Estimates the transmission and recovery parameters using nonlinear least-squares.
- Loads real weekly case counts from Excel. 
- Compares the model’s fitted trajectory to observed data.
- Produces residual diagnostics and system-evolution plots.

## Outputs

This repository includes:

- Observed vs. fitted infection curve 
- Residual plot showing model error structure
- Evolution of the infected and recovered populations under the fitted parameters