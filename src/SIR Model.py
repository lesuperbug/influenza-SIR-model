# %% Import Packages

import scipy.integrate as solver
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lmfit as lmf

# %% Define Model and Solver

def SIR_model(y0, tspan, population, beta, gamma):
    '''
    This function is an auxiliary function for the ODE solver - It defines what is to be solved.
    The 'tspan' argument is only necessary becuase scipy.integrate.odeint requires that to be 
    an argument. Note that 'tspan' is a required argument, but that it is not used to calculate
    anything in this function.
    '''
    
    # IVP
    S, I, R = y0
    N = population
    
    # Equations
    dS_dt = -beta*I*S/N
    dI_dt = beta*I*S/N - gamma*I
    dR_dt = gamma*I
    
    return [dS_dt, dI_dt, dR_dt]

# %% Import Data

inc_data = pd.read_excel("../data/2017 - 2018 Flu Data (Manual Extract).xlsx")
inc_data = inc_data.loc[:, ['Week', 'Cases']]

# %% Initial Conditions and Other Preliminaries

# Parameters
N = 4480486
beta = 1.5                   # Transmission Rate
gamma = .5                   # Recovery Rate
arg_space = (N, beta, gamma)
params = lmf.Parameters()
params.add('beta', value = beta, min = 0)
params.add('gamma', value = gamma, min = 0)   # Problems arise if you set vary = False

# Initial Values
I0 = 1                       # Initial Condition - Infectious
R0 = 0                       # Initial Condition - Recovered
S0 = N - I0 - R0             # Initial Condition - Susceptible
y0 = [S0, I0, R0]

# Time Domain / Interval of Problem
period = 52
tspan = np.linspace(0, period - 1, period)   # np.arange also performs this function well

# %% Non-Linear Parameter Estimation

def error(params, model, y0, population, tspan, data):
    '''
    This function simply returns the error in the model for the given parameters
    '''
    
    # Retrieve Parameters - they need to be passed as a tuple into .odeint
    beta, gamma = params['beta'].value, params['gamma'].value
    unpacked_params = (population, beta, gamma)
    
    # Defining and Returning the Residuals
    solution = solver.odeint(model, y0, tspan, args = unpacked_params)
    return (solution[:, 1] - data).ravel()   # Column 2 contains the modelled infections

result = lmf.minimize(error, params, args = (SIR_model, y0, N, tspan, inc_data.iloc[:, 1]), method = 'leastsq')

# %% Solving the ODEs

# Re-Specifying our Parameters
beta = result.params['beta'].value            
gamma = result.params['gamma'].value            
arg_space = (N, beta, gamma)
solution = solver.odeint(SIR_model, y0, tspan, args = arg_space)

# %% Plot our fit

plt.scatter(tspan, inc_data.iloc[:, 1], label = "Observed Cases", color = 'b')
plt.plot(tspan, solution[:, 1], label = "Fitted Cases", color = 'r')
plt.grid()
plt.legend()
plt.xlabel("Weeks")
plt.ylabel("Cases")
plt.title("Model Fit: Infected Individuals")
plt.savefig("../output/fit_curve.png", dpi=300, bbox_inches="tight")   
plt.show()

# %% Plotting our Residuals

zero_line = [0 for i in range(52)]
residuals = (solution[:, 1] - inc_data.iloc[:, 1])
plt.plot(solution[:, 1], zero_line, linewidth = 0.5, color = 'black')
plt.scatter(solution[:, 1], residuals, color = 'r')
plt.xlabel("Fitted Value")
plt.ylabel("Residuals")
plt.title("Residuals of Fitted Infection Values")
plt.savefig("../output/residuals_fit.png", dpi=300, bbox_inches="tight")   
plt.show()

# %% Plotting the Results

plt.plot(tspan, solution[:, 1], label = "I(t)")
plt.plot(tspan, solution[:, 2], label = "R(t)")
plt.grid()
plt.legend()
plt.xlabel("Weeks")
plt.ylabel("Count")
plt.title("# of Infected / Recovered")
plt.savefig("../output/sir_evolution.png", dpi=300, bbox_inches="tight")
plt.show()