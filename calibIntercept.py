"""
Python script that calibrates 
rainfall canopy interception parameters 
of Rutter model (1975)
------------------------
Author: Vaclav Steinbach
Date: 03.07.2025
Dissertation Work
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution 
import os
from rutter_intercept_calib import rutterIntercept

data_dir = 'data/'
out_dir = "out/"
os.makedirs(out_dir, exist_ok=True)

# load precipitation
R_free = np.loadtxt(data_dir+'rain_free.in')
R_tree = np.loadtxt(data_dir+'rain_tree.in')
# load air temperature
T = np.loadtxt(data_dir+'temp.in')
# load wind speed
u = np.loadtxt(data_dir+'wind.in')
# load relative humidity
rh = np.loadtxt(data_dir+'rh.in')
# load solar radiation
Rn = np.loadtxt(data_dir+'solar.in')

# define time and time step
time = R_free[:,0]
delta_t = time[1] - time[0]

def optIntercept(par):
    S = par[0]
    b = par[1]
    a = par[2]
    calib_vars = [S, b, a]
    # Compute Interception
    C,Tt = rutterIntercept(time,calib_vars,R_free,T,u,rh,Rn)
    # error in root mean square manner
    diff = Tt - R_tree[:,1]
    err = np.sum(diff**2)/len(diff)
    print(f"ERROR: {err}")
    return err

S_bnd = (0.1e-3,5e-3)
b_bnd = (-10000, 10000) # bounds for thermal convection parameter
a_bnd = (-50.0, 50.0)  # Bounds for thermal conductivity parameters
bounds = [S_bnd, b_bnd, a_bnd] 
result = differential_evolution(optIntercept, bounds, 
                                workers=-1, 
                                updating='deferred',
                                tol=1e-8,
                                atol = 1e-8,
                                maxiter=100000)

# Output the optimized parameter values and error
print("Optimized values:\n", result.x, '\n', result.fun)
