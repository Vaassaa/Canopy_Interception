"""
Code for computation of precipitation interception
based on Rutter model (1971)
Author: Vaclav Steinbach
Date: 13.08.2025
Disseration Work
"""
import numpy as np
import matplotlib.pyplot as plt
import os
# from rutter_intercept import rutterIntercept
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

# define time
time = R_free[:,0]

# Read the pass optimalized values
file_path = out_dir+"intercept_opt.txt"
calib_vars = np.loadtxt(file_path, comments="#", max_rows=1)

# calib_vars = [5.20578360e-04,
              # 3.69054261e+03,
              # 5.76918252e+00]
# calib_vars = [5.42257250e-04,
              # 8.64141567e+03,
              # 2.60761163e+00]

# Compute interception
C,Tt = rutterIntercept(time,calib_vars,R_free,T,u,rh,Rn)
print(f"Computed interception with these parameters (S,b,a): \n {calib_vars}")

# Plot results
plt.plot(time, R_free[:,1], color = 'cornflowerblue', label="Rainfall")
plt.plot(time, R_tree[:,1], color = 'mediumblue', linestyle = '--', label="Measurement")
plt.plot(time, Tt, color = 'firebrick', label="Throughfall")
plt.plot(time, C, color = 'darksalmon', linestyle = '-.', label="Water on Canopy")
plt.grid()
plt.legend()
plt.xlabel("time [s]")
plt.ylabel("precipitation [m]")
plt.savefig(out_dir+'intercept.png')
plt.show()
