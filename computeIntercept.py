"""
Code for computation of precipitation interception
based on Rutter model (1971)
Author: Vaclav Steinbach
Date: 13.08.2025
Disseration Work
"""
import numpy as np
import matplotlib.pyplot as plt
from rutter_intercept import rutterIntercept

data_dir = 'data/'
# load precipitation
R = np.loadtxt(data_dir+'rain.in')
# load air temperature
T = np.loadtxt(data_dir+'temp.in')
# load wind speed
u = np.loadtxt(data_dir+'wind.in')
# load relative humidity
rh = np.loadtxt(data_dir+'rh.in')
# load solar radiation
Rn = np.loadtxt(data_dir+'solar.in')

# define time
time = R[:,0]

# Compute Interception
Tt = rutterIntercept(time,R,T,u,rh,Rn)

# plot results
plt.plot(time, R[:,1])
plt.plot(time, Tt)
plt.show()
