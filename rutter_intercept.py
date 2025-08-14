"""
Computer implementation of Rutter canopy interception model
Author: Vaclav Steinbach
Date: 13.08.2025
Dissertation work
"""
import numpy as np

def rutterIntercept(time,R,T,u,rh,Rn):
    # --- Canopy Parameters ---
    p = 0.25         # free throughfall [-]
    S = 1.4e-3          # canopy storage capacity [m] (typically 0.5–2 mm)
    b = 3.7          # first drainage exponential coefficient [?]
    a = -18          # second drainage exponential coefficient [?]
    z = 10           # reference height [m]
    h = 12           # tree heigh [m]
    d = 0.75*h       # displacement height [m]
    z_0 = 0.1*h      # roughness length [m]
    ra = 1/(0.4**2 * np.mean(u)) \
         * (np.log((z - d) / z_0))**2 # aerodynamic resistance [s/m]
    # ra = 9.2 / np.mean(u)    
    # a = np.log(0.002)-b*S         # drainage coefficient a in D=exp(a+bC)

    # --- Physical constants ---
    cp = 1006        # specific heat of air [J/kg/°C]
    rho = 1.204      # air density [kg/m³]
    lambda_v = 2.45e6 # latent heat [J/kg]
    P = 101325      # atmospheric pressure [Pa]
    MW = 0.622      # ratio molecular weight of water vapor/dry air [-]
    gamma = (cp * P) / (lambda_v * MW)  # psychrometric constant [kPa/°C]
    gamma = gamma / 1000 # kPa -> Pa

    # time step
    delta_t = time[1] - time[0]
    # Allocation
    C = 0.0
    D = 0.0
    Tt = np.zeros_like(R)
    for t in range(len(R)):
        print(f'TIME STEP: {time[t]}')
        
        # Compute potential evaporation at this time step
        # Slope of saturation vapor pressure curve, Δ [kPa/°C]
        es = 0.6108 * np.exp(17.27 * T[t,1] / (T[t,1] + 243.04))
        es = es / 1000 # kPa -> Pa
        delta = 4098 * es / ((T[t,1] + 237.3) ** 2)

        ea = es * rh[t,1]
        VPD = es - ea
        Ep = (Rn[t,1] + rho*cp*VPD/ra) / (lambda_v*(delta + gamma))
        # Ep = Ep/3600 # kgm^-2s^-1 -> mm/h
        Ep = Ep/1000 # mm/s -> m/s
        
        # Compute evaporation from wet surface
        if C > S:
            E = Ep
        else:
            E = Ep * (C / S)

        print(f'evap: {E}')
        print(f'rain: {R[t,1]}')

        # Compute drainage from canopy
        if C > 0:
            D = np.exp(a + b * C)
        else:
            D = 0.0
        D = D/(1000*3600) # mm/h -> m/s

        print(f'drain: {D}')

        # Water balance on canopy
        dC = ((1 - p) * R[t,1] - D - E) * delta_t
        print(f'dC: {dC}')

        # Update canopy storage with physical bounds
        C = max(0.0, min(S, C + dC))
        print(f'C: {C}')

        # Throughfall reaching the soil
        Tt[t] = p * R[t,1] + D

    return Tt
