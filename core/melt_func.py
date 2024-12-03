# Pressure dependant melting temperature models
import numpy as np

# Alfe et al. 2002 50 - 350 GPa

# Sixtrude Silicate liquidus curve 2014
def T_ms2014(P, x_light=0):
    T_ms = 5400 * (P*1e-9 / 140)**0.480 / (1 - np.log(1 - x_light))
    return T_ms

# Stevenson et al. 1983 
# Todo: find out span of P in model 
def T_m1983(P, x_light=0): 
    Tm0 = 2060              # Tm0 [K], Tm1 [K TPa^-1], Tm2 [K TPa1-2] are factors of the quadradic Pressure 
    Tm1 = 6.14*1e-12        # dependent liq. function (Stacey 1977b) (values for Earth with x = 0.1)
    Tm2 = -4.5*1e-24   
    ac = 0.1                # light element correction factor
    T_melt = Tm0 * (1 - ac * x_light) * (1 + Tm1 * P + Tm2 * P**2) 
    return T_melt

# Sixtrude 2014 based on Morard et al. 2011 with possible light element correction
def T_m2014(P, x_light=0): 
    T_melt = 6500 * (P*1e-9 / 340)**0.515 * 1 / (1 - np.log(1 - x_light))
    return T_melt

# DAC experiments by Fischer er al. 2013
def T_m2013(P): 
    T_melt = 1623 * ( (P*1e-9 - 1e-4)/23.6 + 1)**(1/1.89)
    return T_melt

# ab initio Li et al. 2020 120 - 256 GPa
def T_m2020(P): 
    T_melt = 4242 * ( (P - 120.6)/98 + 1)**(1/3.35)
    return T_melt

# ab initio by Kraus et al. 2022
# lasers and in situ x-ray diffraction for 260 - 1000 GPa (hcp)structure.
def T_m2022(P): 
    T_melt = 5530 * ((P - 260)/293 + 1)**0.552
    return T_melt
    
# ab initio by Gonzalez-Cataldo et al. 2023
# 300 - 5000 GPa
def T_m2023(P): 
    T_melt = 6469 * (1 + (P*1e-9 - 300)/434.82)**0.54369
    return T_melt