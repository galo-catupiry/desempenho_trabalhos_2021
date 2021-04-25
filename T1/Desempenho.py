"""
T1 - SAA0183
Grupo E
"""


#%% Packages

import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import numpy as np
import matplotlib.pyplot as plt

from Interpolacao import polar, param
from aircraft import jetstar
from Utilities.isa_atmosphere import atmo_isa as ISA



#%%

jet = jetstar(1)

velo_cruzeiro = 811 / 3.6 # [m/s]
altitude = 13105 # [m]
densidade = ISA(altitude)[3] #[kg/m^3]
velo_som = ISA(altitude)[2] # [m/s]
mach = velo_cruzeiro/velo_som

CL = jet.W / (0.5 * jet.S * (velo_cruzeiro**2) * densidade)
CD = polar(param, CL, mach)

E = CL/CD
gamma = -1/E
h_dot = velo_cruzeiro * gamma




    
    
    




beta = 9296
sigma = lambda h: np.exp((-h/beta))


h2 = altitude # [m]
h1 = 0 # [m]



#%% Alcance

#%%%  CL constante





#%%% V constante
















#%% Autonomia

#%%%  CL constante



#%%% V constante