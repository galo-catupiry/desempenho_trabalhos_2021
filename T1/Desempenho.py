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

from Interpolacao import CLp, Mp
from aircraft import jetstar



#%%

jet = jetstar(1)

beta = 9296
sigma = lambda h: np.exp((-h/beta))


#%% Alcance

#%%%  CL constante



#%%% V constante
















#%% Autonomia

#%%%  CL constante



#%%% V constante