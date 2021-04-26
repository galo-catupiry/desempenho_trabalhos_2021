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



#%% Dados Gerais

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


def hdot_v_graph():
    fig = plt.figure(figsize=(6,6))
    
    h_dot_list = []
    V_list = []
    gamma_line = []
    
    for i in np.linspace(1, velo_cruzeiro,500):
        densidade = ISA(altitude)[3] #[kg/m^3]
        velo_som = ISA(altitude)[2] # [m/s]
        mach = i/velo_som

        CL = jet.W / (0.5 * jet.S * (i**2) * densidade)
        CD = polar(param, CL, mach)

        E = CL/CD
        gamma = -1/E
        h_dot = i * gamma
        
        h_dot_list.append(h_dot)
        V_list.append(i)
        
    plt.plot(V_list, h_dot_list, color = 'k')
    plt.grid()
    plt.xlabel("V [m/s]", fontsize = 14)
    plt.ylabel("$\dot{h}$", fontsize = 14)
    plt.show()
        
    
beta = 9296
sigma = lambda h: np.exp((-h/beta))


h2 = altitude # [m]
h1 = 0 # [m]

#%% Alcance e Autonomia com CL constante

#%%% Alcance e Autonomia

#Alcance
def range_CL(h1, h2, CL):
    E = CL #alguma funcao com CL
    deltaX_CL = E * (h2 - h1)
    
    return deltaX_CL

#Autonomia
def endurance_CL(h1,h2,CL):
    #formula 7.24
    E = CL #alguma funcao com CL
    rho_ssl = ISA(2) #conferir isso
    expo = np.e**(-h1/(2*beta)) - np.e**(-h2/(2*beta))
    t_CL = (2 * beta* E * (rho_ssl * CL/(2*jet.WL))**.5) * expo
    
    return t_CL

#%%% Melhot Alcance e melhor Autonomia




    
    
















#%% Alcance e Autonomia com V constante

#%%%  Alcance e Autonomia

#Alcance
def range_V(h1, h2, V):
    #formula 7.20
    rho_ssl = 0
    CD0 = 0
    K = 0
    A = (rho_ssl * CD0 * V**2)/(2 * jet.WL)
    B = (2 * K * jet.WL)/(rho_ssl * V**2)
    tan1 = np.arctan(B**-1 * A * np.e**(-h1/beta))
    tan2 = np.arctan(B**-1 * A * np.e**(-h2/beta))
    
    deltaX_V = (beta / B) * (tan1 - tan2)
    
    return deltaX_V

#Autonomia
def endurance_V(h1, h2, V):
    rho_ssl = 1
    V = 1
    CD0 = 0
    K = 0
    
    a = (CD0 * rho_ssl**2 * np.e**(-h1/beta))/(4*K * jet.WL**2)
    b = (CD0 * rho_ssl**2 * np.e**(-h2/beta))/(4*K * jet.WL**2)
    
    arc_tan = np.arctan(((a-b)*(V**4)) / (1 + a*b*(V**8))
                        
    tv = ((beta*rho_ssl* V)/(2*K * jet.WL)) * arc_tan
    
    
    
#%%% Melhor Alcance e Autonomia 
    