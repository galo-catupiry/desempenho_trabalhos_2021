"""
Análise de Decolagem - T2

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""
# =============================================  
import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import numpy as np
import matplotlib.pyplot as plt

from aircraft import JetStar
from ambiance import Atmosphere
from Interpolacao import DragPolar
from Cruzeiro import jet_buoyancy

# =============================================
jet = JetStar('power approach')
sealevel = Atmosphere(0)
rho0 = sealevel.density[0]
g = 9.81
ft_to_m = 0.3048
# =============================================

CLmax = 2*1*jet.W/(64**2*Atmosphere(0).density[0]*jet.S)
# Estol

R_gas = 287.058 #J/(kg*K) - constante dos gases

# Propulsao
n = 0.85
T0 = 64000
c = 0.5/3600

MTOW = 18500*g #temporário

def air_density(Temp, h):
    '''densidade fora da ISA; Temp em Kelvin'''
    P=Atmosphere(h).pressure
    Temp_isa_mod = Atmosphere(h).temperature + Temp-(273.15 + 15)
    rho = P/(R_gas*Temp_isa_mod) 
    return (rho/rho0)[0]

#rho = air_density(273.15-20, 1000)
#print(rho)

drag = DragPolar()


def ground_acceleration(rho, V, W, CLmax, gamma_pista=0, mi=0.04):
    
    L = 0.5*CLmax*jet.S*rho*V**2   
    
    D = 0.5*jet.CD*jet.S*rho*V**2
    T = T0*((rho/rho0)**n)
    Fric = mi*(W*np.cos(gamma_pista) - L)

    a = (T - D - Fric - jet.W*np.sin(gamma_pista))*(g/W)
    gamma = (T-D)/W

    return a, gamma
    
def running(rho,V_S, W, CLmax, gamma_pista=0, mi=0.04):
    global V_Lo
    
    V_Lo = 1.2*V_S
    a_bar = ground_acceleration(rho, 0.7*V_Lo, W, CLmax, gamma_pista=0, mi=0.04)[0]
    x_G = (V_Lo**2)/(2*a_bar)
    return x_G, V_Lo

def rotation(V_S):
    V_R = 1.1*V_S
    t_R = 3  
    x_R =t_R*V_R  
    return x_R
   
def transition_and_climbing(V_Lo, W, rho):
    #Transição e Subida
    R_Tr = (V_Lo**2)/(0.15*g)
    h_Sc = 35*ft_to_m
    gamma = ground_acceleration(rho, 0.7*V_Lo, W ,CLmax, 0, 0.04)[1]
    h_Tr = (1-np.cos(gamma))*R_Tr
    if h_Sc > h_Tr:
        x_Tr = R_Tr*np.sin(gamma)
        x_Cl = (h_Sc - h_Tr)/np.tan(gamma)    
    else:
        gamma= np.arccos(-h_Sc/R_Tr +1)
        x_Tr = R_Tr*np.sin(gamma)
        x_Cl = 0

    return x_Tr, x_Cl
