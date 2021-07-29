"""
AnÃ¡lise de Pouso - T2

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

from aircraft import JetStar
from ambiance import Atmosphere

# =============================================
jet = JetStar('power approach')
sealevel = Atmosphere(0)
rho0 = sealevel.density[0]
g = 9.81
ft_to_m = 0.3048
R_gas = 287.058 #J/(kg*K) - constante dos gases

# Propulsao
n = 0.85
T0 = 64000
c = 0.5/3600
# =============================================


def air_density(Temp, h):
    '''densidade fora da ISA; Temp em Celcius'''
    P=Atmosphere(h).pressure[0]
    Temp_isa_mod = Atmosphere(h).temperature[0] + Temp-(15)
    rho = P/(R_gas*Temp_isa_mod) 
    return rho

def approach_and_flare(rho, W_ap, V_S):
    V_ap = 1.3*V_S
    
    D = 0.5*jet.CD*jet.S*rho*V_ap**2
    T = T0*((rho/rho0)**n)
    
    R_f = V_ap**2/(0.08*g)   
    h_Sc = 50*ft_to_m
    gamma_ap = (T-D)/W_ap #gamma máximo pela potência 
    h_f = (1-np.cos(gamma_ap))*R_f

    if h_f<h_Sc:                #com altos gammas a aeronave atinge a altura sem a fase de appoach, 
        x_f = R_f*np.sin(gamma_ap)            #limita-se o gamma por h_f=h_Sc
        x_ap = (h_Sc - h_f)/np.tan(gamma_ap)    
    else:
        gamma_ap= np.arccos(-h_Sc/R_f +1)
        x_f = R_f*np.sin(gamma_ap)
        x_ap = 0
    return x_f, x_ap

def rotation(V_S):
    V_R = 1.1*V_S
    t_R = 3  
    x_R =t_R*V_R  
    return x_R

def landing_roll(V_S):
    V_Td = 1.1*V_S 
    d_bar = 0.55*g
    x_gLa = V_Td**2/(2*d_bar) 
    return x_gLa
