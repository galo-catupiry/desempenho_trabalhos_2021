"""
International Standard Atmosphere
Grupo E
"""

import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)


import math
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt

#constantes
rho0=1.225;# [kg/m^3]
p0=101325;# [Pa]
t0=288.15;# [K]
beta=-0.0065;# [K/m]
a0=334.1; # [m/s]
g=9.81; #[m/s^2]
R=287.053; # [J/kg.K]
gama = 1.4 # adimensional
referencia = 11000 # [m] separação entre as duas camadas

#Função que forneça [P,a,rho,T] dada uma altitude em metros
def atmo_isa(alti):
     
    if alti <= referencia: #trecho linear
        t=t0+beta*alti 
        p = (t/t0)**(-g/(beta*R)); p=p*p0 
        a_=math.sqrt(gama*R*t) 
        rho_ = (t/t0)**((-g/(beta*R))-1); rho_=rho_*rho0
                          
    else: #trecho constante
        t=216.65
        a_=math.sqrt(gama*R*t)
        p= math.exp((-g*(alti-referencia))/(R*t));p=p*22620.47
        rho_=math.exp((-g*(alti-referencia))/(R*t));rho_=rho_*0.3639
            
    return round(t,3),round(p,3),round(a_,3),round(rho_,3)

#Função que fornece a altitude dado uma entrada de 
#temperatura, pressao, velocidade do som ou densidade
def reverse_isa(t = None, p = None, a = None, rho = None):
    
    rho0 = 1.2250

    reverse_isa_data = pd.read_csv(current_dir + '\\reverse_isa.csv', sep = ';')
    
    reverse_isa_data.dropna(inplace = True)
    reverse_isa_data.columns = ['h', 'T', 'P', 'rho', 'a']
    
    reverse_isa_data['P'] = reverse_isa_data['P'].apply(lambda p: p * 100000)
    reverse_isa_data['rho'] = reverse_isa_data['rho'].apply(lambda rho: rho * rho0)
    
    
    if t != None:
        h_T = interp1d(reverse_isa_data['T'],reverse_isa_data['h'], 
                       kind = 'nearest-up', fill_value='extrapolate')
        
        return h_T(t)
    
    if p != None:
        h_P = interp1d(reverse_isa_data['P'], reverse_isa_data['h'],
                       fill_value='extrapolate')
        return h_P(p)
    
    
    if a != None:
        h_a = interp1d(reverse_isa_data['a'], reverse_isa_data['h'],
                       kind = 'nearest-up', fill_value='extrapolate')
        return h_a(a)
    
    if rho != None:
        h_rho = interp1d(reverse_isa_data['rho'], reverse_isa_data['h'],
                         fill_value='extrapolate')
        return h_rho(rho)
    
    

    

