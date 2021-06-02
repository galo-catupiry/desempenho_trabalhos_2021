"""
Cálculo dos parâmetros de desempenho - T2

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""

#%% Packages

import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from Interpolacao import DragPolar
from aircraft import JetStar
from scipy.optimize import fsolve
from scipy.optimize import brentq

plt.close('all')

#%% Dados Gerais

jet = JetStar(1)

sealevel = Atmosphere(0)
beta = 9296
rho0 = sealevel.density[0]

#%% Functions for Cruise Flight

# DragPolar object
drag_object = DragPolar()
drag_object.CLp = 0

# Total Drag
def total_drag(V,h):
    
    D_resp, Dmin_resp = [], []
    
    for i in h:
        sigma = Atmosphere(i).density[0]/rho0
        drag_object.Mp = V/Atmosphere(i).speed_of_sound[0]
        
        drag_polar = drag_object.polar()
        CD0 = drag_object.CD0
        K = drag_object.K
        
        Dp = (1/2)*sigma*rho0*(V**2)*jet.S*CD0
        Di = (2*K*(jet.W**2))/(sigma*rho0*jet.S*V**2)
        D_total = Dp + Di
        
        Emax = np.sqrt(CD0/K)/(2*CD0)
        Dmin = jet.W/Emax
        
        D_resp.append(D_total)
        Dmin_resp.append(Dmin)
    
    return D_resp, Dmin_resp

# Cruise Velocity equations
def cruise_velocity_eq(x,h):
    
    CD0 = x[0]
    K = x[1]
    V = x[2]
    
    sigma = Atmosphere(h).density[0]/rho0
    drag_polar = drag_object.polar()
    drag_object.Mp = V/Atmosphere(h).speed_of_sound[0]
    Emax = np.sqrt(CD0/K)/(2*CD0)
    
    # Equations to be solved
    CD0_exp = drag_object.CD0 - CD0
    K_exp = drag_object.K - K
    V_exp = (1/2)*sigma*rho0*jet.S*CD0*V**4 - T0*sigma**n*V**2 + (2*K*jet.W**2)/(sigma*rho0*jet.S) 

    equations = [V_exp,CD0_exp,K_exp]
    
    return equations   
    
    
# Cruise velocity solver
def cruise_velocity_solver(V,h,V_type):
    
    Vresp = []
    D_total = total_drag(V,h)[0]
    T = jet_buoyancy(h, T0)
    
    for j in np.arange(0,len(h)):
        
        d = T[j] - D_total[j]
            
        for i in np.arange(0,len(D_total[0])-1):
            if (d[i] < 0 and d[i+1] > 0):
                V1_0 = V[i]
            elif (d[i] > 0 and d[i+1] < 0):
                V2_0 = V[i]
        
        
        if (V_type == 'V1'):
                initial = (0.01,0.01,V1_0)
        elif (V_type == 'V2'):
                initial = (0.01,0.01,V2_0)
        
        [CD0_resp,K_resp,V_resp] = fsolve(cruise_velocity_eq, initial, args = (h[j]))
        
        Vresp.append(V_resp)
    
    return Vresp
    
# Buoyancy (jet)
def jet_buoyancy(h,T0):
    
    T = []
    
    for i in h:
        sigma = Atmosphere(i).density[0]/rho0
        Ti = T0*(sigma**n) 
        T.append(Ti)
    
    return T


#%% Plots for Cruise Flight

def TD_vs_V(V,D_total,T, Dmin):
    
    plt.figure(1)
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("T and D")
    plt.grid(True)
    color=iter(plt.cm.brg(np.linspace(0,1,len(h))))
    
    if (len(h) > 1):    
        
        for i in np.arange(0,len(T)):
            c=next(color)
            plt.plot(V,D_total[i], color = c)
            plt.plot(V,[T[i]]*len(V), label = 'T,D (h = {} [m])'.format(h[i]),color = c)
            #plt.plot(V,Dmin[i],'--k')
    else:
        
        D_total = D_total[0]
        Dmin = Dmin[0]
        
        plt.plot(V,D_total,'k', label = 'Drag')
        plt.plot(V,[T]*len(V), label = 'Thrust')
        plt.plot(V,Dmin,'--k',label = 'Minimum Drag')
        

    plt.legend(loc = 'best', framealpha = 1)
    plt.ylim(top = 40000)
    return

def h_vs_V(V1,V2):
    # EM CONSTRUÇÃO!!
    
    plt.figure(2)
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("h [m]")
    plt.grid(True)
    plt.plot(V1,h,'k')
    plt.plot(V2,h,'k')
    #plt.legend(loc = 'best')
    return

#%% MAIN

n = 0.85
h = np.arange(0,14430,100).tolist()
T0 = 64000

V = np.linspace(0,320,200)
[D_total,Dmin] = total_drag(V,h)
T = jet_buoyancy(h,T0)

V1 = cruise_velocity_solver(V,h,'V1')
V2 = cruise_velocity_solver(V,h,'V2')

# Figures:
    
#figure_1 = TD_vs_V(V,D_total,T, Dmin)
figure_2 =  h_vs_V(V1,V2)

