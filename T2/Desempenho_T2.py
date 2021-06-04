"""
Cálculo dos parâmetros de desempenho - T2

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""

#%% Bibliotecas
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

#%% Funções para Voo em Cruzeiro

# Polar de arrasto
drag_object = DragPolar()
drag_object.CLp = 0

def total_drag(V,h):
    '''
    
    Parâmetros
    ----------
    V : Lista de velocidades definida em "MAIN".
    h : Lista de altitudes definida em "MAIN".

    Returns
    -------
    D_resp : Arrasto total em cruzeiro (lista)
    Dmin_resp : Arrasto mínimo em cruzeiro (lista)

    '''
    
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

def cruise_velocity_eq(x,h):
    '''

    Parâmetros
    ----------
    x : Variável do sistema
    h : Altura analisada (lista)

    Retorna
    -------
    equations : Lista de equações a serem resolvidas

    '''
    
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
    
def cruise_velocity_solver(V,h,V_type):
    '''
    
    Parâmetros
    ----------
    V : Lista de velocidades definida em "MAIN".
    h : Lista de altitudes definida em "MAIN".
    V_type : Define a análise para V1 ou V2.

    Returns
    -------
    Vresp : V1 ou V2 (lista, pois variam com h)

    '''
    
    Vresp = []
    D_total = total_drag(V,h)[0]
    T = jet_buoyancy(h, T0)
    
    CD0_resp = 0.01
    K_resp = 0.01
    
    for j in np.arange(0,len(h)):
        
        d = T[j] - D_total[j]
            
        for i in np.arange(0,len(D_total[0])-1):
            if (d[i] < 0 and d[i+1] > 0):
                V1_0 = V[i]
            elif (d[i] > 0 and d[i+1] < 0):
                V2_0 = V[i]
        
        
        if (V_type == 'V1'):
                initial = (CD0_resp,K_resp,V1_0)
        elif (V_type == 'V2'):
                initial = (CD0_resp,K_resp,V2_0)
        
        [CD0_resp,K_resp,V_resp] = fsolve(cruise_velocity_eq, initial, args = (h[j]))
        
        Vresp.append(V_resp)
    
    return Vresp
    
def jet_buoyancy(h,T0):
    '''

    Parâmetros
    ----------
    h : Lista de altitudes definida em "MAIN".
    T0 : Empuxo dos quatro motores da aeronave ao nível do mar

    Returns
    -------
    T : Lista de empuxo para cada altitude (h).

    '''
    
    T = []
    
    for i in h:
        sigma = Atmosphere(i).density[0]/rho0
        Ti = T0*(sigma**n) 
        T.append(Ti)
    
    return T

def cruise_range(cond,V1,h,c, zeta):
    drag_cru = DragPolar()
    rho = Atmosphere(h).density[0]
    
    CL_cru = (2*jet.W)/(rho*V1**2*jet.S)
    drag_cru.CLp = CL_cru
    drag_cru.Mp = V1/Atmosphere(h).speed_of_sound[0]
    CD1 = drag_cru.polar()
    
    E1 = CL_cru/CD1
    Em = 4*drag_cru.K/drag_cru.CD0
    
    if(cond == 'h_CL'):
        
        x = (2*V1*E1)/c*(1-np.sqrt(1 - zeta))
    
    elif(cond == 'V_CL'):
        
        x = E1*V1/c*np.log(1/(1 - zeta))
        
    elif(cond == 'V_h'):
       
        x = (2*V1*Em)/c*np.arctan(E1*zeta/(2*Em*(1 - drag_cru.K*E1*CL_cru*zeta)))
        
    return x

#%% Plots for Cruise Flight

def TD_vs_V(h,V,D_total,T, Dmin):
    
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

def h_vs_V(h,V1,V2):
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

# -/----------------- Propulsão ------------------/-
n = 0.85
T0 = 64000
c = 0.45  # (Verificar!)
zeta = 0.45  #  (Verificar!)

# -/-------- Análise de Alcance (Cruzeiro) -------/-
V1_cru =  811 / 3.6
h_cru = 13105

x1 = cruise_range('h_CL',V1_cru,h_cru,c, zeta)
x2 = cruise_range('V_CL',V1_cru,h_cru,c, zeta)
x3 = cruise_range('V_h',V1_cru,h_cru,c, zeta)

# -/--------- Parâmetros para diagramas ----------/-

# Diagrama T,D vs. V
Diagrama1 = False
if(Diagrama1):
    h_fig1 = [10000, 8000]
    V_fig1 = np.linspace(70,320,200)
    [D_total_fig1,Dmin_fig1] = total_drag(V_fig1,h_fig1)
    T_fig1 = jet_buoyancy(h_fig1,T0)
    
    figure_1 = TD_vs_V(h_fig1,V_fig1,D_total_fig1,T_fig1, Dmin_fig1)

# Diagrama h-V
Diagrama2 = True
if(Diagrama2):
    h_fig2 = np.arange(0,14800,10).tolist()
    V_fig2 = np.linspace(0,320,200)
    [D_total_fig2,Dmin_fig2] = total_drag(V_fig2,h_fig2)
    T_fig2 = jet_buoyancy(h_fig2,T0)
    
    V1_fig2 = cruise_velocity_solver(V_fig2,h_fig2,'V1')
    V2_fig2 = cruise_velocity_solver(V_fig2,h_fig2,'V2')
    
    figure_2 =  h_vs_V(h_fig2,V1_fig2,V2_fig2)


