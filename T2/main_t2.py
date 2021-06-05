"""
T2 - SAA0183 - Código Principal

Grupo E:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""

#%%Bibliotecas
import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import numpy as np
import Desempenho_T2 as DT2

#%% MAIN

# -/----------------- Propulsão ------------------/-
n = 0.85
T0 = 64000
c = 0.45  # (Verificar!)
zeta = 0.45  #  (Verificar!)

# -/-------- Análise de Alcance (Cruzeiro) -------/-
V1_cru =  811 / 3.6
h_cru = 13105

x1 = DT2.cruise_range('h_CL',V1_cru,h_cru,c, zeta)
x2 = DT2.cruise_range('V_CL',V1_cru,h_cru,c, zeta)
x3 = DT2.cruise_range('V_h',V1_cru,h_cru,c, zeta)

# -/--------- Parâmetros para diagramas ----------/-

# Diagrama T,D vs. V
Diagrama1 = True
if(Diagrama1):
    h_fig1 = [10000]
    V_fig1 = np.linspace(70,320,200)
    [D_total_fig1,Dmin_fig1] = DT2.total_drag(V_fig1,h_fig1)
    T_fig1 = DT2.jet_buoyancy(h_fig1,T0,n)
    
    figure_1 = DT2.TD_vs_V(h_fig1,V_fig1,D_total_fig1,T_fig1, Dmin_fig1)

# Diagrama h-V
Diagrama2 = False
if(Diagrama2):
    h_fig2 = np.arange(0,14400,10).tolist()
    V_fig2 = np.linspace(0,320,200)
    [D_total_fig2,Dmin_fig2] = DT2.total_drag(V_fig2,h_fig2)
    T_fig2 = DT2.jet_buoyancy(h_fig2,T0,n)
    
    V1_fig2 = DT2.cruise_velocity_solver(V_fig2,h_fig2,'V1',T0,n)
    V2_fig2 = DT2.cruise_velocity_solver(V_fig2,h_fig2,'V2',T0,n)
    
    figure_2 =  DT2.h_vs_V(h_fig2,V1_fig2,V2_fig2)
