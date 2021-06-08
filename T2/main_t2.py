"""
T2 - SAA0183 - C�digo Principal

Grupo E:
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
import Desempenho_T2 as DT2
from aircraft import JetStar

#%% MAIN
jet = JetStar(1)

# -/----------------- Propuls�o ------------------/-

n = 0.85
T0 = 64000
c = 0.0133  # (Verificar!)

# -/--------------- Pesos -----------------------/-

POV = 11566*9.81             # Peso vazio operacional, [N] 
MTOW = 20071.446*9.81        # Peso m�ximo de decolagem, [N]  
max_payload = 907.2*9.81     # M�xima carga paga, [N]
max_fuel = 8149.233872*9.81  # M�xima qtde. de combust�vel, [N]

# -/-------- An�lise de Alcance (Cruzeiro) -------/-

x1 = DT2.cruise_range_new('h_CL',MTOW, c, max_fuel/MTOW)
x2 = DT2.cruise_range_new('V_CL',MTOW, c, max_fuel/MTOW)
x3 = DT2.cruise_range_new('V_h', MTOW, c, max_fuel/MTOW)

# -/------------- Tetos da Aeronave ------------/-

V = np.linspace(50,320,600)
h = np.arange(0,14400,10).tolist()
tol = 0.015
resp = DT2.ceiling(h, T0, n, jet.W, V,tol)

#%% -/--------- Diagramas de Desempenho ----------/-

# Diagrama T,D vs. V
Diagrama1 = False
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

# Carga Paga vs. Alcance
Diagrama3 = False
if(Diagrama3):
    figure3 = DT2.payload_vs_range(c,POV,MTOW,max_payload,max_fuel)
    
# �ngulo de subida vs. Velocidade
Diagrama4 = False
if(Diagrama4):
    V_fig3 = np.linspace(0,320,200)
    h = [10000]
    figure4 = DT2.gamma_graph(h, T0, n, jet.W,V_fig3)

# Raz�o de subida vs. Velocidade
Diagrama5 = True
if(Diagrama5):
    h = [10000]
    figure5 = DT2.h_dot_vs_velocity(h, T0, n, jet.W)
     