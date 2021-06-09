"""
T2 - SAA0183 - Código Principal

Grupo E:
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
import Cruzeiro as cr
import Voo_ascendente as v_asc
from aircraft import JetStar

# =============================================  
jet = JetStar(1)

# Propulsão
n = 0.85
T0 = 64000
c = 0.0133

# Pesos
POV = 11566*9.81             # Peso vazio operacional, [N] 
MTOW = 20071.446*9.81        # Peso máximo de decolagem, [N]  
max_payload = 907.2*9.81     # Máxima carga paga, [N]
max_fuel = 8149.233872*9.81  # Máxima qtde. de combustível, [N]

# Análise de Alcance (Cruzeiro)
x1 = cr.cruise_range_new('h_CL',MTOW, c, max_fuel/MTOW)
x2 = cr.cruise_range_new('V_CL',MTOW, c, max_fuel/MTOW)
x3 = cr.cruise_range_new('V_h', MTOW, c, max_fuel/MTOW)

# Tetos
V = np.linspace(50,320,600)
h = np.arange(0,14400,10).tolist()
tol = 0.015
resp = v_asc.ceiling(h, T0, n, jet.W, V,tol)

# =============================================  
# Gráficos

# Diagrama T,D vs. V
fig1 = True
if(fig1):
    h_fig1 = [10000]
    V_fig1 = np.linspace(70,320,200)
    [D_total_fig1,Dmin_fig1] = cr.total_drag(V_fig1,h_fig1)
    T_fig1 = cr.jet_buoyancy(h_fig1,T0,n)
    
    figure_1 = cr.TD_vs_V(h_fig1,V_fig1,D_total_fig1,T_fig1, Dmin_fig1)

# Diagrama h-V
fig2 = False
if(fig2):
    h_fig2 = np.arange(0,14400,10).tolist()
    V_fig2 = np.linspace(0,320,200)
    [D_total_fig2,Dmin_fig2] = cr.total_drag(V_fig2,h_fig2)
    T_fig2 = cr.jet_buoyancy(h_fig2,T0,n)
    
    V1_fig2 = cr.cruise_velocity_solver(V_fig2,h_fig2,'V1',T0,n)
    V2_fig2 = cr.cruise_velocity_solver(V_fig2,h_fig2,'V2',T0,n)
    
    figure_2 =  cr.h_vs_V(h_fig2,V1_fig2,V2_fig2)

# Carga Paga vs. Alcance
fig3 = False
if(fig3):
    figure3 = cr.payload_vs_range(c,POV,MTOW,max_payload,max_fuel)
    
# Ângulo de subida vs. Velocidade
fig4 = False
if(fig4):
    V_fig3 = np.linspace(0,320,200)
    h = [10000]
    figure4 = v_asc.gamma_graph(h, T0, n, jet.W,V_fig3)

# Razão de subida vs. Velocidade
fig5 = True
if(fig5):
    h = [12000]
    figure5 = v_asc.h_dot_vs_velocity(h, T0, n, jet.W)
     