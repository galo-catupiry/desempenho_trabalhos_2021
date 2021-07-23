"""
T2 - SAA0183 - Codigo Principal

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
import matplotlib.pyplot as plt

import Cruzeiro as cr
import Voo_ascendente as v_asc
import Manobras as manobras
from aircraft import JetStar
from ambiance import Atmosphere

plt.close('all')
# =============================================  
jet = JetStar(1)

# Estol
CLmax = 2*1*jet.W/(64**2*Atmosphere(0).density[0]*jet.S)

# Propulsao
n = 0.85
T0 = 64000
c = 0.5/3600

# Pesos
POV = 11566*9.81             # Peso vazio operacional, [N] 
MTOW = 20071.446*9.81        # Peso maximo de decolagem, [N]  
max_payload = 907.2*9.81     # Maxima carga paga, [N]
max_fuel = 8149.233872*9.81  # Maxima qtde. de combustivel, [N]

# Analise de Alcance (Cruzeiro)

# Tetos
V = np.linspace(50,320,600)
h = np.arange(0,14400,10).tolist()
tol = 0.015
resp = v_asc.ceiling(h, T0, n, jet.W, V,tol)


# =============================================  
# Graficos

# Diagrama T,D vs. V
fig1 = False
if(fig1):
    
    h_fig1 = [10000]
    V_fig1 = np.linspace(70,320,200)
    T_fig1 = cr.jet_buoyancy(h_fig1,T0,n)
    
    [D_total_fig1,Dmin_fig1] = cr.total_drag(V_fig1,h_fig1)

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
    V_s_fig2 = cr.estol(jet.W, jet.S, h_fig2, CLmax)
    print(V_s_fig2)
    figure_2 =  cr.h_vs_V(h_fig2,V1_fig2,V2_fig2, V_s_fig2)

# Carga Paga vs. Alcance
fig3 = False
if(fig3):
    figure_3 = cr.payload_vs_range(c,POV,MTOW,max_payload,max_fuel,200, 10000)
    
# Angulo de subida vs. Velocidade
fig4 = False
if(fig4):
    V_fig3 = np.linspace(0,320,200)
    h = [10000]
    figure_4 = v_asc.gamma_graph(h, T0, n, jet.W,V_fig3)

# Razao de subida vs. Velocidade
fig5 = False
if(fig5):
    h = [10000, 12000]
    figure_5 = v_asc.h_dot_vs_velocity(h, T0, n, jet.W)

# Diagrama T,D vs. V para curva coordenada
fig6 = False
if(fig6):
    
    h_fig6 = [10000]
    fc_fig6 = np.linspace(1,10, 1000, endpoint = True)
    V_fig6 = np.linspace(20,320,200)
    fc_max = manobras.fc_maximo(fc_fig6, V_fig6, h_fig6, T0, n)
    
    fc_fig6 = np.linspace(1, fc_max, 5, endpoint = True)
    T_manobra = manobras.T(h_fig6, T0, n, fc_fig6)
    D_manobra = manobras.drag(V_fig6, h_fig6, fc_fig6)
    
    figure_6 = manobras.TD_vs_V_manobras(h_fig6, fc_fig6, V_fig6, D_manobra, T_manobra)

# Diagrama de desempenho em curva
fig7 = True
if(fig7):

    # Dados
    h_fig7 = [0]
    V_fig7 = np.linspace(20,320,200)
    fc_fig7 = np.linspace(1,10,1000, endpoint = True)
    fc_max = manobras.fc_maximo(fc_fig7, V_fig7, h_fig7, T0, n)
    
    # Velocidades de Cruzeiro
    V1_fig7 = cr.cruise_velocity_solver(V_fig7,h_fig7,'V1',T0,n)
    V2_fig7 = cr.cruise_velocity_solver(V_fig7,h_fig7,'V2',T0,n)

    # Velocidades de Manobra
    fc_fig7 = np.linspace(1, fc_max, 100, endpoint = True)
    [V1_manobras,V2_manobras] = manobras.velocidades_manobra_solver(V1_fig7[0], V2_fig7[0], h_fig7, T0, n , fc_fig7)
    omega1_manobras = manobras.omega(fc_fig7, V1_manobras)
    omega2_manobras = manobras.omega(fc_fig7, V2_manobras)
    omega_V = (9.81/V_fig7)*np.sqrt(fc_max**2 - 1)

    # Velocidade de Estol
    fc_s = np.linspace(1,5, 100, endpoint = True)
    [Vs, omega_s] = manobras. estol(fc_s, jet.W, jet.S, h_fig7, CLmax)

    # Velocidade de Mergulho
    Vd = 1.4*V2_fig7[0]

    # Figura
    figure_7 = manobras.omega_vs_V(V1_manobras, V2_manobras, omega1_manobras, omega2_manobras, 
                                    fc_fig7, h_fig7, V_fig7, omega_V, fc_max, Vs, omega_s, Vd)
