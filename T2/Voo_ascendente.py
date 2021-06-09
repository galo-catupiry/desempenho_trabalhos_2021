"""
Análise de Voo Ascendente - T2

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""
# =============================================
from aircraft import JetStar
from Interpolacao import DragPolar
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from ambiance import Atmosphere

import Cruzeiro as cr

# =============================================
jet = JetStar(1)

sealevel = Atmosphere(0)
beta = 9296
rho0 = sealevel.density[0]
# =============================================

def gamma(h,T0,n,W,V):
    
    T = cr.jet_buoyancy(h, T0, n)[0]
    D = cr.total_drag(V, h)[0][0]
    return (T - D)/W

def h_dot(h,T0,n,W,V):
    
    T = cr.jet_buoyancy(h, T0, n)[0]
    D = cr.total_drag(V, h)[0][0]
    h_dot = (T*V - D*V)/W
    
    return h_dot

def ceiling(h,T0,n,W,V,tol):
    
    for i in h:
        aux = h_dot([i], T0, n, W, V)
        h_dot_max = max(aux)
        
        if (abs(h_dot_max - 1.524) <= tol):
            print("Teto operacional: h = {:.2f}".format(i))
            #operating_ceiling = i
            
        elif (abs(h_dot_max - 0.508) <= tol):
            print("Teto de serviço: h = {:.2f}".format(i))
            #service_ceiling = i
            
        elif(abs(h_dot_max) <= tol):
            print("Teto absoluto: h = {:.2f}".format(i))
            #absolute_ceiling = i
    
    return 

# TODO: parâmetros ótimos
def optimal_parameters():
    
    # Vel. para máxima razão de subida
     
    return

# =============================================
# Gráficos

def gamma_graph(h,T0,n,W,V):
    
    plt.style.use('ggplot')
    
    plt.figure(4)
    plt.ylabel("$\\gamma$", fontsize = 12)
    plt.xlabel("Velocity [m/s]", fontsize = 12)
    plt.grid(True)
    
    V1 = cr.cruise_velocity_solver(V, h, 'V1', T0, n)
    V2 = cr.cruise_velocity_solver(V, h, 'V2', T0, n)
    
    V_plot = np.linspace(V1,V2,300,endpoint=True)
    plt.plot(V_plot,gamma(h,T0,n,W,V_plot),color = 'purple')
    plt.show()

    return

def h_dot_vs_velocity(h, T0, n, W):
    
    plt.style.use('ggplot')
    
    plt.figure(5)
    plt.ylabel("$\dot{h} \:\: [m/s]$",fontsize = 12)
    plt.xlabel("Velocity  [m/s]",fontsize = 12)
    plt.grid(True)
    
    V = np.linspace(0,320,200)
    V1 = cr.cruise_velocity_solver(V, h, 'V1', T0, n)
    V2 = cr.cruise_velocity_solver(V, h, 'V2', T0, n)
    
    gamma_list = gamma(h,T0,n,W,V)
    gamma_max = max(gamma_list)   
    
    plt.plot(V[:170], V[:170]*gamma_max,"--k", label = '$\gamma_{máx}$')
      
    V_plot = np.linspace(V1,V2,300,endpoint=True)
    plt.plot(V_plot, h_dot(h, T0, n, W, V_plot),color = 'purple')
    legend = plt.legend(loc = 'best', framealpha = 1)
    plt.setp(legend.get_texts(), color='k')
    plt.show()
    
    return