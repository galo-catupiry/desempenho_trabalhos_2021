"""
Analise de Voo Ascendente - T2

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

from aircraft import JetStar
import numpy as np
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
            print("Teto operacional: h = {:.2f} [ft]".format(i*3.28084))
            #operating_ceiling = i
            
        elif (abs(h_dot_max - 0.508) <= tol):
            print("Teto de servico: h = {:.2f} [ft]".format(i*3.28084))
            #service_ceiling = i
            
        elif(abs(h_dot_max) <= tol):
            print("Teto absoluto: h = {:.2f} [ft]".format(i*3.28084))
            #absolute_ceiling = i
    
    return 

# TODO: Parametros otimos
def optimal_parameters():
    
    # Vel. para maxima razao de subida
     
    return

# =============================================
# Graficos

def gamma_graph(h,T0,n,W,V):
    
    plt.style.use('default')
    
    plt.figure()
    plt.ylabel("$\\gamma \: [\\degree]$", fontsize = 12)
    plt.xlabel("Velocity [m/s]", fontsize = 12)
    plt.grid(False)
    
    V1 = cr.cruise_velocity_solver(V, h, 'V1', T0, n)
    V2 = cr.cruise_velocity_solver(V, h, 'V2', T0, n)
    
    V_plot = np.linspace(V1,V2,300,endpoint=True)
    gamma_plot = [i*180/np.pi for i in gamma(h,T0,n,W,V_plot)]
    plt.plot(V_plot,gamma_plot,color = 'purple', label = 'h = {:.1f} km'.format(h[0]/1000))
    plt.legend(loc = 'best', framealpha = 1)
    plt.savefig('gamma_ascendente.svg')
    plt.show()

    return

def h_dot_vs_velocity(h, T0, n, W):
    
    plt.style.use('default')
    
    plt.figure()
    plt.ylabel("$\dot{h} \:\: [m/s]$",fontsize = 12)
    plt.xlabel("Velocity  [m/s]",fontsize = 12)
    plt.grid(False)
    
    for i in h:
        V = np.linspace(0,320,200)
        V1 = cr.cruise_velocity_solver(V, [i], 'V1', T0, n)
        V2 = cr.cruise_velocity_solver(V, [i], 'V2', T0, n)
        
        gamma_list = gamma([i],T0,n,W,V)
        gamma_max = max(gamma_list)   
        
        #plt.plot(V[:170], V[:170]*gamma_max,"--k", label = '$\gamma_{máx}$')
          
        V_plot = np.linspace(V1,V2,300,endpoint=True)
        plt.plot(V_plot, h_dot([i], T0, n, W, V_plot), label = 'h = {:.1f} [km]'.format(i/1000))
    legend = plt.legend(loc = 'best', framealpha = 1)
    plt.setp(legend.get_texts(), color='k')
    plt.ylim(bottom = 0)
    plt.savefig('razaoSubida_velocidade.svg')
    plt.show()

    
    return