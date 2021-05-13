"""
Cálculo dos parâmetros de desempenho

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
from Interpolacao import polar, param
from aircraft import JetStar


#%% Dados Gerais

jet = JetStar(1)
sealevel = Atmosphere(0)

beta = 9296
rho0 = sealevel.density[0]

#%% Alcance e Autonomia com CL constante
def alcance_autonomia_CL(altitude, V, 
                         graph_E_V = False, save_graph_E_V = False,
                         graph_h_V = False, save_graph_h_V = False):
    
    h1 = 0 # [m]
    h2 = altitude # [m]
    
    rho = Atmosphere(h2).density[0]
    
    CL = jet.W / (0.5 * jet.S * (V**2) * rho) #esse cl é mantido constante
    
    h_linspace = np.linspace(h2,h1,1000, retstep = True)
    range_h = h_linspace[0]
    range_dh = abs(h_linspace[1])
    
    deltaX_CL = 0
    t_CL = 0
    E_list = []
    V_list = []
    hdot_list = []
    
    for h in range_h:
        
        rho_i = Atmosphere(h).density[0]
        velo_som_i = Atmosphere(h).speed_of_sound[0]
        V_i = (jet.W / (0.5 * CL * rho_i * jet.S))**.5
                
        mach_i = V_i / velo_som_i
        
        CL_i = CL #mantendo o CL igual o CL inicial
        CD_i = polar(param, CL_i, mach_i)[0]
        E_i = CL_i / CD_i
        gamma_i = - 1 / E_i
        hdot_i = V_i * np.sin(gamma_i)
                
        #Alcance
        deltaX_CL_i = E_i * range_dh # E * [h(i+1) - h(i)]
        deltaX_CL += deltaX_CL_i
        
        E_list.append(E_i)
        V_list.append(V_i)
        hdot_list.append(hdot_i)
        
        
        #Autonomia
        exp_t =  np.e**(- (h - range_dh)/(2*beta)) - np.e**(-h/(2*beta))
        t_CL_i = 2*beta*E_i * (((rho0 * CL_i)/(2 * jet.WL))**.5) * exp_t
        t_CL += t_CL_i
        
    if graph_E_V == True:
        
        fig_CL= plt.figure(figsize=(10,4))
        plt.subplots_adjust(wspace = 0.3, hspace = 0.4)
        
        #Eficiencia pela altitude
        fig_E = fig_CL.add_subplot(121)
        plt.plot(range_h, E_list, label = "E")
        plt.xlabel("Altitude [m]", fontsize = 12)
        plt.ylabel("Eficiência ", fontsize = 12)
        plt.grid()
        ax = plt.gca()
        ax.invert_xaxis()
        
        #Velocidade pela altitude
        fig_V = fig_CL.add_subplot(122)
        plt.plot(range_h, V_list, label = "V")
        plt.grid()
        plt.xlabel("Altitude [m]", fontsize = 12)
        plt.ylabel("Velocidade [m/s] ", fontsize = 12)
        ax = plt.gca()
        ax.invert_xaxis()
        
        if save_graph_E_V == True:
            plt.savefig("alcance_autonomia_CL.pdf")
        
        plt.show()
        
    
    if graph_h_V == True:
        
        
        fig_h_v = plt.figure(figsize=(10,4))
        plt.plot(V_list, hdot_list)
        plt.grid()
        plt.xlabel("Velocidade [m/s]", fontsize = 12)
        plt.ylabel("Razão de descida $\dot{h}$", fontsize = 12)
        
        if save_graph_h_V == True:
            plt.savefig("h_V.pdf")
        plt.show()
    
    return deltaX_CL, t_CL


#%% Alcance e Autonomia com V constante
def alcance_autonomia_V(altitude, V, graph = False, save_graph = False):
    
    ## Considerando velocidade e altitude de cruzeiro
    h1 = 0 # [m]
    h2 = altitude # [m]

    #velo_som = ISA(h2)[2] # [m/s]
    velo_som = Atmosphere(h2).speed_of_sound[0]
    mach_cru = V / velo_som
    
    h_linspace = np.linspace(h2,h1,1000, retstep = True) 
    range_h = h_linspace[0]
    range_dh = abs(h_linspace[1])
    
    deltaX_V = 0
    t_V = 0
    E_list = []
    CL_list = []
    
    for h in range_h:
        
        rho_i = Atmosphere(h).density[0]
        CL_i = (2 * jet.WL)/(rho_i * V**2)
        CD_i = polar(param, CL_i, mach_cru)[0]
        E_i = CL_i / CD_i
        CD0_i = polar(param, CL_i, mach_cru)[1]
        K_i = polar(param, CL_i, mach_cru)[2]
        
        A_i = (rho0 * CD0_i * V**2)/(2 * jet.WL)
        B_i = (2 * K_i * jet.WL)/(rho0 * V**2)
        
        tan1_i = np.arctan(B_i**-1 * A_i * np.e**(-(h - range_dh)/beta))
        tan2_i = np.arctan(B_i**-1 * A_i * np.e**(-h/beta))
        
        #Alcance
        deltaX_V_i = (beta/B_i)* (tan1_i - tan2_i)
        deltaX_V += deltaX_V_i
        
        #Autonomia
        t_V_i = beta/(B_i * V) * (tan1_i - tan2_i)
        t_V += t_V_i
        
        E_list.append(E_i)
        CL_list.append(CL_i)
        
        
    if graph == True:
        
        fig_V= plt.figure(figsize=(10,4))
        plt.subplots_adjust(wspace = 0.3, hspace = 0.4)
        
        
        #Eficiencia pela altitude
        fig_E = fig_V.add_subplot(121)
        plt.plot(range_h, E_list, label = "E")
        plt.xlabel("Altitude [m]", fontsize = 12)
        plt.ylabel("Eficiência ", fontsize = 12)
        plt.grid()
        ax = plt.gca()
        ax.invert_xaxis()
        
        #CL pela altitude
        fig_CL = fig_V.add_subplot(122)
        plt.plot(range_h, CL_list, label = "CL")
        plt.grid()
        plt.xlabel("Altitude [m]", fontsize = 12)
        plt.ylabel("$C_L$", fontsize = 12)
        ax = plt.gca()
        ax.invert_xaxis()
        
        if save_graph == True:
            plt.savefig("alcance_autonomia_V.pdf")
        
        plt.show()
        
    return deltaX_V, t_V
