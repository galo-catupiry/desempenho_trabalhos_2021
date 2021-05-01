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

from Interpolacao import polar, param
from aircraft import JetStar
from Utilities.isa_atmosphere import atmo_isa as ISA

#%% Dados Gerais

jet = JetStar(1)

beta = 9296
rho0 = ISA(0)[3] #densidade ao nível do mar
sigma = lambda h: np.exp((-h/beta))

#%% Alcance e Autonomia com CL constante

#%%% Alcance e Autonomia

#Alcance
def alcance_autonomia_CL(altitude, V, graph = False, save_graph = False):
    
    ## Considerando velocidade e altitude de cruzeiro
    h1 = 0 # [m]
    h2 = altitude # [m]
    
    rho = ISA(h2)[3] #[ kg/m^3]
    
    CL = jet.W / (0.5 * jet.S * (V**2) * rho) #esse cl é mantido constante
    
    range_h = np.linspace(h2,h1,10000)
    range_dh = range_h[0] - range_h[1] # da um dh próximo de 13.11 m
    
    deltaX_CL = 0
    t_CL = 0
    E_list = []
    V_list = []
    
    for h in range_h:
        
        rho_i = ISA(h)[3]
        velo_som_i = ISA(h)[2]
        V_i = (jet.W / (0.5 * CL * rho_i * jet.S))**.5
        mach_i = V_i / velo_som_i
        
        CL_i = CL #mantendo o CL igual o CL inicial
        CD_i = polar(param, CL_i, mach_i)[0]
        E_i = CL_i / CD_i
                
        #Alcance
        deltaX_CL_i = E_i * range_dh # E * [h(i+1) - h(i)]
        deltaX_CL += deltaX_CL_i
        
        E_list.append(E_i)
        V_list.append(V_i)
        
        
        #Autonomia
        exp_t =  np.e**(- (h - range_dh)/(2*beta)) - np.e**(-h/(2*beta))
        t_CL_i = 2*beta*E_i * (((rho0 * CL_i)/(2 * jet.WL))**.5) * exp_t
        t_CL += t_CL_i
        
    if graph == True:
        
        fig_CL = plt.figure(figsize=(10,7))
        
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
        
        if save_graph == True:
            plt.savefig("alcance_autonomia_CL.pdf")
        
        plt.show()
        
    return deltaX_CL, t_CL









#%%% Melhot Alcance e melhor Autonomia



#%% Alcance e Autonomia com V constante

#%%%  Alcance e Autonomia

#Alcance
def alcance_autonomia_V(altitude, V, graph = False, save_graph = False):
    
    ## Considerando velocidade e altitude de cruzeiro
    h1 = 0 # [m]
    h2 = altitude # [m]

    velo_som = ISA(h2)[2] # [m/s]
    mach_cru = V / velo_som
    
    range_h = np.linspace(h2,h1,10000)
    range_dh = range_h[0] - range_h[1] # da um dh próximo de 13.11 m
    
    deltaX_V = 0
    t_V = 0
    E_list = []
    CL_list = []
    
    for h in range_h:
        
        rho_i = ISA(h)[3]
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
        
        fig_V= plt.figure(figsize=(10,7))
        
        
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
        plt.ylabel("Velocidade [m/s] ", fontsize = 12)
        ax = plt.gca()
        ax.invert_xaxis()
        
        if save_graph == True:
            plt.savefig("alcance_autonomia_V.pdf")
        
        plt.show()
    
    return deltaX_V, t_V

#%%% Melhor Alcance e Autonomia 
    