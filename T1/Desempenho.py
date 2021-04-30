"""
T1 - SAA0183
Grupo E
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
def alcance_autonomia_CL(graphs = 0):
    
    ## Considerando velocidade e altitude de cruzeiro
    h1 = 0 # [m]
    h2 = 13105 # [m]
    
    velo_cru = 811 / 3.6 # [m/s]
    densidade_cru = ISA(h2)[3] #[kg/m^3]
    velo_som_cru = ISA(h2)[2] # [m/s]
    mach_cru = velo_cru / velo_som_cru
    
    CL_cru = jet.W / (0.5 * jet.S * (velo_cru**2) * densidade_cru) #esse cl é mantido constante
    #print("CL cruzeiro: {}".format(round(CL_cru,2)))
    
    
    range_h = np.linspace(h2,h1,10000)
    range_dh = range_h[0] - range_h[1] # da um dh próximo de 13.11 m
    
    deltaX_CL = 0
    t_CL = 0
    E_list = []
    V_list = []
    
    
    for h in range_h:
        
        
        rho_i = ISA(h)[3]
        velo_som_i = ISA(h)[2]
        V_i = (jet.W / (0.5 * CL_cru * rho_i * jet.S))**.5
        mach_i = V_i / velo_som_i
        
        CL_i = CL_cru #mantendo o CL igual o CL de cruzeiro
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
        
    if graphs == 1:
        
        
        fig_CL = plt.figure(figsize=(10,7))
        
        #Eficiencia pela altitude
        fig_E = fig_CL.add_subplot(121)
        plt.plot(range_h, E_list, label = "E")
        #plt.legend()
        plt.xlabel("Altitude [m]", fontsize = 12)
        plt.ylabel("Eficiência ", fontsize = 12)
        plt.grid()
        ax = plt.gca()
        ax.invert_xaxis()
        
        #Velocidade pela altitude
        fig_V = fig_CL.add_subplot(122)
        plt.plot(range_h, V_list, label = "V")
        #plt.legend()
        plt.grid()
        plt.xlabel("Altitude [m]", fontsize = 12)
        plt.ylabel("Velocidade [m/s] ", fontsize = 12)
        ax = plt.gca()
        ax.invert_xaxis()
        
        plt.show()
        
    

    return deltaX_CL, t_CL





# print("----- Alcance e autonomia [Caso CL constante] -----" )
# #print("h2: 13105, h1: 0 [m]")
# print("Delta_X = {} [m]".format(round(alcance_autonomia_CL()[0],2)))
# print("t = {} [s]".format(round(alcance_autonomia_CL()[1],2)))
# alcance_autonomia_CL(1)



#%%% Melhot Alcance e melhor Autonomia




    
    
















#%% Alcance e Autonomia com V constante

#%%%  Alcance e Autonomia

#Alcance
def alcance_autonomia_V(graphs = 0):
    
    ## Considerando velocidade e altitude de cruzeiro
    h1 = 0 # [m]
    h2 = 13105 # [m]
    
    velo_cru = 811 / 3.6 # [m/s]
    densidade_cru = ISA(h2)[3] #[kg/m^3]
    velo_som_cru = ISA(h2)[2] # [m/s]
    mach_cru = velo_cru / velo_som_cru
    
    range_h = np.linspace(h2,h1,10000)
    range_dh = range_h[0] - range_h[1] # da um dh próximo de 13.11 m
    
    deltaX_V = 0
    t_V = 0
    E_list = []
    CL_list = []
    
    
    for h in range_h:
        
        rho_i = ISA(h)[3]
        CL_i = (2 * jet.WL)/(rho_i * velo_cru**2)
        CD_i = polar(param, CL_i, mach_cru)[0]
        E_i = CL_i / CD_i
        CD0_i = polar(param, CL_i, mach_cru)[1]
        K_i = polar(param, CL_i, mach_cru)[2]
        
        A_i = (rho0 * CD0_i * velo_cru**2)/(2 * jet.WL)
        B_i = (2 * K_i * jet.WL)/(rho0 * velo_cru**2)
        
        tan1_i = np.arctan(B_i**-1 * A_i * np.e**(-(h - range_dh)/beta))
        tan2_i = np.arctan(B_i**-1 * A_i * np.e**(-h/beta))
        
        deltaX_V_i = (beta/B_i)* (tan1_i - tan2_i)
        deltaX_V += deltaX_V_i
        
        E_list.append(E_i)
        CL_list.append(CL_i)
        
        
    if graphs == 1:
        
        fig_V = plt.figure(figsize=(10,7))
        
        #Eficiencia pela altitude
        fig_E = fig_V.add_subplot(121)
        plt.plot(range_h, E_list, label = "E")
        #plt.legend()
        plt.xlabel("Altitude [m]", fontsize = 12)
        plt.ylabel("Eficiência ", fontsize = 12)
        plt.grid()
        ax = plt.gca()
        ax.invert_xaxis()
        
        #CL pela altitude
        fig_CL = fig_V.add_subplot(122)
        plt.plot(range_h, CL_list, label = "CL")
        #plt.legend()
        plt.grid()
        plt.xlabel("Altitude [m]", fontsize = 12)
        plt.ylabel("Velocidade [m/s] ", fontsize = 12)
        ax = plt.gca()
        ax.invert_xaxis()
        
        plt.show()
    
    
    
    

    # #formula 7.20
    # rho_ssl = 0
    # CD0 = 0
    # K = 0
    # A = (rho_ssl * CD0 * V**2)/(2 * jet.WL)
    # B = (2 * K * jet.WL)/(rho_ssl * V**2)
    # tan1 = np.arctan(B**-1 * A * np.e**(-h1/beta))
    # tan2 = np.arctan(B**-1 * A * np.e**(-h2/beta))
    
    # deltaX_V = (beta / B) * (tan1 - tan2)
    
    return deltaX_V


print("----- Alcance e autonomia [Caso V constante] -----" )
#print("h2: 13105, h1: 0 [m]")
print("Delta_X = {} [m]".format(round(alcance_autonomia_V(),2)))
#print("t = {} [s]".format(round(alcance_autonomia_V(),2)))
alcance_autonomia_V(1)

#Autonomia
# def endurance_V(h1, h2, V):
#     rho_ssl = 1
#     V = 1
#     CD0 = 0
#     K = 0
    
#     a = (CD0 * rho_ssl**2 * np.e**(-h1/beta))/(4*K * jet.WL**2)
#     b = (CD0 * rho_ssl**2 * np.e**(-h2/beta))/(4*K * jet.WL**2)
    
#     arc_tan = np.arctan(((a-b)*(V**4)) / (1 + a*b*(V**8)))
                        
#     tv = ((beta*rho_ssl* V)/(2*K * jet.WL)) * arc_tan
    
#     return tv
    
    
    
#%%% Melhor Alcance e Autonomia 
    