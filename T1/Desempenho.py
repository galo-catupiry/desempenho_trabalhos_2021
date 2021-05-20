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
from Interpolacao import DragPolar
from aircraft import JetStar
from datetime import datetime

#%% Dados Gerais

jet = JetStar(1)

sealevel = Atmosphere(0)
beta = 9296
rho0 = sealevel.density[0]

#%% Alcance e Autonomia com CL constante
def alcance_autonomia_CL(V, graph_E_V = False, save_graph_E_V = False, *altitudes):
    
    print("----- Alcance e Autonomia [Caso CL constante] -----")
    
    if graph_E_V == True:
        
        fig_E ,ax_E = plt.subplots(figsize=(6,4))
        ax_E.set_xlabel("Altitude [m]", fontsize = 12)
        ax_E.set_ylabel("Eficiência ", fontsize = 12)
        ax_E.grid()
        ax = plt.gca()
        ax.invert_xaxis()
        
        
        
        fig_V, ax_V = plt.subplots(figsize=(6,4))
        ax_V.grid()
        ax_V.set_xlabel("Altitude [m]", fontsize = 12)
        ax_V.set_ylabel("Velocidade [m/s] ", fontsize = 12)
        ax = plt.gca()
        ax.invert_xaxis()
        
    for altitude in altitudes:
    
        h1 = 0 # [m]
        h2 = altitude # [m]
        
        rho = Atmosphere(h2).density[0]
        
        CL = jet.W / (0.5 * jet.S * (V**2) * rho) #esse cl é mantido constante
        
        h_linspace = np.linspace(h2,h1,800, retstep = True)
        range_h = h_linspace[0]
        range_dh = abs(h_linspace[1])
        
        deltaX_CL = 0
        t_CL = 0
        E_list = []
        V_list = []
        hdot_list = []
        
        drag = DragPolar()
        
        for h in range_h:
            
            rho_i = Atmosphere(h).density[0]
            velo_som_i = Atmosphere(h).speed_of_sound[0]
            V_i = (jet.W / (0.5 * CL * rho_i * jet.S))**.5
                    
            mach_i = V_i / velo_som_i
            drag.Mp = mach_i
            
            CL_i = CL #mantendo o CL igual o CL inicial
            drag.CLp = CL_i
            CD_i = drag.polar()
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
            
          
        print("Altitude : {} m".format(altitude))
        print("Velocidade inicial {} m/s".format(velocidade))
        print("Alcance: {} m ".format(deltaX_CL))
        print("Autonomia: {} s".format(t_CL))
        print("----------------------")  
        
        try:
            ax_E.plot(range_h, E_list, label = "Altitude :{} m".format(altitude))
            ax_V.plot(range_h, V_list, label = "Altitude :{} m".format(altitude))
        except:
            None
            
    try:
        ax_E.legend()
        ax_V.legend()
        plt.tight_layout()
        
        if save_graph_E_V == True:
            fig_E.savefig("eficiencia_cl_constante.pdf")
            fig_V.savefig("velocidade_cl_constante.pdf")
    except:
        None
            



#%% Alcance e Autonomia com V constante
def alcance_autonomia_V(V, graph = False, save_graph = False, *altitudes):
    
    
    print("----- Alcance e Autonomia [Caso V constante] -----")
    
    if graph == True:
        
        fig_E ,ax_E = plt.subplots(figsize=(6,4))
        ax_E.set_xlabel("Altitude [m]", fontsize = 12)
        ax_E.set_ylabel("Eficiência ", fontsize = 12)
        ax_E.grid()
        ax = plt.gca()
        ax.invert_xaxis()
        
    
        fig_CL, ax_CL = plt.subplots(figsize=(6,4))
        ax_CL.grid()
        ax_CL.set_xlabel("Altitude [m]", fontsize = 12)
        ax_CL.set_ylabel("$C_L$", fontsize = 12)
        ax = plt.gca()
        ax.invert_xaxis()
    
    for altitude in altitudes:
    
        ## Considerando velocidade e altitude de cruzeiro
        h1 = 0 # [m]
        h2 = altitude # [m]

        velo_som = Atmosphere(h2).speed_of_sound[0]
        mach_cru = V / velo_som
    
        h_linspace = np.linspace(h2,h1,800, retstep = True) 
        range_h = h_linspace[0]
        range_dh = abs(h_linspace[1])
    
        deltaX_V = 0
        t_V = 0
        E_list = []
        CL_list = []
    
        drag = DragPolar()
    
        for h in range_h:
        
            rho_i = Atmosphere(h).density[0]
            velo_som_i = Atmosphere(h).speed_of_sound[0]
            mach_i = V/velo_som_i
            drag.Mp = mach_i
        
            CL_i = (2 * jet.WL)/(rho_i * V**2)
            drag.CLp = CL_i
            CD_i = drag.polar()
            E_i = CL_i / CD_i
        
        
            A_i = (rho0 * drag.CD0 * V**2)/(2 * jet.WL)
            B_i = (2 * drag.K * jet.WL)/(rho0 * V**2)
        
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
            
            
            
        print("Altitude : {} m".format(altitude))
        print("Velocidade inicial {} m/s".format(velocidade))
        print("Alcance: {} m ".format(deltaX_V))
        print("Autonomia: {} s".format(t_V))
        print("----------------------")  
        
        try:
            ax_E.plot(range_h, E_list, label = "Altitude :{} m".format(altitude))
            ax_CL.plot(range_h, CL_list, label = "Altitude :{} m".format(altitude))
        except:
            None
            
    try:
        ax_E.legend()
        ax_CL.legend()
        plt.tight_layout()
        
        if save_graph == True:
            fig_E.savefig("eficiencia_v_constante.pdf")
            fig_CL.savefig("cl_v_constante.pdf")
    except:
        None
#%% Gráfico de razão de descida por velocidade
def hdot_V(velocidade, save_graph = False, *altitudes):
    
    CL_max = 1
    
    V_stall = np.sqrt(jet.W / (CL_max * 0.5 * jet.S * Atmosphere(13105).density[0]))
    
    fig_hdotV = plt.figure()
    
    # E_max = 13.25
    # gamma_min = -1/E_max 
        
    # rho_cru = Atmosphere(13105).density[0]
    # V_cru = 811 / 3.6
    # mach_cru = V_cru / Atmosphere(13105).speed_of_sound[0]
    # CL_cru = jet.W / (0.5 * rho_cru * (V_cru**2.0) * jet.S)
    
    # drag_cru = DragPolar(CL_cru, mach_cru)
    # drag_cru.CLp = CL_cru
    # drag_cru.Mp = mach_cru
    
    V_list = np.linspace(0.4*V_stall, 1.5*V_stall, 100)
    
    #plt.plot(V_list, [i*gamma_min for i in V_list], color = 'k' , label = "$\gamma$ min= {}".format(round(gamma_min,5)))
    
    for altitude in altitudes:
        
        hdot_list = []
        
        rho = Atmosphere(altitude).density[0]
        velo_som = Atmosphere(altitude).speed_of_sound[0]
    
        drag = DragPolar()
        
        for Vi in V_list:
            
            mach = Vi/velo_som
            
            CL = jet.W / (0.5 * rho * (Vi ** 2) * jet.S)
            
            drag.CLp = CL
            drag.Mp = mach
            drag.polar(False)
        
            CD0 = drag.CD0
            K = drag.K
            
            CD = CD0 + K * (CL**2)

            E = CL/CD
            gamma = -1/E
            hdot = Vi*gamma
            
            hdot_list.append(hdot)
            
        plt.plot(V_list, hdot_list, label = "Altitude: {} m".format(altitude))
        
    plt.axvline(V_stall, ymin = -60, ymax = 1, ls = '--', color = 'k', label = "Velocidade de Stall")
    plt.grid()
    plt.ylim([-30,1])
    plt.xlim([100,250])
    plt.xlabel("Velocidade [m/s]", fontsize = 12)
    plt.ylabel("Razão de descida $\dot{h}$", fontsize = 12)
    plt.legend()
    plt.tight_layout()
    
    if save_graph == True:
        plt.savefig("hdot_V.pdf")
    
    plt.show()


#%% MAIN
#%% ------ MAIN -------
start = datetime.now()

altitude = 13105 # [m]
velocidade = 811/3.6  # [m/s]


alcance_autonomia_CL(velocidade, False, False, 
                                        altitude, altitude - 3000, altitude - 6000)

alcance_autonomia_V(velocidade, False, True,
                    altitude, altitude - 3000, altitude - 6000)


hdot_V(velocidade, False, altitude, altitude - 3000, 
                      altitude - 6000, altitude - 8000)

print(datetime.now() - start)
