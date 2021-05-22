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
from scipy.optimize import fsolve

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

        
        drag = DragPolar()
        
        for h in range_h:
            

            rho_i = (Atmosphere(h).density[0] + Atmosphere(h-range_dh).density[0])/2
            velo_som_i = (Atmosphere(h).speed_of_sound[0] + Atmosphere(h-range_dh).speed_of_sound[0])/2

            V_i = (jet.W / (0.5 * CL * rho_i * jet.S))**.5
                    
            mach_i = V_i / velo_som_i
            drag.Mp = mach_i
            
            CL_i = CL #mantendo o CL igual o CL inicial
            drag.CLp = CL_i
            CD_i = drag.polar()
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
            
          
        print("Altitude : {} m".format(altitude))
        print("Velocidade inicial {} m/s".format(round(V,2)))
        print("CL : {} ".format(round(CL,2)))
        print("Alcance: {} m ".format(round(deltaX_CL,2)))
        print("Autonomia: {} s".format(round(t_CL,2)))
        print("---------------------------------------------------")  
        
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

    
        h_linspace = np.linspace(h2,h1,800, retstep = True) 
        range_h = h_linspace[0]
        range_dh = abs(h_linspace[1])
    
        deltaX_V = 0
        t_V = 0
        E_list = []
        CL_list = []
    
        drag = DragPolar()
    
        for h in range_h:
                    
            rho_i = (Atmosphere(h).density[0] + Atmosphere(h-range_dh).density[0])/2
            velo_som_i = (Atmosphere(h).speed_of_sound[0] + Atmosphere(h-range_dh).speed_of_sound[0])/2
            
            
            mach_i = V/velo_som_i
            drag.Mp = mach_i
        
            CL_i = (2 * jet.WL)/(rho_i * V**2)
            drag.CLp = CL_i
            CD_i = drag.polar()
            E_i = CL_i / CD_i
        
        
            A_i = (rho0 * drag.CD0 * V**2)/(2 * jet.WL)
            B_i = (2 * drag.K * jet.WL)/(rho0 * V**2)
        
            tan1_i = np.arctan((A_i/B_i) * np.e**(-(h - range_dh)/beta))
            tan2_i = np.arctan((A_i/B_i) * np.e**(-h/beta))
        
            #Alcance
            deltaX_V_i = (beta/B_i)* (tan1_i - tan2_i)
            deltaX_V += deltaX_V_i
            
            #Autonomia
            t_V_i = beta/(B_i * V) * (tan1_i - tan2_i)
            t_V += t_V_i
            
            E_list.append(E_i)
            CL_list.append(CL_i)
            
            
            
        print("Altitude : {} m".format(altitude))
        print("Velocidade {} m/s".format(round(V,2)))
        print("Alcance: {} m".format(round(deltaX_V,2)))
        print("Autonomia: {} s".format(round(t_V,2)))
        print("---------------------------------------------------")  
        
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
        
    V_list = np.linspace(0.4*V_stall, 1.5*V_stall, 100)
    
    
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
    
#%% Parâmetros ótimos para CL constante

drag_CL = DragPolar()
drag_CL.CLp = 0


def max_CL_cte(x,h1,h2,cond):
    
    rho = (Atmosphere(h1).density[0] + Atmosphere(h2).density[0])/2
    velo_som = (Atmosphere(h1).speed_of_sound[0] +  Atmosphere(h2).speed_of_sound[0])/2
    
    # Variáveis desconhecidas:
    CD0 = x[0]
    k = x[1]
    V = x[2]
    
    drag_CL.Mp = V/velo_som

    if (cond.get('condicao') == 'max_range'):
        
        # CL de otimização do alcance:
        CL = np.sqrt(CD0/k)
            
        # Equações do sistema:
        eq1 = (jet.W / (0.5 * CL * rho * jet.S))**.5 - V
        eq2 = drag_CL.CD0 - CD0
        eq3 = drag_CL.K - k
            
        return [eq1,eq2,eq3] 

    elif (cond.get('condicao') == 'max_endurance'):
        
        # CL de otimização da autonomia:
        CL = np.sqrt(3*CD0/k)
        
        # Equações do sistema:
        eq1 = (jet.W / (0.5 * CL * rho * jet.S))**.5 - V
        eq2 = drag_CL.CD0 - CD0
        eq3 = drag_CL.K - k
        
        return [eq1,eq2,eq3]


def parametros_otimos_CL(cond):
    
    velocidade = cond.get('v')
    altitude = cond.get('h')
    
    xmax_CL = []
    CL_resp1 = []
    tmax_CL = []
    CL_resp2 = []
    E_resp = []
    
    h1 = 0 # [m]
    h2 = altitude # [m]
    
    h_linspace = np.linspace(h2,h1,200, retstep = True)
    range_h = h_linspace[0]
    
    initial = [0.08, 0.01, velocidade]
    
    if (cond.get('condicao') == 'max_range'):                       # Case: Máximo range
        
        for i in np.arange(0,len(range_h) - 1):
            h1 = range_h[i+1]
            h2 = range_h[i]
        
            [CD0_resp,k_resp,V_resp] = fsolve(max_CL_cte, initial, args = (h1,h2,cond))
                
            initial = [CD0_resp,k_resp,V_resp]
              
            # Parâmetros de interesse:   
            
            CL = np.sqrt(CD0_resp/k_resp)
            CD = 2*CD0_resp
            E = CL/CD
            range_max = E*(h2 - h1)  
            
            CL_resp1.append(CL)
            xmax_CL.append(range_max)
    
    
        plt.figure()
        plt.plot(range_h[:len(range_h) - 1], CL_resp1)
        plt.gca().invert_xaxis()
        plt.ylabel(" CL para máximo alcance")
        plt.xlabel ("Altitude [m]")
        plt.grid(True)
            
        return xmax_CL
    
    elif (cond.get('condicao') == 'max_endurance'):
        
        for i in np.arange(0,len(range_h) - 1):
            h1 = range_h[i+1]
            h2 = range_h[i]
        
            [CD0_resp,k_resp,V_resp] = fsolve(max_CL_cte, initial, args = (h1,h2,cond))
                
            initial = [CD0_resp,k_resp,V_resp]
              
            CL = np.sqrt(3*CD0_resp/k_resp)
            CD = 4*CD0_resp
            E = CL/CD
            endurance_max = (2*beta*E)*np.sqrt((rho0*jet.S*CL)/(2*jet.W))*(np.exp(-h1/(2*beta)) - np.exp(-h2/(2*beta)))
            
            CL_resp2.append(CL)
            E_resp.append(E)
            tmax_CL.append(endurance_max)
        
        plt.figure()
        #plt.plot(range_h[:len(range_h) - 1], CL_resp2)
        plt.plot(range_h[:len(range_h) - 1], E_resp)
        plt.gca().invert_xaxis()
        plt.ylabel(" Eficiência máxima (E)")
        #plt.ylabel(" CL para máxima autonomia")
        plt.xlabel ("Altitude [m]")
        plt.grid(True)
        
        return tmax_CL

#%% Parâmetros ótimos para V constante 

drag_V = DragPolar()
drag_V.CLp = 0


def max_V_cte(x,h1, h2, cond):
    
    CD0 = x[0]
    k = x[1]
    V = x[2]
    
    # Parâmetros comuns à ambos os casos:
        
    a1 = rho0**2*CD0*jet.S**2/(4*jet.W**2*k)*np.exp(-h1/beta)
    a2 = rho0**2*CD0*jet.S**2/(4*jet.W**2*k)*np.exp(-h2/beta)
    
    #A = rho0*CD0*(V**2)*jet.S/(2*jet.W) 
    #B = 2*jet.W*k/(rho0*(V**2)*jet.S)
    
    velo_som1 = Atmosphere(h1).speed_of_sound[0]
    velo_som2 = Atmosphere(h2).speed_of_sound[0]
    velo_som = (velo_som1 + velo_som2)/2
    
    drag_V.Mp = V/velo_som
    
    CD0_exp = drag_V.CD0 - CD0
    k_exp = drag_V.K - k
    
    # Velocidade para cada caso de análise:
        
    if (cond.get('condicao') == 'max_range'):
        V_exp = (3/(a1*a2))**(1/8) - V
    elif(cond.get('condicao') == 'max_endurance'):
        V_exp = (5/(3*a1*a2))**(1/8) - V

    return [V_exp, CD0_exp, k_exp]

def parametros_otimos_V(cond):
    
    xmax = []
    Vresp1 = []
    tmax = []
    Vresp2 = []
    
    
    altitude = cond.get('h')
    velocidade = cond.get('v')
    
    
    h1 = 0 # [m]
    h2 = altitude # [m]
    initial = [0.08, 0.01, velocidade]
    
    h_linspace = np.linspace(h2,h1,200, retstep = True,endpoint = True) 
    range_h = h_linspace[0]
    
    if (cond.get('condicao') == 'max_range'):                       # Case: Máximo range
        
        for i in np.arange(0,len(range_h) - 1):
            h1 = range_h[i+1]
            h2 = range_h[i]
            
            [CD0_resp,k_resp,V_resp] = fsolve(max_V_cte, initial, args = (h1,h2,cond))
            
            initial = [CD0_resp,k_resp,V_resp]
            
            a1 = rho0**2*CD0_resp*jet.S**2/(4*jet.W**2*k_resp)*np.exp(-h1/beta)
            a2 = rho0**2*CD0_resp*jet.S**2/(4*jet.W**2*k_resp)*np.exp(-h2/beta)
            V = (3/(a1*a2))**(1/8)
            range_max = beta*rho0*jet.S/(2*jet.W*k_resp)*(a1 - a2)*V**6/(1 + a1*a2*V**8)
    
            xmax.append(range_max)
            Vresp1.append(V)
 
        
        plt.figure()
        plt.plot(Vresp1,range_h[:len(range_h) - 1])
        plt.xlabel(" Velocidade para máximo alcance [m/s]")
        plt.ylabel ("Altitude [m]")
        plt.grid(True)
        
        return xmax 
    
    elif (cond.get('condicao') == 'max_endurance'):                 # Case: Máximo endurance
         
        for j in np.arange(0,len(range_h) - 1):
            h1 = range_h[j+1]
            h2 = range_h[j]
            
            [CD0_resp,k_resp,V_resp] = fsolve(max_V_cte, initial, args = (h1,h2,cond))
            
            initial = [CD0_resp,k_resp,V_resp]
            
            a1 = (rho0**2)*CD0_resp*jet.S**2/(4*jet.W**2*k_resp)*np.exp(-h1/beta)
            a2 = (rho0**2)*CD0_resp*jet.S**2/(4*jet.W**2*k_resp)*np.exp(-h2/beta)
            
            V = (5/(3*a1*a2))**(1/8)
            
            A = rho0*CD0_resp*(V**2)*jet.S/(2*jet.W) 
            B = 2*jet.W*k_resp/(rho0*(V**2)*jet.S)
    
            endurance_max = beta/(B*V_resp)*(np.arctan(A/B*np.exp(-h1/beta)) - np.arctan(A/B*np.exp(-h2/beta)))
            
            tmax.append(endurance_max)
            Vresp2.append(V)
            
        plt.figure()
        plt.plot(Vresp2,range_h[:len(range_h) - 1], 'r')
        plt.xlabel(" Velocidade para máxima autonomia [m/s]")
        plt.ylabel ("Altitude [m]")
        plt.grid(True)
        
        return tmax  