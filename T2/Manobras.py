"""
T2 - SAA0183 - Analise em Manobras

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
from scipy.optimize import fsolve

from aircraft import JetStar
from ambiance import Atmosphere
from Interpolacao import DragPolar
from Cruzeiro import jet_buoyancy

# =============================================
jet = JetStar(1)
sealevel = Atmosphere(0)
rho0 = sealevel.density[0]
g = 9.81
# =============================================

# Polar de arrasto
drag_manobra = DragPolar()

def CL(fc, V, h):
    '''
    Coeficiente de sustentação em curva coordenada.

    Entradas:

    h: Altitude (inteiro)
    V: Velocidade (lista)
    
    '''
    rho = Atmosphere(h).density[0]
    CL = 2*fc*jet.W/(rho*(V**2)*jet.S)
    return CL

def drag(V, h, fc):
    ''' 
    Arrasto em curva coordenada.
    '''
    
    D_list = []
    rho = Atmosphere(h).density[0]
    drag_manobra.Mp = V/Atmosphere(h).speed_of_sound[0] 
    
    for i in fc:
        drag_manobra.CLp = CL(i, V, h)
        CD = drag_manobra.polar()
<<<<<<< HEAD
=======
        
>>>>>>> Abner
        D = (1/2)*rho*(V**2)*jet.S*CD
        D_list.append(D)
        
    return D_list

def T(h, T0, n, fc):
    '''
    Empuxo em curva coordenada.
    '''
    
    T = []
    
    for i in fc:
        sigma = Atmosphere(h).density[0]/rho0
        Ti = T0*(sigma**n) 
        T.append(Ti)
    
    return T

def omega(fc, V):
    '''
    Velocidade angular (cte) em curva coordenada.
    '''
    omega_resp = []

    for i in np.arange(len(fc)):
        omega = (g/V[i])*np.sqrt(fc[i]**2 - 1)
        omega_resp.append(omega)

    return omega_resp 

def R(V, omega, fc):
    ''' 
    Raio da curva coordenada.
    '''
    R_list = []
    
    for i in np.arange(0, len(fc)):
        R = V/omega[i] 
        R_list.append(R)
    return R_list

def fc_maximo(fc, V, h, T0, n):
        
        T_aux = T(h, T0, n, fc)
        D_aux = drag(V, h, fc)
        
        i = 0
        while(max([T_aux[i]]*len(V) - D_aux[i]) > 0):
            i += 1   
        
        fc_max = fc[i + 1]
        
        return fc_max

def velocidades_manobra_eq(x, h, T0, n , fc):
    ''' 
    Montagem do sistema 3x3 para determinação
    das velocidades possíveis de manobra 
    '''

    CD0 = x[0]
    K = x[1]
    V = x[2]

    Empuxo = T(h, T0, n, [fc])
    Arrasto = drag(V, h, [fc])

    # Equations to be solved
    CD0_exp = drag_manobra.CD0 - CD0
    K_exp = drag_manobra.K - K
    V_exp = Empuxo[0] - Arrasto[0]

    eq = [CD0_exp, K_exp, V_exp]

    return eq

def velocidades_manobra_solver(V1, V2, h, T0, n , fc):
    '''
    Resolução do sistema 3x3 para
    as velocidades de manobra '''

    V1_manobras, V2_manobras = [], []

    V1_0 = V1
    V2_0 = V2

    CD0_0 = 0.01
    K_0 = 0.01

    for i in fc:

        [CD0_resp, K_resp, V1_man] = fsolve(velocidades_manobra_eq, (CD0_0,K_0,V1_0), args=(h, T0, n , i))
        [CD0_resp, K_resp, V2_man] = fsolve(velocidades_manobra_eq, (CD0_0,K_0,V2_0), args=(h, T0, n , i))

        V1_manobras.append(V1_man)
        V2_manobras.append(V2_man)

        V1_0 = V1_man 
        V2_0 = V2_man 


    return V1_manobras, V2_manobras

def estol(fc_s, W, S, h, CLmax):

    sigma = Atmosphere(h).density[0]/Atmosphere(0).density[0]
    omega_resp, V_resp = [],[]

    for i in fc_s:
        V_s = (2*i*(W/S)/(Atmosphere(0).density[0]*sigma*CLmax))**(0.5)
        omega_s = 9.81/V_s*np.sqrt(i**2 - 1)

        V_resp.append(V_s)
        omega_resp.append(omega_s)

    return V_resp, omega_resp

# ============================================= 
# Gráficos

def TD_vs_V_manobras(h, fc, V, D_manobra, T_manobra):
    
    plt.figure()
    plt.style.use('default')
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("T and D [N]")
    plt.grid(False)
    color=iter(plt.cm.rainbow(np.linspace(0,1,len(fc))))
    
    for i in np.arange(0, len(fc)):
        c = next(color)
        plt.plot(V, [T_manobra[i]]*len(V), 'r')
        plt.plot(V, D_manobra[i], color = c, label = 'n = {:.2f}'.format(fc[i]))
    
    plt.legend(loc = 'best', framealpha = 1)
    plt.text(200, 35000,"h = {:.1f} km".format(h[0]/1000))
    plt.ylim(bottom = 0,top = 40000)
    plt.savefig('TD_vs_V_manobras.svg')
    plt.show()
    return

def omega_vs_V(V1_manobras, V2_manobras, omega1_manobras, omega2_manobras, 
               fc, h, V, omega_V, fc_max, Vs, omega_s, Vd):
    
    plt.figure()
    plt.style.use('seaborn-bright')
    plt.xlabel('V [m/s]')
    plt.ylabel('$\\omega$ [rad/s]')
    plt.grid(False)
    #color=iter(plt.cm.jet(np.linspace(0,1,len(fc))))

    plt.plot(V1_manobras, omega1_manobras, 'k', label = 'h = {:.1f} km'.format(h[0]/1000))
    plt.plot(V2_manobras, omega2_manobras, 'k')
    plt.plot(V, omega_V, 'b', label = '$n_{{{}}} = {:.2f}$'.format("max", fc_max))
    plt.plot(Vs, omega_s, 'r', label = 'Estol')
    #plt.vlines(Vd, 0 , 0.15, label = 'Aeroelasticidade')
    plt.ylim(bottom = 0, top = 0.5)
  
    plt.legend(loc = 'best', framealpha = 1)
    plt.savefig("desempenho_em_curva.svg")
    
    plt.show()
    return

