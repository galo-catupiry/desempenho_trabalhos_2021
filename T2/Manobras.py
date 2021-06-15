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
    omega_list = []
    
    for j in fc:
        omega = (g/V)*np.sqrt(j**2 - 1)
        omega_list.append(omega)
    return omega_list

def R(V, omega, fc):
    ''' 
    Raio da curva coordenada.
    '''
    R_list = []
    
    for i in np.arange(0, len(fc)):
        R = V/omega[i] 
        R_list.append(R)
    return R_list

# ============================================= 
# Gráficos

def TD_vs_V_manobras(h, fc, V, D_manobra, T_manobra):
    
    plt.figure()
    plt.style.use('default')
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("T and D [N]")
    plt.grid(True)
    color=iter(plt.cm.rainbow(np.linspace(0,1,len(fc))))
    
    for i in np.arange(0, len(fc)):
        c = next(color)
        plt.plot(V, [T_manobra[i]]*len(V), 'r')
        plt.plot(V, D_manobra[i], color = c, label = 'n = {:.2f}'.format(fc[i]))
    
    plt.legend(loc = 'best', framealpha = 1)
    plt.ylim(top = 40000)
    return

def omega_vs_V(V, omega, fc, R):
    
    plt.figure()
    plt.style.use('seaborn-bright')
    plt.xlabel('V [m/s]')
    plt.ylabel('$\\omega$ [rad/s]')
    plt.grid(False)
    color=iter(plt.cm.jet(np.linspace(0,1,len(fc))))
    offset_x = 5
    offset_y = 0
    
    for i in np.arange(0, len(fc)):
        c = next(color)
        plt.plot(V, omega[i], color = c, label = 'n = {:.2f}'.format(fc[i]))
        
    for j in np.arange(0, len(R)):
        aux = V/R[j]
        plt.plot(V, aux, '-.k', linewidth = 0.9)
        plt.text(V[-1] + offset_x, aux[-1] + offset_y, 'R = {:.1f} [km]'.format(R[j]/1000))
    
    plt.legend(loc = 'best', framealpha = 1)
    
    plt.xlim(right = 400)
    
    return

