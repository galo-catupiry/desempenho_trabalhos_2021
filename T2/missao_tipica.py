"""
Analise de Missao Tipica - T2

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""

import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from aircraft import JetStar
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from ambiance import Atmosphere as Atmo
from Interpolacao import DragPolar

#%% dados de entrada
jet = JetStar(None)
h_cruz = 42500*0.3048 # ft
M = 20000 # kg
M_pouso = 13000 # kg
rho_SL = Atmo(0).density[0]
g = 9.81
h_cgh = 802 # m
h_bsb = 1066 # m
Mach_cruz = 0.85
V_cruz = Atmo(h_cruz).speed_of_sound[0]*Mach_cruz
dist_voo = 1000 # km
# potencia
n = 0.85
T0 = 64000
# arrasto
drag = DragPolar()
# CL
CL_max = 1.34

#%% razao de subida
def P_req(h,V):
    rho = Atmo(h).density[0]
    cl = 2*M*g/(jet.S*rho*V**2)
    drag.CLp = cl
    drag.Mp = V/Atmo(h).speed_of_sound[0]
    cd = drag.polar()
    return (0.5*rho*pow(V,2)*cd*jet.S)*V

def P_disp(h,V):
    T = T0*pow(Atmo(h).density[0]/rho_SL,n)
    return T*V

def h_dot(h,V):
    return (P_disp(h,V)-P_req(h,V))/(M*g)

def h_dot_max(h):
    V = fmin(lambda v: -h_dot(h,v),100,xtol=1e-8,disp=False)[0]
    return h_dot(h,V),V
    
## Trecho ascendente
h_asc = np.linspace(h_cgh,h_cruz,500)
dh = h_asc[1]-h_asc[0]
v_asc = []
x_asc = []

for i,h in enumerate(h_asc):
    if i == 0:
        x_i = 0
    else:
        x_i = x_asc[-1]
    hdot_i,v_i = h_dot_max(h)
    v_asc.append(v_i)
    dt = dh/hdot_i
    dx = v_i*dt/1000+x_i
    x_asc.append(dx)

# trecho descendente - iremos supor que a aeronave se aproxima com um angulo
# de 5 graus constante
dist_desc = (h_cruz - h_bsb)/np.tan(np.deg2rad(5))/1000 # em km
V_estol_pouso = np.sqrt(2*M_pouso*g/(Atmo(h_bsb).density[0]*jet.S*CL_max))

#%% plot
fig,[ax,ax2] = plt.subplots(2,1,sharex=True,figsize=(6,3.5))
ax.grid()
ax2.grid()
# altitude
ax.plot(x_asc,h_asc,'b')
ax.plot((x_asc[-1],dist_voo-dist_desc),(h_cruz,h_cruz),'b')
ax.plot((dist_voo-dist_desc,dist_voo),(h_cruz,h_bsb),'b')
ax.set_ylabel('Altitude [m]')

# velocidade
ax2.plot(x_asc,v_asc,'r')
ax2.plot((x_asc[-1],x_asc[-1]),(v_asc[-1],V_cruz),'r')
ax2.plot((x_asc[-1],dist_voo-dist_desc),(V_cruz,V_cruz),'r')
ax2.plot(((dist_voo-dist_desc,dist_voo)),(V_cruz,1.3*V_estol_pouso),'r')
ax2.set_ylabel('Velocidade [m/s]')
ax2.set_xlabel('Distancia percorrida [km]')

plt.tight_layout()
plt.savefig('missao_tipica.pdf')
















