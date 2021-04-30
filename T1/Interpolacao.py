"""
T1 - SAA0183
Grupo E
"""

# Bibliotecas:
import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import leastsq
from mpl_toolkits import mplot3d

# Planilhas:

xlsx = pd.ExcelFile(parent_dir + '\Planilhas\Dados - JetStar.xlsx')
CL_mach = pd.read_excel(xlsx, 'CL_mach')
CD_mach = pd.read_excel(xlsx, 'CD_mach')

# Dados experimentais:

CL_exp = CL_mach.iloc[:,1]
CD_exp = CD_mach.iloc[:,1]
M_exp  = CL_mach.iloc[:,0]
     
# Funções:

def func_residuo(p):
    
    R = (p[0] + p[1]*CL_exp**2)*(p[2] + p[3]*M_exp + p[4]*M_exp**2) - CD_exp
    
    return R

def polar(x, CLp, Mp):
    
    CL = CLp
    M  = Mp
    resp = (x[0] + x[1]*CL**2)*(x[2] + x[3]*M + x[4]*M**2)
    
    CD_0 = x[0]*(x[2] + x[3]*M + x[4]*M**2)
    K = x[1]*(x[2] + x[3]*M + x[4]*M**2)
    
    return resp, CD_0,K

param = 0.1*np.ones(5)
param = leastsq(func_residuo,param)
param = np.array(param[0])


# Dados de entrada:
    
CLp = np.linspace(0.05,1,20)
Mp  = np.linspace(0.05,1,20)

[resp, CD_0,K] = polar(param,CL_exp,M_exp)

print("     CD_exp - CD_previsto \n {}".format(abs(CD_exp-resp)))

# ----------------------------- GRÁFICOS ------------------------------------#

# Superfície 3D: Mach variável
'''
fig = plt.figure()
ax = plt.axes(projection="3d")
CLp, Mp = np.meshgrid(CLp, Mp)
resp_graf = polar(param,CLp,Mp)[0]
ax.plot_surface(Mp,CLp, resp_graf, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none', alpha = 0.9)
ax.scatter3D(M_exp,CL_exp, CD_exp,c = 'red',alpha=1)
#ax.set_title("Polar de arrasto")
ax.set_xlabel("CL")
ax.set_ylabel("M")
ax.set_zlabel("CD")
plt.savefig("polar_arrasto.svg")
plt.show()
'''

# Curva CL x CD: Mach fixo
''' 
fig = plt.figure()
resp_2d = polar(param,CLp,Mp)[0]
plt.plot(CLp,resp_2d)
#plt.legend(loc = 'best', framealpha = 1)
plt.xlabel("CL")
plt.ylabel("CD")
plt.grid(True)
plt.title("Curva CL x CD")
plt.show()
'''