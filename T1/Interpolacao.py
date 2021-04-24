# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 15:48:47 2021

"""

# Bibliotecas:
    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import leastsq
from mpl_toolkits import mplot3d

# Planilhas:

xlsx = pd.ExcelFile('C:\\Users\Administrador\Documents\GitHub\desempenho_trabalhos_2021\Planilhas\Dados - JetStar.xlsx')
CL_mach = pd.read_excel(xlsx, 'CL_mach')
CD_mach = pd.read_excel(xlsx, 'CD_mach')

# Dados experimentais:

CL_exp = CL_mach.iloc[:,1]
CD_exp = CD_mach.iloc[:,1]
M_exp  = CL_mach.iloc[:,0]
     
# Funções:

def func_residuo(p):
    
    R = (p[0] + p[1]*CL_exp + p[2]*CL_exp**2)*(p[3] + p[4]*M_exp + p[5]*M_exp**2) - CD_exp
    
    return R

def polar(x, CLp, Mp):
    
    CL = CLp
    M  = Mp
    resp = (x[0] + x[1]*CL + x[2]*CL**2)*(x[3] + x[4]*M + x[5]*M**2)
    
    return resp

param = 0.1*np.ones(6)
param = leastsq(func_residuo,param)
param = np.array(param[0])

CLp = np.linspace(0.1,1,20)
Mp  = np.linspace(0.1,0.9,20)

# Gráfico:
    
fig = plt.figure()
ax = plt.axes(projection="3d")
CLp, Mp = np.meshgrid(CLp, Mp)
resp = polar(param,CLp,Mp)
ax.plot_surface(CLp,Mp, resp, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.scatter3D(CL_exp,M_exp, CD_exp,c = M_exp ,cmap='flag',alpha=1)
ax.set_title("Polar de arrasto")
ax.set_xlabel("CL")
ax.set_ylabel("M")
ax.set_zlabel("CD")
plt.show()