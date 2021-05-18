"""
Cálculo dos parâmetros aerodinâmicos

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
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

xlsx = pd.ExcelFile(parent_dir + '/Planilhas/Dados-JetStar.xlsx')
CL_mach = pd.read_excel(xlsx, 'CL_mach')
CD_mach = pd.read_excel(xlsx, 'CD_mach')

# Dados experimentais:

CL_exp = CL_mach.iloc[:,1]
CD_exp = CD_mach.iloc[:,1]
M_exp  = CL_mach.iloc[:,0]
     
# Funções:

class drag_polar():
    
    def __init__(self,CL_exp, CD_exp, M_exp, CLp, Mp):
        
        self.CL_exp = CL_exp
        self.CD_exp = CD_exp
        self.M_exp = M_exp
        
        self.CLp = CLp
        self.Mp = Mp
    
        return

    def func_residuo(self, p):
    
        self.R = (p[0] + p[1]*self.CL_exp**2)*(p[2] + p[3]*self.M_exp + p[4]*self.M_exp**2) - self.CD_exp
    
        return self.R
    
    def params(self):
        
        self.param = 0.1*np.ones(5)
        self.param = leastsq(self.func_residuo,self.param)
        self.param = np.array(self.param[0])
        
        return self.param

    def polar(self,x):
    
        self.CL = self.CLp
        self.M  = self.Mp
        self.resp = (x[0] + x[1]*self.CL**2)*(x[2] + x[3]*self.M + x[4]*self.M**2)
        
        self.CD_0 = x[0]*(x[2] + x[3]*self.M + x[4]*self.M**2)
        self.K = x[1]*(x[2] + x[3]*self.M + x[4]*self.M**2)
    
        return self.resp, self.CD_0, self.K
    
    def extra(self,x):  # Extração apenas de CDO e K 
        
        self.CD_0_e = x[0]*(x[2] + x[3]*self.M + x[4]*self.M**2)
        self.K_e = x[1]*(x[2] + x[3]*self.M + x[4]*self.M**2)
        
        return self.CD_0_e, self.K_e




# Dados de entrada:
    
CLp = np.linspace(0.05,1,20)
Mp  = np.linspace(0.05,1,20)

# Polar de arrasto:
problem = drag_polar(CL_exp, CD_exp, M_exp, CLp, Mp)
params = problem.params()
[resp, CD_0, K] = problem.polar(params)

# Extração de CD0 e K, apenas:
[CD_0_e,K_e] = problem.extra(params)



# ----------------------------- GRÁFICOS ------------------------------------#

# Superfície 3D: Mach variável
'''
fig = plt.figure()
ax = plt.axes(projection="3d")
CLp, Mp = np.meshgrid(CLp, Mp)

plot = drag_polar(CL_exp, CD_exp, M_exp, CLp, Mp)
params_plot = plot.params()
resp_graf = plot.polar(params_plot)[0]

ax.plot_surface(Mp,CLp, resp_graf, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none', alpha = 0.9)
ax.scatter3D(M_exp,CL_exp, CD_exp,c = 'red',alpha=1)
#ax.set_title("Polar de arrasto")
ax.set_xlabel("CL")
ax.set_ylabel("M")
ax.set_zlabel("CD")
#plt.savefig("polar_arrasto.svg")
plt.show()
'''

# Curva CL x CD: Mach fixo
'''
fig = plt.figure()
#resp_2d = polar(params,CLp,Mp)[0]
plt.plot(CLp,resp)
#plt.legend(loc = 'best', framealpha = 1)
plt.xlabel("CL")
plt.ylabel("CD")
plt.grid(True)
plt.title("Curva CL x CD")
plt.show()
'''
