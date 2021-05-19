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
from ambiance import Atmosphere

# Planilhas:

xlsx = pd.ExcelFile(parent_dir + '/Planilhas/Dados-JetStar.xlsx')
CL_mach = pd.read_excel(xlsx, 'CL_mach')
CD_mach = pd.read_excel(xlsx, 'CD_mach')

# Dados experimentais:

CL_expe = CL_mach.iloc[:,1]
CD_expe = CD_mach.iloc[:,1]
M_expe  = CL_mach.iloc[:,0]
     
# Funções:

class DragPolar:

    def __init__(self,CL_exp=CL_expe, CD_exp=CD_expe, M_exp=M_expe):
        
        self.CL_exp = CL_exp
        self.CD_exp = CD_exp
        self.M_exp = M_exp
        
        # valores de CL e Mach para serem avaliados
        self.CLp = None # pode ser um vetor
        self.Mp = None # deve ser int ou float
        
        self.params()

    def func_residuo(self, p):
    
        self.R = (p[0] + p[1]*self.CL_exp**2)*(p[2] + p[3]*self.M_exp + p[4]*self.M_exp**2) - self.CD_exp
    
        return self.R
    
    def params(self):
        
        self.param = 0.1*np.ones(5)
        self.param = leastsq(self.func_residuo,self.param)
        self.param = np.array(self.param[0])

    def polar(self,ret=True):
        
        x = self.param
        self.CL = self.CLp
        self.M  = self.Mp
        cd = (x[0] + x[1]*self.CL**2)*(x[2] + x[3]*self.M + x[4]*self.M**2)
        
        self.CD0 = x[0]*(x[2] + x[3]*self.M + x[4]*self.M**2)
        self.K = x[1]*(x[2] + x[3]*self.M + x[4]*self.M**2)
        
        if ret:
            return cd
    

# Dados de entrada:
    
CLp = np.linspace(0.05,1,20)
Mp = np.linspace(0.05,1,20)


# ----------------------------- GRÁFICOS ------------------------------------#

# Superfície 3D: Mach variável

fig = plt.figure()
ax = plt.axes(projection="3d")
CLp_mesh, Mp_mesh = np.meshgrid(CLp, Mp)

drag_plot = DragPolar()
drag_plot.CLp = CLp_mesh
drag_plot.Mp = Mp_mesh
resp_graf = drag_plot.polar()

ax.plot_surface(Mp_mesh,CLp_mesh, resp_graf, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none', alpha = 0.9)
ax.scatter3D(M_expe,CL_expe, CD_expe,c = 'red',alpha=1)
#ax.set_title("Polar de arrasto")
ax.set_xlabel("CL")
ax.set_ylabel("M")
ax.set_zlabel("CD")
#plt.savefig("polar_arrasto.svg")
plt.show()


# Curva CL x CD: Mach fixo
drag_plot.Mp = 0.8
drag_plot.CLp = CLp

fig = plt.figure()
plt.plot(CLp,drag_plot.polar())
plt.xlabel("CL")
plt.ylabel("CD")
plt.grid(True)
plt.title("Curva CL x CD")
plt.show()


## Curvas CDO, K por mach

CD0_list = []
K_list = []

alti_otima = 3800
velo_otima = 128

mach_otimo = velo_otima / Atmosphere(alti_otima).speed_of_sound[0]

mach_list = np.linspace(0.1, 0.8, 10)
drag_plot.CLp = 0.506892114858062
for mach in mach_list:
    
    
    drag_plot.Mp = mach
    drag_plot.polar(False)
    CD0 = drag_plot.CD0
    K = drag_plot.K
    
    CD0_list.append(CD0)
    K_list.append(K)
    
fig_CD0K = plt.figure(figsize = (6,4))
plt.plot(mach_list, CD0_list, label = "CD0")
plt.plot(mach_list, K_list, label = "K")
#plt.vlines(mach_otimo, 0.01, 0.09, color = 'k', ls = '--')
plt.xlabel("Mach", fontsize = 12)
plt.ylabel("$CD_0$, $K$", fontsize = 12)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("cd0_k_mach.pdf")
plt.show()
