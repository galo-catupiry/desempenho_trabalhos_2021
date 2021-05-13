"""
Analise da estabilidade da aeronave Lockheed JetStar em condicao de planagem
Codigo referente ao primeiro trabalho da disciplina SAA0183 - Desempenho
de Aeronaves (1o sem 2021)

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""
# file path
import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

# bibliotecas
from ambiance import Atmosphere
from aircraft import JetStar
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from scipy.stats import linregress

# conversao de unidades
lb_to_kg = 0.453592
ft_to_m = 0.3048
psf_to_Pa = 47.880258888889

# ler o arquivo contendo os dados da tabela do CR2144
xlsx = pd.ExcelFile(os.path.join(parent_dir,'Planilhas','Dados-JetStar.xlsx'))

# cria o objeto aircraft
aircraft = JetStar(0,xlsx)
aircraft.read_aerodynamic_data()
aircraft.read_stability_data()
CL_data = aircraft.stability_data

# condicoes da analise
atmo = Atmosphere(aircraft.h_cruise)
Mach = aircraft.V_cruise/atmo.speed_of_sound[0]

#%% tratamento dos dados da tabela
CL_data['H'] = CL_data['H (ft)']*ft_to_m
CL_data['V'] = [CL_data['M'].iloc[i]*Atmosphere(CL_data['H'].iloc[0]).speed_of_sound[0] 
                for i in range(len(CL_data))]
CL_data['CL'] = [CL_data['W (lbs)'].iloc[i]*lb_to_kg*9.81/(CL_data['Q (psf)'].iloc[i]*psf_to_Pa*aircraft.S) 
                 for i in range(len(CL_data))]

# regressao linear de CL_trim vs alpha_trim
mmq = linregress(np.deg2rad(CL_data['Alpha (deg)']),CL_data['CL'])
CL_alpha = mmq.slope
CL_zero = mmq.intercept

# dados de CL e CM da aeronave no Mach desejado
CL_deltaE = np.interp(Mach,CL_data['M'],CL_data['CldeltaE'])
CM_alpha = np.interp(Mach,CL_data['M'],CL_data['Cmalpha'])
CM_deltaE = np.interp(Mach,CL_data['M'],CL_data['CmdeltaE'])

# determinando CM_zero
# sera escolhido deltaE = 0 na condicao de menor alpha trimado
CM_zero = -CM_alpha*np.deg2rad(min(CL_data['Alpha (deg)']))

def CMcg(alpha,deltaE):
    return CM_zero + CM_alpha*alpha + CM_deltaE*deltaE

def CL(alpha,deltaE):
    return CL_zero + CL_alpha*alpha + CL_deltaE*deltaE
  
#%% imprime valores importantes para colocar no relatorio
print(f'CL_alpha = {round(CL_alpha,4)}')
print(f'CL_zero = {round(CL_zero,4)}')
print(f'CL_deltaE = {round(CL_deltaE,4)}')
print(f'CM_alpha = {round(CM_alpha,4)}')
print(f'CM_deltaE = {round(CM_deltaE,4)}')
print(f'CM_zero = {round(CM_zero,4)}')

#%% plots  
deltas = [-10,-5,0,5,10] # graus
alphas = np.linspace(-5,12) # graus
machs = [0.25,0.4,0.5,0.6,0.8]
CMcgs = []
CLs = []
cores = list(pl.cm.jet(np.linspace(0,1,len(deltas))))

for i in range(len(deltas)):
    CMcgs.append(CMcg(np.deg2rad(alphas),np.deg2rad(deltas[i])))
    CLs.append(CL(np.deg2rad(alphas),np.deg2rad(deltas[i])))

# CMcg
fig,ax = plt.subplots(figsize=(6,4))
ax.grid()
for i in range(len(CMcgs)):
    plt.plot(alphas,CMcgs[i],label=f'$\delta_E$ = {deltas[i]}$\degree$',
             color=cores[i])
ax.legend()
plt.xlabel('AoA (deg)')
plt.ylabel('CMcg')
plt.tight_layout()
#plt.savefig('CMcg_alpha.eps')

# CL
fig,ax = plt.subplots(figsize=(6,4))
ax.grid()
for i in range(len(CMcgs)):
    plt.plot(CLs[i],CMcgs[i],label=f'$\delta_E$ = {deltas[i]}$\degree$',
             color=cores[i])
ax.legend()
plt.xlabel('CL')
plt.ylabel('CMcg')
plt.tight_layout()
#plt.savefig('CMcg_CL.eps')    

# CL vs alpha
fig,ax = plt.subplots()
ax.grid()
ax.plot(CL_data['Alpha (deg)'],np.deg2rad(CL_data['Alpha (deg)'])*CL_alpha+CL_zero,
        color='k',ls='dashdot',label='MMQ')
ax.scatter(CL_data['Alpha (deg)'],CL_data['CL'],color='r',label='Experimental')
ax.set_xlabel('$AoA_{trim}$ [deg]')
ax.set_ylabel('CL')
ax.legend()


# CL x alpha x Mach
alpha_Lzero = -CL_zero/CL_alpha
fig,ax = plt.subplots(figsize=(8,4))
ax.grid()
for i in range(len(machs)):
    cl_alphai = np.interp(machs[i],aircraft.CL_alpha_mach['mach'],
                          aircraft.CL_alpha_mach['CL_alpha'])
    plt.plot(alphas,cl_alphai*np.deg2rad(alphas-alpha_Lzero),
             color=cores[i],label=f'Mach = {machs[i]}')
ax.set_xlabel('AoA [deg]')
ax.set_ylabel('CL')
ax.legend()
plt.tight_layout()
#plt.savefig('CL_mach.eps')