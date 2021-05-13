"""
Arquivo principal do trabalho T1 de Desempenho.

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

from Desempenho import alcance_autonomia_CL, alcance_autonomia_V

altitude = 4000 # [m]
velocidade = 400 / 3.6 # [m/s]

print("----- Alcance e Autonomia [Caso CL constante] -----")
deltaX_CL, t_CL = alcance_autonomia_CL(altitude, velocidade, True, False, False, False)
print("Delta_X = {} [m]".format(round(deltaX_CL,2)))
print("t = {} [s]".format(round(t_CL,2)))



print("----- Alcance e autonomia [Caso V constante] -----")
deltaX_V, t_V = alcance_autonomia_V(altitude, velocidade, True, False)
print("Delta_X = {} [m]".format(round(deltaX_V,2)))
print("t = {} [s]".format(round(t_V,2)))

