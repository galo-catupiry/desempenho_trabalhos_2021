"""
Arquivo principal do trabalho T1 de Desempenho.

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""

from Desempenho import alcance_autonomia_CL, alcance_autonomia_V

altitude = 13105 # [m]
velocidade = 811 / 3.6 # [m/s]

print("----- Alcance e autonomia [Caso CL constante] -----")
print("Delta_X = {} [m]".format(round(alcance_autonomia_CL(altitude, velocidade)[0],2)))
print("t = {} [s]".format(round(alcance_autonomia_CL(altitude, velocidade)[1],2)))
alcance_autonomia_CL(altitude, velocidade, True, True)


print("----- Alcance e autonomia [Caso V constante] -----")
print("Delta_X = {} [m]".format(round(alcance_autonomia_V(altitude, velocidade)[0],2)))
print("t = {} [s]".format(round(alcance_autonomia_V(altitude, velocidade)[1],2)))
alcance_autonomia_V(altitude, velocidade, True, True)
