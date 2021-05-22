"""
Arquivo principal do trabalho T1 de Desempenho.

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
    
    
As funções abaixo permitem analisar as condições gerais (CL constante ou V constante),
o gráfico de razão de descida por velocidade, e as condições envolvendo os parâmetros ótimos

Descomente a função que deseja usar, alterando os valores de interesse na variável 'cond'.
A chave 'condicao' deve ser "max_range" ou "max_endurance".

Em caso de dúvidas, contate um membro do grupo.    
    
"""

import os, sys
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from datetime import datetime
from Desempenho import (alcance_autonomia_CL, alcance_autonomia_V, 
                        hdot_V, parametros_otimos_V, parametros_otimos_CL)

#%% ------ MAIN -------
start = datetime.now()

altitude = 13105 # [m]
velocidade = 811 / 3.6  # [m/s]

cond = {'condicao': 'max_endurance', 'h': altitude, 'v': velocidade}

# =========== Caso Geral =========== #
# alcance_autonomia_CL(cond.get('v'), False, False, 
#                                         cond.get('h'))

# alcance_autonomia_V(cond.get('v'), False, False,
#                     cond.get('h'))


# hdot_V(cond.get('v'), False, 
#         cond.get('h'))

# =========== Parâmetros Ótimos =========== #

# resp_V = parametros_otimos_V(cond)
# resp_CL = parametros_otimos_CL(cond)

# if (cond.get('condicao') == 'max_range'):
#     print(" ----- Máximo alcance para V constante -----")
#     print("Delta_X = {} [m]".format(round(sum(resp_V),2)))
#     print(" ----- Máximo alcance para CL constante -----")
#     print("Delta_X = {} [m]".format(round(sum(resp_CL),2)))
    
# elif (cond.get('condicao') == 'max_endurance'):
#     print(" ----- Máxima autonomia para V constante -----")
#     print("Delta_t = {} [s]".format(round(sum(resp_V),2)))
#     print(" ----- Máxima autonomia para CL constante -----")
#     print("Delta_t = {} [s]".format(round(sum(resp_CL),2)))

print(datetime.now() - start)

