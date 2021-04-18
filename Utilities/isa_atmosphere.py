"""
International Standard Atmosphere
Grupo E
"""

import math

#constantes
rho0=1.225;# [kg/m^3]
p0=101325;# [Pa]
t0=288.15;# [K]
beta=-0.0065;# [K/m]
a0=334.1; # [m/s]
g=9.81; #[m/s^2]
R=287.053; # [J/kg.K]
gama = 1.4 # adimensional
referencia = 11000 # [m] separação entre as duas camadas

#Função que forneça [P,a,rho,T] dada uma altitude em metros
def atmo_isa(alti):
     
    if alti <= referencia: #trecho linear
        t=t0+beta*alti 
        p = (t/t0)**(-g/(beta*R)); p=p*p0 
        a_=math.sqrt(gama*R*t) 
        rho_ = (t/t0)**((-g/(beta*R))-1); rho_=rho_*rho0
                          
    else: #trecho constante
        t=216.65
        a_=math.sqrt(gama*R*t)
        p= math.exp((-g*(alti-referencia))/(R*t));p=p*22620.47
        rho_=math.exp((-g*(alti-referencia))/(R*t));rho_=rho_*0.3639
            
    return round(t,3),round(p,3),round(a_,3),round(rho_,3)

