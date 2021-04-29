"""
Parametros geometricos e aerodinamicos da aeronave Lockheed JetStar
Codigo referente ao primeiro trabalho da disciplina SAA0183 - Desempenho
de Aeronaves (1o sem 2021)

Integrantes:
    Abner Micael de Paula Souza - 10788676
    Alessandro Melo de Oliveira - 10788662
    Guilherme Beppu de Souza    - 10696681
    Thiago Buchignani De Amicis - 10277418
"""
from numpy import cos, sqrt, deg2rad
import pandas as pd
# Aircraft: Lockheed JetStar

# Conversion of units
lb_to_kg = 0.453592
slug_ft2_to_kg_m2 = 1.3558179619 
ft_to_m = 0.3048
ft2_to_m2 = 0.092903
deg_to_rad = 0.0174533
kt_to_ms = 0.51444

class JetStar():

    def __init__(self, condition, xlsx):
        self.xlsx = xlsx
        self.condition = condition
        
        # Fixed features
        
        self.P = 38204 * lb_to_kg
        self.W = self.P * 9.81
        self.Ix = 118773  * slug_ft2_to_kg_m2
        self.Iy = 135869 * slug_ft2_to_kg_m2
        self.Iz = 243504 * slug_ft2_to_kg_m2
        self.Ixz = 5061 * slug_ft2_to_kg_m2
        # wing
        self.S = 542.5 * ft2_to_m2 #wing area
        self.b = 53.75 * ft_to_m #wing span
        self.c_barra = 10.93 * ft_to_m 
        self.c_tip = 1.62
        self.c_root = 4.44
        self.AR = self.b**2 / self.S #aspct ratio
        self.WL = self.W / self.S #wing loading
        # self.e = 0.7 #oswald number
        # self.K = (3.1415 * self.e * self.AR)**-1 #induced  drag factor
        self.lambd = self.c_tip/self.c_root
        self.sweep = deg2rad(28.7) #sweep back
        # tailpane
        self.S_t = 27.75
        self.i_t = 0 # incidence
        
        # positions relative to DATUM
        # wing
        self.x_w = 9.96 # aerodinamic center at 0.25c_barra
        self.z_w = 0
        self.x_t = 15.63 # aerodinamic center at 0.25c_t
        self.z_t = 2.31
        
        # cg
        self.x_cg = self.x_w
        self.z_cg = self.z_w
        
        # velocity
        self.V_t0 = 132.5 * kt_to_ms
        
        # downwash
        Kh = (1 - abs(self.z_t/self.b))/pow(2*(self.x_t-self.x_w)/self.b,1/3)
        Klambda = (10 - 3*self.lambd)/7
        Kar = 1/self.AR - 1/(1+pow(self.AR,1.7))
        self.dEdalpha = 4.44*pow(Kh*Klambda*Kar*sqrt(cos(self.sweep/4)),1.19)
        
        # Power Approach Non-Dimensional Stability Derivatives
        if self.condition == 1:
            self.alpha_zero = 6.5 * deg_to_rad
            self.V_t0 = 132.5 * kt_to_ms
            
            #Longitudinal
            self.CL = 0.737
            self.CD = 0.095
            self.CL_alpha = 5.0
            self.CD_alpha = .75 
            self.Cm_alpha = -.80
            self.Cm_alpha_dot = -3.0
            self.Cm_q = -8.0
            self.CL_delta_e = .4 
            self.Cm_delta_e = -.81
            
            #Lateral-Directional
            self.Cy_beta = -.72 
            self.Cn_beta = .137
            self.Cl_beta = -.103
            self.Cl_p = -.37
            self.Cn_p = -.14
            self.Cl_r = .11
            self.Cn_r = -.16
            self.Cn_delta_a = -.0075
            self.Cl_delta_a = .054
            self.Cy_delta_r = .175
            self.Cn_delta_r = -.063
            self.Cl_delta_r = .029
            
            
    def read_aerodynamic_data(self):
        self.alpha_zero_mach = pd.read_excel(self.xlsx,'alpha_zero_mach')
        self.CL_mach = pd.read_excel(self.xlsx, 'CL_mach')
        self.CD_mach = pd.read_excel(self.xlsx, 'CD_mach')
        self.CL_alpha_mach = pd.read_excel(self.xlsx, 'CL_alpha_mach')
        self.CD_alpha_mach = pd.read_excel(self.xlsx, 'CD_alpha_mach')
        self.Cm_alpha_mach = pd.read_excel(self.xlsx, 'Cm_alpha_mach')
        self.Cm_alphadot_mach = pd.read_excel(self.xlsx,'Cm_alphadot_mach')
        self.Cm_q_mach = pd.read_excel(self.xlsx,'Cm_q_mach')
        self.CD_M_mach = pd.read_excel(self.xlsx,'CD_M_mach')
        self.CM_M_mach = pd.read_excel(self.xlsx,'CM_M_mach')
        self.CL_delta_e = pd.read_excel(self.xlsx,'CL_delta_e')
        self.Cm_delta_e = pd.read_excel(self.xlsx,'Cm_delta_e')