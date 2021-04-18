# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 13:12:16 2021

"""

#Aircraft: Lockheed JetStar


#Conversion of units
lb_to_kg = 0.453592
slug_ft2_to_kg_m2 = 1.3558179619 
ft_to_m = 0.3048
ft2_to_m2 = 0.092903
deg_to_rad = 0.0174533
kt_to_ms = 0.51444

class jetstar():
    
    def __init__(self, condition):
        
        self.condition = condition
        
        #Fixed features
        self.W = 382204 * lb_to_kg
        self.Ix = 118773  * slug_ft2_to_kg_m2
        self.Iy = 135869 * slug_ft2_to_kg_m2
        self.Iz = 243504 * slug_ft2_to_kg_m2
        self.Ixz = 5061 * slug_ft2_to_kg_m2
        self.S = 542.5 * ft2_to_m2
        self.b = 53.75 * ft_to_m
        self.c_barra = 10.93 * ft_to_m
        
        #Power Approach Non-Dimensional Stability Derivatives
        if self.condition == 1:
            self.alpha_zero = 6.5 * deg_to_rad
            self.V_t0 = 132.5 * kt_to_ms
            
            #Longitudinal
            self.CL = 0.737
            self.CD = 0.95
            self.CL_alpha = 5.0
            self.CD_alpha = .75 
            self.Cm_alpha = -.80
            self.Cm_alpha_dot = -.30
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
            
        