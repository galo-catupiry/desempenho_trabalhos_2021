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
            self.v_t0 = 132.5 * kt_to_ms
            
            #Longitudinal
            self.cl = 0.737
            self.cd = 0.95
            self.cl_alpha = 5.0
            self.cd_alpha = .75 
            self.cm_alpha = -.80
            self.cm_alpha_dot = -.30
            self.cm_q = -8.0
            self.cl_delta_e = .4 
            self.cm_delta_e = -.81
            
            #Lateral-Directional
            self.cy_beta = -.72 
            self.cn_beta = .137
            self.cl_beta = -.103
            self.cl_p = -.37
            self.cn_p = -.14
            self.cl_r = .11
            self.cn_r = -.16
            self.cn_delta_a = -.0075
            self.cl_delta_a = .054
            self.cy_delta_r = .175
            self.cn_delta_r = -.063
            self.cl_delta_r = .029
            
        