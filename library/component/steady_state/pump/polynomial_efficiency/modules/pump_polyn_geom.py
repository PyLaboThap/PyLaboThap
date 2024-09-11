# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:02:31 2024

@author: Basile
"""

import numpy as np
import matplotlib.pyplot as plt

def EtaPumpFun(V_dot):
    "V_dot : Pump volumetric flowrate in m^3/s"
    "eta_is_pump = 70%"
    "eta_pump : Total efficiency of the pump (from electrical power to power transmitted to the wf"
    
    coefs_pump = [2.85653,2.45967,-0.02657,0.00009]
    
    eta_pump = 0
    
    for i in range(len(coefs_pump)):
        eta_pump = eta_pump + coefs_pump[i]*V_dot**i    
    
    return eta_pump

class GeometryPolynPump(object):
    
    def __init__(self, **kwargs):

        "Pump data"
        
        self.min_flowrate = None
        self.rated_flowrate = None  
        self.max_flowrate = None               
        
        self.N_pp_rated = None
        self.pump_voltage = None
        self.pump_phases = None
        
        "Swept volume"
        self.V_swept = None
        
        "Pump total efficiency"    
        self.coefs_pump = None
        self.eta_polyn = None

    def set_parameters(self, name, **kwargs):
        
        if name == "DECAGONE_PUMP":
            
            "Pump data"
            
            self.min_flowrate = 21 # m^3/h    
            self.rated_flowrate = 67.8 # m^3/h    
            self.max_flowrate = 74.6 # m^3/h                
            
            self.N_pp_rated = 50 # Hz : !!! Assumption
            self.pump_voltage = 400 # V 
            self.pump_phases = 3
            
            "Volumetric efficiency"
            self.eta_v = 0.95 # !!: arbitrary
            
            "Swept volume"
            self.V_swept = self.rated_flowrate/(3600*self.N_pp_rated*self.eta_v) # m^3 / rev
            
            "Electrical Motor efficiency"
            self.eta_max_motor = 0.95 # !!! Assumption
            self.W_dot_el_rated = 84800
            
            "Pump total efficiency"    
            self.coefs_pump = [0.00009,-0.02657,2.45967,2.85653]
            self.eta_tot = np.poly1d(self.coefs_pump)
            
            "Mechanical efficiency"
            self.eta_m = 0.95 # !!! arbitrary
            
            # "Isentropic Efficiency"
            # self.eta_is = 0.7
            
        else: # manual setting
        
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger geometry.")
