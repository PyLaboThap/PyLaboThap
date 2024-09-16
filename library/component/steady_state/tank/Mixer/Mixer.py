# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile
"""

from Port.Mass_connector import Mass_connector
from CoolProp.CoolProp import PropsSI
import numpy as np

class Mixer(object):
    
    class geom():
            pass 

    def __init__(self, geom):
        
        
        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.defined = False
        
        "Geometry sub-object"
        self.geom = geom # parameters 
        
        "Input"
        self.point_su = np.full(geom.n_inlet, None)

        "Outputs"
        self.point_ex = [None]
                
#%%    
    def check_calculable(self):
        if self.point_su[0] != None and self.point_su[1] != None:
            if (self.point_su[0].completely_known) and (self.point_su[1].completely_known):
                self.calculable = True
        
    def check_parametrized(self):
        if self.geom != None:
            self.parametrized = True
           
#%%
    def set_parameters(self, **kwargs):
            """
            Set parameters of the heat exchanger.
    
            Parameters
            ----------
            **kwargs : dict
                Key-value pairs representing parameters and their values.
                
                Example of call : heat_exchanger.set_parameters(D_o=0.05, t=0.002, H_core=2.0)
            """
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger.")
            
            self.check_parametrized()
            
#%%

    def solve(self):
        
        "1) Check if pressures are the same"
        
        for i in range(len(self.point_su)):
            if self.point_su[0].p != self.point_su[i].p:
                print("Not the same pressure at all mixer inputs")
                return
        
        "2) Compute output"
        
        self.point_ex[0] = Mass_connector()
        self.point_ex[0].set_fluid(self.point_su[0].fluid)
        self.point_ex[0].set_p(self.point_su[0].p)
        
        m_dot_out = 0
        
        for i in range(len(self.point_su)):
            m_dot_out = m_dot_out + self.point_su[i].m_dot
            
        self.point_ex[0].set_m_dot(m_dot_out)
        
        h_out_mean = 0
        
        for i in range(len(self.point_su)):
            h_out_mean = h_out_mean + self.point_su[i].m_dot*self.point_su[i].h
        
        h_out_mean = h_out_mean/m_dot_out
        
        self.point_ex[0].set_h(h_out_mean)
        
        self.defined = True
        
        return
        
        
            