#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:53:29 2023

@author: olivierthome
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.interpolate import interp1d

class GeometryCrossFlowFins(object):
    
    def __init__(self, **kwargs):

        "Performance Data - Tube Side"
        
        # Heat transfer area  
        self.A_unfinned = None # m^2
        self.A_finned = None # m^2

        # Fooling factor
        self.fooling = None # (m^2*K)/W
        
        "Performance Data - Air Unit Side"
        
        self.A_flow = None # From design point : from geometry, can be estimated by 2*7.3*6.88
        self.n_fan = None
        
        self.altitude = None # m
        self.stat_P = None # Pa
        
        self.Fan_power = None # kW
        self.Min_design_out_T = None # K
        
        "Construction of one bundle"
        self.Finned_tube_flag = None # 0 if no fins considered / 1 if fins considered
        
        # HTX dimensions
        self.L = None # [m] : length
        self.w = None # [m] : width
        self.h = None # [m] : height
                 
        # Tube carac
        self.n_tubes = None # m
        self.Tube_OD = None # m
        self.Tube_t = None # m
        self.Tube_L = None # m
        self.Tube_cond = None # [W/(m*K)] : tube conduction

        # Tube bundle carac
        self.n_passes = None
        self.n_rows = None
        self.pitch_ratio = None
        self.pitch = None
        
        # Fins values 
        self.Fin_OD = None # m
        self.Fin_per_m = None
        self.Fin_t = None # m
                
        # Motor Value
        self.mot_power_per_drive = None # [W]
        self.mot_V = None # [V]
        self.alim_N = None # [V]
        self.mot_N = None # [RPM]

        # Fan Value
        self.Fan_D = None # m
        self.Fan_n_blades = None # aluminum blades
        self.Fan_N = None # RPM : HTD-belt frequecy drive of 7.71 ratio
        self.Power_per_fan = None # kW

    def set_parameters(self, name, **kwargs):
        
        if name == "DECAGONE_ACC":
            
            # Mass for the HTX on a 6M structure (structure included) = 49500 kg
            
            "Performance Data - Tube Side"
            
            # Heat transfer area  
            self.A_unfinned = 801.3 # m^2
            self.A_finned = 16501 # m^2

            # Fouling factor
            self.fouling = 0.000344 # (m^2*K)/W
            
            "Performance Data - Air Unit Side"
            
            self.A_flow = 256.4/(2*2.64) # From design point : from geometry, can be estimated by 2*7.3*6.88
            self.n_fan = 2
            
            self.altitude = 300 # m
            self.stat_P = 152.5 # Pa
            
            self.Fan_power = 63.1 # kW
            self.Min_design_out_T = 273.15 - 20 # K
            
            "Construction of one bundle"
            self.Finned_tube_flag = 1 # 0 if no fins considered / 1 if fins considered
            self.Fin_type = "Annular"

            # HTX dimensions
            self.L = 15 # [m] : length
            self.w = 6.88/2 # [m] : width
            self.h = 2.55 # [m] : height
                     
            # Tube carac
            self.n_tubes = 249 # m
            self.Tube_OD = 0.0381 # m
            self.Tube_t = 0.00165 # m
            self.Tube_L = self.L # m
            self.Tube_cond = 50 # [W/(m*K)] : tube conduction

            # Tube bundle carac
            self.n_passes = 2
            self.n_rows = 6
            self.pitch_ratio = 2
            self.pitch = self.pitch_ratio*self.Tube_OD
            self.tube_arrang = "Inline"
            
            # Fins values 
            self.Fin_OD = 0.06985 # m
            self.Fin_per_m = 433
            self.Fin_t = 0.0004 # m
            self.k_fin = 200 # W/(m*K) : Aluminum
            
            # Area
            N_fins = self.Tube_L*self.Fin_per_m - 1                                # Number of Fins

            A_out_fin = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes
            A_out_tube = 2*np.pi*(self.Tube_OD/2)*self.n_tubes*(self.Tube_L - N_fins*self.Fin_t)            #
            A_out_plate_fin = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes
            
            self.A_out_tot = A_out_fin + A_out_tube + A_out_plate_fin
            self.A_in_tot = np.pi*(self.Tube_OD - 2*self.Tube_t)*self.Tube_L*self.n_tubes
            
            # print(self.A_out_tot)
            # print(self.A_in_tot)
            
            "Fan characteristics"

            # Motor Value
            self.mot_power_per_drive = 45*1e3 # [W]
            self.mot_V = 400 # [V]
            self.alim_N = 50 # [V]
            self.mot_N = 1480 # [RPM]

            # Fan Value
            self.Fan_D = 5.182 # m
            self.Fan_n_blades = 4 # aluminum blades
            self.Fan_N = 192 # RPM : HTD-belt frequecy drive of 7.71 ratio
            self.Power_per_fan = 31.5 # kW
            
            A_in_tube = np.pi*((self.Tube_OD - 2*self.Tube_t)/2)**2
            
            self.B_V_tot = self.L*self.w*self.h
            self.T_V_tot = A_in_tube*self.n_tubes*self.n_passes*self.Tube_L
            
        if name == "DECAGONE_ACC_5Bundle":
            
            # Mass for the HTX on a 6M structure (structure included) = 49500 kg
            
            "Performance Data - Tube Side"
            
            # Heat transfer area  
            self.A_unfinned = 667.75 # m^2
            self.A_finned = 13750.83 # m^2

            # Fouling factor
            self.fouling = 0.000344 # (m^2*K)/W
            
            "Performance Data - Air Unit Side"
            
            self.A_flow = 256.4/(2*2.64) # From design point : from geometry, can be estimated by 2*7.3*6.88
            self.n_fan = 2
            
            self.altitude = 300 # m
            self.stat_P = 152.5 # Pa
            
            self.Fan_power = 63.1 # kW
            self.Min_design_out_T = 273.15 - 20 # K
            
            "Construction of one bundle"
            self.Finned_tube_flag = 1 # 0 if no fins considered / 1 if fins considered
            self.Fin_type = "Annular"

            # HTX dimensions
            self.L = 15 # [m] : length
            self.w = 6.88/2 # [m] : width
            self.h = 2.55 # [m] : height
                     
            # Tube carac
            self.n_tubes = 207 # m
            self.Tube_OD = 0.0381 # m
            self.Tube_t = 0.00165 # m
            self.Tube_L = self.L # m
            self.Tube_cond = 50 # [W/(m*K)] : tube conduction

            # Tube bundle carac
            self.n_passes = 1
            self.n_rows = 5
            self.pitch_ratio = 2
            self.pitch = self.pitch_ratio*self.Tube_OD
            self.tube_arrang = "Inline"
            
            # Fins values 
            self.Fin_OD = 0.06985 # m
            self.Fin_per_m = 433
            self.Fin_t = 0.0004 # m
            self.k_fin = 200 # W/(m*K) : Aluminum
            
            # Area
            N_fins = self.Tube_L*self.Fin_per_m - 1                                # Number of Fins

            A_out_fin = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes
            A_out_tube = 2*np.pi*(self.Tube_OD/2)*self.n_tubes*(self.Tube_L - N_fins*self.Fin_t)            #
            A_out_plate_fin = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes
            
            self.A_out_tot = A_out_fin + A_out_tube + A_out_plate_fin
            self.A_in_tot = np.pi*(self.Tube_OD - 2*self.Tube_t)*self.Tube_L*self.n_tubes
            
            # print(self.A_out_tot)
            # print(self.A_in_tot)
                     
            "Fan characteristics"

            # Motor Value
            self.mot_power_per_drive = 45*1e3 # [W]
            self.mot_V = 400 # [V]
            self.alim_N = 50 # [V]
            self.mot_N = 1480 # [RPM]

            # Fan Value
            self.Fan_D = 5.182 # m
            self.Fan_n_blades = 4 # aluminum blades
            self.Fan_N = 192 # RPM : HTD-belt frequecy drive of 7.71 ratio
            self.Power_per_fan = 31.5 # kW
            
            A_in_tube = np.pi*((self.Tube_OD - 2*self.Tube_t)/2)**2
            
            self.B_V_tot = self.L*self.w*self.h
            self.T_V_tot = A_in_tube*self.n_tubes*self.n_passes*self.Tube_L

        if name == "DECAGONE_ACC_1Bundle":
            
            # Mass for the HTX on a 6M structure (structure included) = 49500 kg
            
            "Performance Data - Tube Side"
            
            # Heat transfer area  
            self.A_unfinned = 133.55 # m^2
            self.A_finned = 2750.17 # m^2

            # Fouling factor
            self.fouling = 0.000344 # (m^2*K)/W
            
            "Performance Data - Air Unit Side"
            
            self.A_flow = 256.4/(2*2.64) # From design point : from geometry, can be estimated by 2*7.3*6.88
            self.n_fan = 2
            
            self.altitude = 300 # m
            self.stat_P = 152.5 # Pa
            
            self.Fan_power = 63.1 # kW
            self.Min_design_out_T = 273.15 - 20 # K
            
            "Construction of one bundle"
            self.Finned_tube_flag = 1 # 0 if no fins considered / 1 if fins considered
            self.Fin_type = "Annular"

            # HTX dimensions
            self.L = 15 # [m] : length
            self.w = 6.88/2 # [m] : width
            self.h = 2.55 # [m] : height
                     
            # Tube carac
            self.n_tubes = 42 # m
            self.Tube_OD = 0.0381 # m
            self.Tube_t = 0.00165 # m
            self.Tube_L = self.L # m
            self.Tube_cond = 50 # [W/(m*K)] : tube conduction

            # Tube bundle carac
            self.n_passes = 1
            self.n_rows = 1
            self.pitch_ratio = 2
            self.pitch = self.pitch_ratio*self.Tube_OD
            self.tube_arrang = "Inline"
            
            # Fins values 
            self.Fin_OD = 0.06985 # m
            self.Fin_per_m = 433
            self.Fin_t = 0.0004 # m
            self.k_fin = 200 # W/(m*K) : Aluminum
            
            # Area
            N_fins = self.Tube_L*self.Fin_per_m - 1                                # Number of Fins

            A_out_fin = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes
            A_out_tube = 2*np.pi*(self.Tube_OD/2)*self.n_tubes*(self.Tube_L - N_fins*self.Fin_t)            #
            A_out_plate_fin = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes
            
            self.A_out_tot = A_out_fin + A_out_tube + A_out_plate_fin
            self.A_in_tot = np.pi*(self.Tube_OD - 2*self.Tube_t)*self.Tube_L*self.n_tubes
            
            # print(self.A_out_tot)
            # print(self.A_in_tot)
            
            "Fan characteristics"

            # Motor Value
            self.mot_power_per_drive = 45*1e3 # [W]
            self.mot_V = 400 # [V]
            self.alim_N = 50 # [V]
            self.mot_N = 1480 # [RPM]

            # Fan Value
            self.Fan_D = 5.182 # m
            self.Fan_n_blades = 4 # aluminum blades
            self.Fan_N = 192 # RPM : HTD-belt frequecy drive of 7.71 ratio
            self.Power_per_fan = 31.5 # kW
            
            A_in_tube = np.pi*((self.Tube_OD - 2*self.Tube_t)/2)**2
            
            self.B_V_tot = self.L*self.w*self.h
            self.T_V_tot = A_in_tube*self.n_tubes*self.n_passes*self.Tube_L
            
            
        else: # manual setting
        
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger geometry.")
