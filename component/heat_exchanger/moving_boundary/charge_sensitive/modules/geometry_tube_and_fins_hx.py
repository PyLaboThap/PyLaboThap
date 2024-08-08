#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:53:29 2023

@author: olivierthome
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.interpolate import interp1d

class TubeAndFinsGeom(object):
    """
    This class aims to be a database for tube and fins heat exchangers geometries
    """
    
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
            
            # Heat transfer area  
            self.A_unfinned = 801.3 # m^2
            self.A_finned = 16501 # m^2
            
            "Performance Data - Tube Side"

            # Fouling factor
            self.fouling = 0.000344 # (m^2*K)/W
            
            "Performance Data - Air Unit Side"
            
            self.A_flow = 256.4/2.64 # From design point : from geometry, can be estimated by 2*7.3*6.88
            self.n_fan = 2
            
            self.altitude = 300 # m
            self.stat_P = 152.5 # Pa
            
            self.Fan_power = 63.1 # kW
            self.Min_design_out_T = 273.15 - 20 # K
            
            "Construction of one bundle"
            self.Finned_tube_flag = 0 # 0 if no fins considered / 1 if fins considered
            
            # HTX dimensions
            self.L = 15 # [m] : length
            self.w = 6.88 # [m] : width
            self.h = 2.55 # [m] : height
                     
            # Tube carac
            self.n_tubes = 249 # m
            self.Tube_OD = 0.0381 # m
            self.Tube_t = 0.00165 # m
            self.Tube_L = self.L # m
            self.Tube_cond = 50 # [W/(m*K)] : tube conduction

            # Tube bundle carac
            self.n_passes = 2 # 2
            self.n_rows = 6
            self.pitch_ratio = 2
            self.pitch = self.pitch_ratio*self.Tube_OD
            self.tube_arrang = "Inline"
            
            # Fins values 
            self.Fin_OD = 0.06985 # m
            self.Fin_per_m = 433
            self.Fin_t = 0.0004 # m
            self.k_fin = 230 # W/(m*K) : Aluminum
            self.Fin_type = "Annular"

            # Area
            N_fins = self.Tube_L*self.Fin_per_m - 1                                # Number of Fins

            A_out_fin = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes*self.n_passes
            A_out_tube = 2*np.pi*(self.Tube_OD/2)*self.n_passes*self.n_tubes*(self.Tube_L - N_fins*self.Fin_t)
            A_out_plate_fin = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes*self.n_passes
            
            self.A_finned = A_out_fin + A_out_tube + A_out_plate_fin
            self.A_unfinned = A_out_tube
            
            self.A_out_tot = self.A_finned
            self.A_in_tot = self.n_tubes*self.n_passes*self.Tube_L*np.pi*(self.Tube_OD-2*self.Tube_t)
            
            # A_out_fin_2 = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes
            # A_out_tube_2 = 2*np.pi*(self.Tube_OD/2)*self.n_tubes*(self.Tube_L - N_fins*self.Fin_t)
            # A_out_plate_fin_2 = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes
            
            # A_out_tot_2 = A_out_fin_2 + A_out_tube_2 + A_out_plate_fin_2
            
            # print(A_out_tot_2)
            
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
            self.T_V_tot = A_in_tube*self.n_tubes*self.Tube_L
            
        elif name == "DECAGONE_ACC_5Bundle":
            
            # Mass for the HTX on a 6M structure (structure included) = 49500 kg
            
            "Performance Data - Tube Side"
            
            # Heat transfer area  
            self.A_unfinned = 667.75 # m^2
            self.A_finned = 13750.83 # m^2

            # Fouling factor
            self.fouling = 0.000344 # (m^2*K)/W
            
            "Performance Data - Air Unit Side"
            
            self.A_flow = 256.4/2.64 # From design point : from geometry, can be estimated by 2*7.3*6.88
            self.n_fan = 2
            
            self.altitude = 300 # m
            self.stat_P = 152.5 # Pa
            
            self.Fan_power = 63.1 # kW
            self.Min_design_out_T = 273.15 - 20 # K
            
            "Construction of one bundle"
            self.Finned_tube_flag = 0 # 0 if no fins considered / 1 if fins considered
            
            # HTX dimensions
            self.L = 15 # [m] : length
            self.w = 6.88 # [m] : width
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
            self.Fin_type = "Annular"

            # Area
            N_fins = self.Tube_L*self.Fin_per_m - 1                                # Number of Fins

            A_out_fin = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes*self.n_passes
            A_out_tube = 2*np.pi*(self.Tube_OD/2)*self.n_passes*self.n_tubes*(self.Tube_L - N_fins*self.Fin_t)
            A_out_plate_fin = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes*self.n_passes
            
            self.A_out_tot = A_out_fin + A_out_tube + A_out_plate_fin
            
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
            self.T_V_tot = A_in_tube*self.n_tubes*self.Tube_L

        elif name == "DECAGONE_ACC_1Bundle":
            
            # Mass for the HTX on a 6M structure (structure included) = 49500 kg
            
            "Performance Data - Tube Side"
            
            # Heat transfer area  
            self.A_unfinned = 133.55 # m^2
            self.A_finned = 2750.17 # m^2

            # Fouling factor
            self.fouling = 0.000344 # (m^2*K)/W
            
            "Performance Data - Air Unit Side"
            
            self.A_flow = 256.4/2.64 # From design point : from geometry, can be estimated by 2*7.3*6.88
            self.n_fan = 2
            
            self.altitude = 300 # m
            self.stat_P = 152.5 # Pa
            
            self.Fan_power = 63.1 # kW
            self.Min_design_out_T = 273.15 - 20 # K
            
            "Construction of one bundle"
            self.Finned_tube_flag = 0 # 0 if no fins considered / 1 if fins considered
            
            # HTX dimensions
            self.L = 15 # [m] : length
            self.w = 6.88 # [m] : width
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
            self.Fin_type = "Annular"
            
            # Area
            N_fins = self.Tube_L*self.Fin_per_m - 1                                # Number of Fins

            A_out_fin = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes*self.n_passes
            A_out_tube = 2*np.pi*(self.Tube_OD/2)*self.n_passes*self.n_tubes*(self.Tube_L - N_fins*self.Fin_t)
            A_out_plate_fin = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes*self.n_passes
            
            self.A_out_tot = A_out_fin + A_out_tube + A_out_plate_fin
            
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
            self.T_V_tot = A_in_tube*self.n_tubes*self.Tube_L
            
        elif name == "DECAGONE_RECUP":
            
            # Mass for the HTX on a 6M structure (structure included) = 49500 kg
            
            "Performance Data - Tube Side"

            # Fouling factor
            self.fouling = 0.9 # (m^2*K)/W
            
            "Performance Data - Air Unit Side"
            
            self.A_flow = 3.36 # !!! From design point : from geometry - considered constant
            
            "Construction of one bundle"
            self.Finned_tube_flag = 1 # 0 if no fins considered / 1 if fins considered
            
            # HTX dimensions
            self.L = 4.5 # [m] : length
            self.w = 1.6*(263.62/536.25) # [m] : bundle width : (263.62/536.25) factor = Bundle width over total HTX diameter  
            self.h = 1.6 # [m] : height
                     
            # Tube carac
            self.n_tubes = 80 # m
            self.Tube_OD = 0.0165 # m
            self.Tube_t = 0.0015 # m
            self.Tube_L = self.L # m
            self.Tube_cond = 50 # [W/(m*K)] : tube conduction

            # Tube bundle carac
            self.n_passes = 6
            self.n_rows = 4
            self.pitch_ratio = 2.125 # !!! assumed from drawings and HT surface area
            self.pitch = self.pitch_ratio*self.Tube_OD
            self.tube_arrang = "Staggered"
            
            # Fins values 
            self.Fin_OD = self.pitch # m
            self.Fin_per_m = 400
            self.Fin_t = 0.00025 # m
            self.k_fin = 50 # W/(m*K) : Steel
            self.Fin_type = "Square"
            
            # Area
            N_fins = self.Tube_L*self.Fin_per_m - 1                                # Number of Fins
            
            if self.Fin_type == "Annular":
                
                A_out_fin = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes*self.n_passes
                A_out_tube = 2*np.pi*(self.Tube_OD/2)*self.n_tubes*self.n_passes*(self.Tube_L - N_fins*self.Fin_t)
                A_out_plate_fin = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes*self.n_passes
                
                self.A_out_tot = A_out_fin + A_out_tube + A_out_plate_fin
                
            elif self.Fin_type == "Square":
                
                "HT Area computations"
                
                # Fin HT area
                N_fins = self.Tube_L*self.Fin_per_m                              # Number of Fins
                Fin_spacing = (self.Tube_L - N_fins*self.Fin_t)/(N_fins)
                
                A_r = 2*(self.Fin_OD**2 - 0.785*self.Tube_OD**2 + 2*self.Fin_OD*self.Fin_t)*(self.Tube_L/Fin_spacing)*self.n_tubes*self.n_passes  # 
                
                L_t = self.Tube_L - N_fins*self.Fin_t
                A_t = np.pi*self.Tube_OD*(self.Tube_L*(1 - self.Fin_t/Fin_spacing)*self.n_tubes*self.n_passes + L_t)
                
                self.A_out_tot = A_r + A_t
            
            else:
                print("Fin geometry is not 'Annular' or 'Square'")
            
            self.A_in_tot = np.pi*(self.Tube_OD - 2*self.Tube_t)*self.Tube_L*self.n_tubes*self.n_passes
                        
            # Heat transfer area  
            self.A_unfinned = self.A_in_tot # m^2
            self.A_finned = self.A_out_tot # m^2
            self.A_finned_DS = 2103 # m^2
            
            # Volume
            A_in_tube = np.pi*((self.Tube_OD - 2*self.Tube_t)/2)**2
            
            self.B_V_tot = self.L*self.w*self.h
            self.T_V_tot = A_in_tube*self.n_tubes*self.n_passes*self.Tube_L
                        
        else: # manual setting
        
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger geometry.")
