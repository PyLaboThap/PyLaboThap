#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:53:29 2023

@author: olivierthome
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.interpolate import interp1d

class ShellAndTubeGeom(object):
    """
    This class aims to be a database for shell and tube heat exchangers geometries
    """
    
    def __init__(self, **kwargs):

        "Process Conditions"
        
        # Fouling factors
        self.foul_s = None # m^2*K/W
        self.foul_t = None # m^2*K/W
        
        "HTX Performance"      
        
        # Nominal Heat transfer coefficients
        
        self.Shell_h = None # W/(m^2*K)
        self.Tube_h = None # W/(m^2*K)
        
        # !!! What is EMTD ? 
        self.EMTD = None # °
        
        # U coefficients (actual and required => HTX oversized)
        self.U_act = None # W/(m^2*K)
        self.U_req = None # W/(m^2*K)
        
        # Total duty and effective area
        self.Q_tot = None # W
        self.A_eff = None # m^2
        
        # Overdesign
        self.OD = None
        
        "Shell Geometry"
        
        # TEMA : Tubular Exchangers Manufacturers Association
        self.TEMA_Type = None
        
        # Shell Internal D
        self.Shell_ID = None # m
        
        # !!! ?
        self.n_series = None
        self.parallel = None
        self.orient = None # !!! : along which axis ? °
        
        "Baffle Geometry"
        
        # Baffle type
        self.Baffle_type = None
        
        # Baffle Cut
        self.Baffle_cut = None # % of diameter
        
        # Baffle Orientation
        self.Baffle_orient = None
        
        # Central spacing
        self.central_spacing = None # m 
        
        # Crosspasses
        self.cross_passes = None
        
        "Tube Geometry"
        
        self.Tube_type = None
        self.Tube_OD = None # m
        self.Tube_t = None # m 
        self.Tube_L = None # m

        # Tube leakage
        self.D_OTL = None # m : Outer Tube Limit diameter diameter in which all the tubes are)
            
        self.pitch_ratio = None
        self.layout = None # °
        self.n_tubes = None
        self.Tube_pass = None
        self.Tubesheet_t = None # Total tubesheet thickness
        self.tube_cond = None # W/(m*K)
        
        "Nozzle Values"
        
        self.shell_D_in = None # m
        self.shell_D_out = None # m
        self.inlet_height = None # m
        self.outlet_height = None # m
        self.Tube_inlet = None # m
        
        "Clearances"
        self.clear_TB = None # m : Tube to baffle clearance
        self.clear_BDLS = None # m : Bundle to shell clearance
        self.clear_BS = None # m : Baffle to shell clearance    
        self.N_strips = None # m : Number of sealing strips
        
        "Volume"
        self.S_V_tot = None # m^3
        self.T_V_tot = None # m^3
        
        "Type heat Exchanger"
        self.Type_HX = "Shell_and_Tube" # "Crossflow" # "Shell_and_Tube"
        
    def set_parameters(self, name, **kwargs):
        
        if name == "DECAGONE_EVAP_Equ":

            "Process Conditions"
            
            # Fouling factors
            self.foul_s = 0.000045 # m^2*K/W
            self.foul_t = 0 # m^2*K/W
            
            "HTX Performance"      
            
            # Nominal Heat transfer coefficients
            
            self.Shell_h = 1268.7 # W/(m^2*K)
            self.Tube_h = 1300.2 # W/(m^2*K)
            
            # !!! What is EMTD ? 
            self.EMTD = 31.5 # °
            
            # U coefficients (actual and required => HTX oversized)
            self.U_act = 542.1 # W/(m^2*K)
            self.U_req = 457.06 # W/(m^2*K)
            
            # Total duty and effective area
            self.Q_tot = 7935*1e3 # W
            self.A_eff = 551 # 275.5 # m^2
            
            # Overdesign
            self.OD = self.U_act/self.U_req
            
            "Shell Geometry"
            
            # TEMA : Tubular Exchangers Manufacturers Association
            self.TEMA_Type = "BFU"
            
            # Shell Internal D
            self.Shell_ID = 742.95*1e-3 # m
            
            # !!! ?
            self.n_series = 2 # 1 # 
            self.parallel = 1
            self.orient = 0 # !!! : along which axis ? °
            
            "Baffle Geometry"
            
            # Baffle type
            self.Baffle_type = "Double_Seg"
            
            # Baffle Cut
            self.Baffle_cut = 24.62 # % of diameter
            
            # Baffle Orientation
            self.Baffle_orient = "Parallel"
            
            # Central spacing
            self.central_spacing = 295.33*1e-3 # m 
            
            # Inlet/Outlet spacing
            self.inlet_spacing = 400.05*1e-3 # m 
            self.outlet_spacing = 400.05*1e-3 # m 
         
            # Crosspasses
            self.cross_passes = 24
            
            "Tube Geometry"
            
            self.Tube_type = "Plain"
            self.Tube_OD = 15.875*1e-3 # m
            self.Tube_t = 1.651*1e-3 # m 
            self.Tube_L = 7.315 # m
            self.D_OTL = 726.26*1e-3 # m : Outer Tube Limit diameter diameter in which all the tubes are)
            
            self.pitch_ratio = 1.4
            self.tube_layout = 45 # °
            self.n_tubes = 389
            self.Tube_pass = self.n_series*2
            self.Tubesheet_t = 96.837*1e-3 # [m] Total tubesheet thickness
            self.tube_cond = 50 # W/(m*K) steel Tube conductivity

            "Nozzle Values"
            
            self.shell_D_in = 146.33*1e-3 # m
            self.shell_D_out = 146.33*1e-3 # m
            self.inlet_height = 31.338*1e-3 # m
            self.outlet_height = 31.338*1e-3 # m
            self.Tube_inlet = 122.25*1e-3 # m

            "Clearances"
            self.clear_TB = 0.7937*1e-3 # m : Tube to baffle clearance
            self.clear_BDLS = 16.687*1e-3 # m : Bundle to shell clearance
            self.clear_BS = 4.7625*1e-3 # m : Baffle to shell clearance  
            self.N_strips = 0 # m : Number of sealing strips
            
            "Volume"
            T_V_tot_out = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*(self.Tube_OD/2)**2 # m^3
            
            self.T_V_tot = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*((self.Tube_OD-2*self.Tube_t)/2)**2 # m^3
            self.S_V_tot = (np.pi*(self.Shell_ID/2)**2 * self.Tube_L) - T_V_tot_out

        elif name == "DECAGONE_EVAP":

            "Process Conditions"
            
            # Fouling factors
            self.foul_s = 0.000045 # m^2*K/W
            self.foul_t = 0 # m^2*K/W
            
            "HTX Performance"      
            
            # Nominal Heat transfer coefficients
            
            self.Shell_h = 1268.7 # W/(m^2*K)
            self.Tube_h = 1300.2 # W/(m^2*K)
            
            # !!! What is EMTD ? 
            self.EMTD = 31.5 # °
            
            # U coefficients (actual and required => HTX oversized)
            self.U_act = 542.1 # W/(m^2*K)
            self.U_req = 457.06 # W/(m^2*K)
            
            # Total duty and effective area
            self.Q_tot = 7935*1e3 # W
            self.A_eff = 275.5 # m^2
            
            # Overdesign
            self.OD = self.U_act/self.U_req
            
            "Shell Geometry"
            
            # TEMA : Tubular Exchangers Manufacturers Association
            self.TEMA_Type = "BFU"
            
            # Shell Internal D
            self.Shell_ID = 742.95*1e-3 # m
            
            # !!! ?
            self.n_series = 1 # 1 # 
            self.parallel = 1
            self.orient = 0 # !!! : along which axis ? °
            
            "Baffle Geometry"
            
            # Baffle type
            self.Baffle_type = "Double_Seg"
            
            # Baffle Cut
            self.Baffle_cut = 24.62 # % of diameter
            
            # Baffle Orientation
            self.Baffle_orient = "Parallel"
            
            # Central spacing
            self.central_spacing = 295.33*1e-3 # m 
            
            # Inlet/Outlet spacing
            self.inlet_spacing = 400.05*1e-3 # m 
            self.outlet_spacing = 400.05*1e-3 # m 
         
            # Crosspasses
            self.cross_passes = 24
            
            "Tube Geometry"
            
            self.Tube_type = "Plain"
            self.Tube_OD = 15.875*1e-3 # m
            self.Tube_t = 1.651*1e-3 # m 
            self.Tube_L = 7.315 # m
            self.D_OTL = 726.26*1e-3 # m : Outer Tube Limit diameter diameter in which all the tubes are)
            
            self.pitch_ratio = 1.4
            self.tube_layout = 45 # °
            self.n_tubes = 389
            self.Tube_pass = self.n_series*2
            self.Tubesheet_t = 96.837*1e-3 # [m] Total tubesheet thickness
            self.tube_cond = 50 # W/(m*K) steel Tube conductivity

            "Nozzle Values"
            
            self.shell_D_in = 146.33*1e-3 # m
            self.shell_D_out = 146.33*1e-3 # m
            self.inlet_height = 31.338*1e-3 # m
            self.outlet_height = 31.338*1e-3 # m
            self.Tube_inlet = 122.25*1e-3 # m

            "Clearances"
            self.clear_TB = 0.7937*1e-3 # m : Tube to baffle clearance
            self.clear_BDLS = 16.687*1e-3 # m : Bundle to shell clearance
            self.clear_BS = 4.7625*1e-3 # m : Baffle to shell clearance  
            self.N_strips = 0 # m : Number of sealing strips
            
            "Volume"
            T_V_tot_out = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*(self.Tube_OD/2)**2 # m^3
            
            self.T_V_tot = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*((self.Tube_OD-2*self.Tube_t)/2)**2 # m^3
            self.S_V_tot = (np.pi*(self.Shell_ID/2)**2 * self.Tube_L) - T_V_tot_out
            
        elif name == "Valid_Sam":
            
            "HTX Performance"      
            
            self.A_eff = 110 # m^2
            
            "Shell Geometry"
            
            # Shell Internal D
            self.Shell_ID = 539.8*1e-3 # m
            
            # !!! ?
            self.n_series = 1 # number of passes
            self.parallel = 1
            self.orient = 0 # !!! : along which axis ? °
            
            "Baffle Geometry"
            
            # Baffle Cut
            self.Baffle_cut = 25 # % of diameter
            
            # Baffle Orientation
            self.Baffle_orient = "Parallel"
            
            # Central spacing
            self.central_spacing = 127*1e-3 # m 
            # Inlet/Outlet spacing
            self.inlet_spacing = 0.738*self.central_spacing # 400.05*1e-3 # m 
            self.outlet_spacing = 0.738*self.central_spacing # 400.05*1e-3 # m 

            # Crosspasses
            self.cross_passes = 38
            
            "Tube Geometry"
            
            self.Tube_type = "Plain"
            self.Tube_OD = 25.4*1e-3 # m 
            self.Tube_t = 2.413*1e-3/2 # m  # 1.651*1e-3 # m 
            self.Tube_L = 4.8768 # m
            self.D_OTL = 0.9775*(self.Shell_ID) # m : Outer Tube Limit diameter diameter in which all the tubes are)
            
            self.pitch_ratio = 1.25
            self.tube_layout = 0 # °
            self.n_tubes = 158
            self.Tube_pass = 2
            self.Tubesheet_t = 0.013*self.Tube_L # [m] Total tubesheet thickness
            self.tube_cond = 50 # W/(m*K) steel Tube conductivity

            "Nozzle Values"
            
            self.shell_D_in = None # 146.33*1e-3 # m
            self.shell_D_out = None # 146.33*1e-3 # m
            self.inlet_height = None # 31.338*1e-3 # m
            self.outlet_height = None # 31.338*1e-3 # m
            self.Tube_inlet = None # 122.25*1e-3 # m

            "Clearances"
            self.clear_TB = 0.8*1e-3 # m : Tube to baffle clearance
            self.clear_BDLS = 35*1e-3 # m : Bundle to shell clearance
            self.clear_BS = 5*1e-3 # m : Baffle to shell clearance  
            self.N_strips = 0 # Number of sealing strips 
            
            "Volume"
            T_V_tot_out = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*(self.Tube_OD/2)**2 # m^3
            
            self.T_V_tot = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*((self.Tube_OD-2*self.Tube_t)/2)**2 # m^3
            self.S_V_tot = (np.pi*(self.Shell_ID/2)**2 * self.Tube_L) - T_V_tot_out

        elif name == "kakac_Ex_9_3":
            
            "HTX Performance"      
            
            self.A_eff = 112 # m^2
            
            "Shell Geometry"
            
            # Shell Internal D
            self.Shell_ID = 580*1e-3 # m
            
            # !!! ?
            self.n_series = 2 # number of passes
            self.parallel = 1
            self.orient = 0 # !!! : along which axis ? °
            
            "Baffle Geometry"
            
            # Baffle Cut
            self.Baffle_cut = 25 # % of diameter
            
            # Baffle Orientation
            self.Baffle_orient = "Parallel"
            
            # Central spacing
            self.central_spacing = 500*1e-3 # m 
            # Inlet/Outlet spacing
            self.inlet_spacing = self.central_spacing/0.738 # 400.05*1e-3 # m 
            self.outlet_spacing = self.central_spacing/0.738 # 400.05*1e-3 # m 

            # Crosspasses
            self.cross_passes = 38
            
            "Tube Geometry"
            
            self.Tube_type = "Plain"
            self.Tube_OD = 19*1e-3 # m 
            self.Tube_t = 1.5*1e-3 # m  # 1.651*1e-3 # m 
            self.Tube_L = 4.8768 # m
            self.D_OTL = 0.9775*(self.Shell_ID) # m : Outer Tube Limit diameter diameter in which all the tubes are)
            
            self.pitch_ratio = 1.33
            self.tube_layout = 0 # °
            self.n_tubes = 374
            self.Tube_pass = 1
            self.Tubesheet_t = 0.013*self.Tube_L # [m] Total tubesheet thickness
            self.tube_cond = 42.3 # W/(m*K) steel Tube conductivity

            "Nozzle Values"
            
            self.shell_D_in = None # 146.33*1e-3 # m
            self.shell_D_out = None # 146.33*1e-3 # m
            self.inlet_height = None # 31.338*1e-3 # m
            self.outlet_height = None # 31.338*1e-3 # m
            self.Tube_inlet = None # 122.25*1e-3 # m

            "Clearances"
            self.clear_TB = 0.8*1e-3 # m : Tube to baffle clearance
            self.clear_BDLS = 35*1e-3 # m : Bundle to shell clearance
            self.clear_BS = 5*1e-3 # m : Baffle to shell clearance  
            self.N_strips = 36 # m : Number of sealing strips        

            "Volume"
            T_V_tot_out = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*(self.Tube_OD/2)**2 # m^3
            
            self.T_V_tot = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*((self.Tube_OD-2*self.Tube_t)/2)**2 # m^3
            self.S_V_tot = (np.pi*(self.Shell_ID/2)**2 * self.Tube_L) - T_V_tot_out

        elif name == "kakac_Ex_12_2":
            
            "HTX Performance"      
            
            self.A_eff = 37.12 # m^2
            
            "Shell Geometry"
            
            # Shell Internal D
            self.Shell_ID = 0.387 # m
            
            # !!! ?
            self.n_series = 1 # number of passes
            self.parallel = 1
            self.orient = 0 # !!! : along which axis ? °
            
            "Baffle Geometry"
            
            # Baffle Cut
            self.Baffle_cut = 25 # % of diameter
            
            # Baffle Orientation
            self.Baffle_orient = "Parallel"
            
            # Central spacing
            self.central_spacing = 0.35 # m 
            # Inlet/Outlet spacing
            self.inlet_spacing = self.central_spacing/0.738 # 400.05*1e-3 # m 
            self.outlet_spacing = self.central_spacing/0.738 # 400.05*1e-3 # m 

            # Crosspasses
            self.cross_passes = None
            
            "Tube Geometry"
            
            self.Tube_type = "Plain"
            self.Tube_OD = 0.01905 # m 
            self.Tube_t = 0.889*1e-3 # m  # 1.651*1e-3 # m 
            self.Tube_L = 4.53 # m
            self.D_OTL = 0.9775*(self.Shell_ID) # m : Outer Tube Limit diameter diameter in which all the tubes are)
            
            self.pitch_ratio = 25.4/self.Tube_OD
            self.tube_layout = 45 # °
            self.n_tubes = 137
            self.Tube_pass = 1
            self.Tubesheet_t = 0.013*self.Tube_L # [m] Total tubesheet thickness
            self.tube_cond = 121 # W/(m*K) Brass (Laiton) Tube conductivity

            "Nozzle Values"
            
            self.shell_D_in = None # 146.33*1e-3 # m
            self.shell_D_out = None # 146.33*1e-3 # m
            self.inlet_height = None # 31.338*1e-3 # m
            self.outlet_height = None # 31.338*1e-3 # m
            self.Tube_inlet = None # 122.25*1e-3 # m

            "Clearances"
            self.clear_TB = 0.6*1e-3 # m : Tube to baffle clearance
            self.clear_BDLS = 30*1e-3 # m : Bundle to shell clearance
            self.clear_BS = 3*1e-3 # m : Baffle to shell clearance  
            self.N_strips = 0 # Number of sealing strips      

            "Volume"
            T_V_tot_out = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*(self.Tube_OD/2)**2 # m^3
            
            self.T_V_tot = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*((self.Tube_OD-2*self.Tube_t)/2)**2 # m^3
            self.S_V_tot = (np.pi*(self.Shell_ID/2)**2 * self.Tube_L) - T_V_tot_out
            
        elif name == 'DECAGONE_RECUP':
            
            # Mass for the HTX on a 6M structure (structure included) = 49500 kg

            "HTX Performance"      
            
            self.A_eff = 37.12 # m^2
            
            "Shell Geometry"
            
            # Shell Internal D
            self.Shell_ID = 0.387 # m
            
            # !!! ?
            self.n_series = 1 # number of passes
            self.parallel = 1
            self.orient = 0 # !!! : along which axis ? °
            
            "Baffle Geometry"
            
            # Baffle Cut
            self.Baffle_cut = 25 # % of diameter
            
            # Baffle Orientation
            self.Baffle_orient = "Parallel"
            
            # Central spacing
            self.central_spacing = 0.35 # m 
            # Inlet/Outlet spacing
            self.inlet_spacing = self.central_spacing/0.738 # 400.05*1e-3 # m 
            self.outlet_spacing = self.central_spacing/0.738 # 400.05*1e-3 # m 

            # Crosspasses
            self.cross_passes = None
            
            "Tube Geometry"
            
            self.Tube_type = "Plain"
            self.Tube_OD = 0.01905 # m 
            self.Tube_t = 0.889*1e-3 # m  # 1.651*1e-3 # m 
            self.Tube_L = 4.53 # m
            self.D_OTL = 0.9775*(self.Shell_ID) # m : Outer Tube Limit diameter diameter in which all the tubes are)
            
            self.pitch_ratio = 25.4/self.Tube_OD
            self.tube_layout = 45 # °
            self.n_tubes = 137
            self.Tube_pass = 1
            self.Tubesheet_t = 0.013*self.Tube_L # [m] Total tubesheet thickness
            self.tube_cond = 121 # W/(m*K) Brass (Laiton) Tube conductivity

            "Nozzle Values"
            
            self.shell_D_in = None # 146.33*1e-3 # m
            self.shell_D_out = None # 146.33*1e-3 # m
            self.inlet_height = None # 31.338*1e-3 # m
            self.outlet_height = None # 31.338*1e-3 # m
            self.Tube_inlet = None # 122.25*1e-3 # m

            "Clearances"
            self.clear_TB = 0.6*1e-3 # m : Tube to baffle clearance
            self.clear_BDLS = 30*1e-3 # m : Bundle to shell clearance
            self.clear_BS = 3*1e-3 # m : Baffle to shell clearance  
            self.N_strips = 0 # Number of sealing strips      

            "Volume"
            T_V_tot_out = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*(self.Tube_OD/2)**2 # m^3
            
            self.T_V_tot = self.Tube_pass*self.n_tubes*self.Tube_L*np.pi*((self.Tube_OD-2*self.Tube_t)/2)**2 # m^3
            self.S_V_tot = (np.pi*(self.Shell_ID/2)**2 * self.Tube_L) - T_V_tot_out

            # ------------------------------------------------------------------------------
            
            "Performance Data - Tube Side"

            # Fouling factor
            self.fouling = 0.9 # (m^2*K)/W
            
            "Performance Data - Air Unit Side"
            
            self.A_flow = 3.36 # !!! From design point : from geometry - considered constant
            
            "Construction of one bundle"
            self.Finned_tube_flag = 0 # 0 if no fins considered / 1 if fins considered
            
            # HTX dimensions
            self.L = 4.5 # [m] : length
            self.w = 1.6 # [m] : width
            self.h = 1.6 # [m] : height
                     
            # Tube carac
            self.n_tubes = 480 # m
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
                
                A_out_fin = 2*np.pi*(self.Fin_OD/2)*self.Fin_t*N_fins*self.n_tubes
                A_out_tube = 2*np.pi*(self.Tube_OD/2)*self.n_tubes*(self.Tube_L - N_fins*self.Fin_t)
                A_out_plate_fin = 2*np.pi*((self.Fin_OD/2)**2 - (self.Tube_OD/2)**2)*N_fins*self.n_tubes
                
                self.A_out_tot = A_out_fin + A_out_tube + A_out_plate_fin
                
            elif self.Fin_type == "Square":
                
                "HT Area computations"
                
                # Fin HT area
                N_fins = self.Tube_L*self.Fin_per_m                              # Number of Fins
                Fin_spacing = (self.Tube_L - N_fins*self.Fin_t)/(N_fins)
                
                A_r = 2*(self.Fin_OD**2 - 0.785*self.Tube_OD**2 + 2*self.Fin_OD*self.Fin_t)*(self.Tube_L/Fin_spacing)*self.n_tubes  # 
                
                L_t = self.Tube_L - N_fins*self.Fin_t
                A_t = np.pi*self.Tube_OD*(self.Tube_L*(1 - self.Fin_t/Fin_spacing)*self.n_tubes + L_t)
                
                self.A_out_tot = A_r + A_t
            
            else:
                print("Fin geometry is not 'Annular' or 'Square'")
            
            self.A_in_tot = np.pi*(self.Tube_OD - 2*self.Tube_t)*self.Tube_L*self.n_tubes
            
            print(self.A_in_tot)
            
            # Heat transfer area  
            self.A_unfinned = self.A_in_tot # m^2
            self.A_finned = self.A_out_tot # m^2
            self.A_finned_DS = 2103 # m^2

            print(self.A_finned)

            
            # Volume
            A_in_tube = np.pi*((self.Tube_OD - 2*self.Tube_t)/2)**2
            
            self.B_V_tot = self.L*self.w*self.h
            self.T_V_tot = A_in_tube*self.n_tubes*self.n_passes*self.Tube_L
                        
            
            
        else: # manual setting
            print("No heat exchanger geometry has this name")
            
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger geometry.")

