
"""
INDEPENDENT VARIABLES
---------------------

D_o possible values : [1/2, 3/4, 1, 1+1/4, 1+1/2]*25.4*1e-3 # [m]

Shell_ID possible values : 
[8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27, 29, 31, 33, 35,
37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120]*25.4*1e-3 [m]

Central Spacing Limited Values : 
[Shell_ID/5; 74*D_o**0.75] # To put back in meters and computed with D_o in inches. 

L_shell has free values

"""

"""
COULD BE VARIABLES BUT FIXED
----------------------------

n_passes = 2 # (or 1)

Tube_layout_angle = 45 # [°] (or 60, 90) : 45 and 90 mean square / 60 means triangular

Pitch_ratio constrainted to values depending on D_o, could be varied from 1.2 to 1.5 on its own
Square:
-------
D_o = 1/2     [in] => Pitch_ratio = 1.25
D_o = 3/4     [in] => Pitch_ratio = (1)/D_o
D_o = 1       [in] => Pitch_ratio = (1+1/4)/D_o
D_o = 1 + 1/4 [in] => Pitch_ratio = (1+9/16)/D_o
D_o = 1 + 1/2 [in] => Pitch_ratio = (1+7/8)/D_o
Triangular:
-----------
D_o = 1/2     [in] => Pitch_ratio = 1.25
D_o = 3/4     [in] => Pitch_ratio = (15/16)/D_o
D_o = 1       [in] => Pitch_ratio = (1+1/4)/D_o
D_o = 1 + 1/4 [in] => Pitch_ratio = (1+9/16)/D_o
D_o = 1 + 1/2 [in] => Pitch_ratio = (1+7/8)/D_o

Baffle_cut = 0.25 # Could be varied from 0.15 to 0.4 but 0.25 is usual value for liquid flow
"""

import __init__
from library.component.sizing.heat_exchanger.shell_and_tube.modules.tubes_toolbox import estimate_number_of_tubes, carbon_steel_pipe_thickness, pitch_ratio_fun
from scipy.interpolate import interp1d
from CoolProp.CoolProp import PropsSI
from library.component.sizing.heat_exchanger.basic_sizing_UA import find_UA
from central_spacing_comp import find_divisors_between_bounds
from library.connector.mass_connector import MassConnector
from component.base_component import BaseComponent
 
import pandas as pd
import random
import numpy as np

#%%

class ShellAndTubeSizingOpt(BaseComponent):
    def __init__(self):
        super().__init__()

        self.params = {}

        # Optimization related parameters/variables
        self.opt_vars = {}
        self.bounds = {}
        self.choice_vectors = {}

        # For tube thickness study
        self.P_max_cycle = None
        self.T_max_cycle = None

        # Will be Mass connectors
        self.su_S = None
        self.ex_S = None

        self.su_T = None
        self.ex_T = None

#%%

    def set_opt_vars(self, opt_vars):
        for opt_var in opt_vars:
            self.opt_vars[opt_var] = None 
        return

    def set_opt_vars_values(self, opt_vars_val):
        for key, value in opt_vars_val.items():
            if key in self.opt_vars:
                self.opt_vars[key] = value
            else:
                print(f"Key {key} is not an optimization variable.")
        return

    def set_bounds(self, **bounds):
        for key, value in bounds:
            self.bounds[key] = value 
        return

    def set_choice_vectors(self, choice_vectors):
        for key, value in choice_vectors.items():
            self.choice_vectors[key] = value 
        return

    def set_max_cycle_prop(self, T_max_cycle = None, p_max_cycle = None):
        self.T_max_cycle = T_max_cycle
        self.P_max_cycle = p_max_cycle
        return

    def set_thermo_BC(self, su_S = None, ex_S = None, su_T = None, ex_T = None):
        self.su_S = su_S
        self.ex_S = ex_S

        self.su_T = su_T
        self.ex_T = ex_T
        return 

#%%

    def compute_geom(self):
        """
        Compute rest of geometry
        """

        # Pipe Thickness
        pipe_wall_T = (self.su_S.T + self.su_T.T)/2
        pipe_thickness = carbon_steel_pipe_thickness(self.choice_vectors['D_o_inch'], self.T_max_cycle, self.su_S.p, self.P_max_cycle)

        Tube_t = pipe_thickness[str(self.opt_vars['D_o_inch'])]

        # pitch_ratio
        pitch_ratio = pitch_ratio_fun(self.opt_vars['D_o_inch'], self.params['Layout_angle'])

        # Pipe length
        L_tube = self.opt_vars['L_shell']*self.params['Tube_pass']

        # Cross Passes
        Cross_Passes = self.opt_vars['L_shell']/self.opt_vars['Central_spac'] - 1

        D_o = self.opt_vars['D_o_inch']*25.4*1e-3
        Shell_ID = self.opt_vars['Shell_ID_inch']*25.4*1e-3

        # Number of tubes 
        min_tubes_in_row = 8
        n_tubes = estimate_number_of_tubes(Shell_ID, D_o, pitch_ratio*D_o, self.params['Layout_angle'], min_tubes_in_row)[0]

        # HT Area and HTX volumes
        A_eff = n_tubes*L_tube*np.pi*D_o
        
        T_V_tot = L_tube*n_tubes*np.pi*((D_o - 2*Tube_t)/2)**2

        T_V_out = np.pi*(D_o/2)**2*L_tube*n_tubes
        S_V_tot = self.opt_vars['L_shell']*np.pi*(Shell_ID/2)**2 - T_V_out
            
        self.set_parameters( 
                        A_eff = A_eff, S_V_tot = S_V_tot, Shell_ID = Shell_ID, T_V_tot = T_V_tot, Tube_L = L_tube, 
                        Tube_OD = D_o, Tube_t = Tube_t, central_spacing = self.opt_vars['Central_spac'], 
                        cross_passes = Cross_Passes, n_tubes = n_tubes, pitch_ratio = pitch_ratio
                        ) 

        return

    def Tube_Mass(self):
        
        rho_carbon_steel = 7850 # kg/m^3
        T_mass = np.pi*((self.opt_vars['D_o']/2)**2 - ((self.opt_vars['D_o'] - self.params['Tube_t'])/2)**2) * self.params['L_tube'] * rho_carbon_steel * self.params['n_tubes'] * self.params['n_series'] * self.params['n_passes']

        return T_mass


#%%

"""
Preliminary Sizing Method : Heat Exchangers (Kakac, Liu, Pramuanjaroenkij) - Section 9.3.1
"""

def ShellAndTube_PrelimSizing(Q_dot, T_c_i, T_c_o, T_h_i, T_h_o, params, h_coefs):

    UA_req = find_UA(Q_dot = Q_dot, T_c_i = T_c_i, T_c_o = T_c_o, T_h_i = T_h_i, T_h_o = T_h_o, flow = 'Shell&Tube', params = params)
    U = 0

    for key, value in h_coefs.items():
        U += 1/value

    U = 1/U

    # Required Area
    A_req = UA_req/U
    A_req_OD = A_req*1.2

    # Number of Tubes
    
    # Tube count constant
    if params['Tube_pass'] == 1:
        CTP = 0.93
    elif params['Tube_pass'] == 2:
        CTP = 0.9
    elif params['Tube_pass'] == 3:
        CTP = 0.85
    else:
        print("Tube_pass number not considered.")

    # Tube layout constant
    if params['tube_layout'] == 45 or params['tube_layout'] == 90:
        CL = 1
    elif params['tube_layout'] == 30 or params['tube_layout'] == 60:
        CL = 0.87

    D_s = 0.637*(CL/CTP)**(0.5)*((A_req_OD*((params['pitch_ratio'])**2)*params['Tube_OD'])/params['L_tube'])**(1/2)
    N_tubes = 0.785*(CTP/CL)*(D_s**2)/(params['pitch_ratio']**2*params['Tube_OD']**2)

    return A_req, A_req_OD, D_s, N_tubes

#%%

HX_test = ShellAndTubeSizingOpt()

"""
Optimization related parameters/variables
"""

HX_test.set_opt_vars(['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac'])

choice_vectors = {
                    'D_o_inch' : ['0.5', '0.75', '1', '1.25', '1.5'],
                    'Shell_ID_inch' : ['8', '10', '12', '13.25', '15.25', '17.25', '19.25', '21.25', '23.25', '25', '27',
                                            '29', '31', '33', '35', '37', '39', '42', '45', '48', '54', '60', '66', '72', '78',
                                            '84', '90', '96', '108', '120']
}

HX_test.set_choice_vectors(choice_vectors)

"""
Max T and P for pipe thickness computation
"""

# Worst Case
P_max_cycle = 1048*1e3 # Pa
T_max_cycle = 273.15+147

HX_test.set_max_cycle_prop(T_max_cycle = T_max_cycle, p_max_cycle = P_max_cycle)

"""
Thermodynamical parameters : Inlet and Outlet Design States
"""

su_S = MassConnector()
su_S.set_properties(T = 273.15 + 24, # K
                    P = 101.3*1e3, # Pa
                    m_dot = 1000, # kg/s
                    fluid = 'Water'
                    )

ex_S = MassConnector()
ex_S.set_properties(T = 273.15 + 27.29, # K
                    P = 101.3*1e3, # Pa
                    m_dot = 1000, # kg/s
                    fluid = 'Water'
                    )

su_T = MassConnector()
su_T.set_properties(T = 273.15 + 38.43, # K
                    P = 51.75*1e3, # Pa
                    m_dot = 32.85, # kg/s
                    fluid = 'Cyclopentane'
                    )

ex_T = MassConnector()
ex_T.set_properties(T = 273.15 + 27.21, # K
                    P = 51.75*1e3, # Pa
                    m_dot = 32.85, # kg/s
                    fluid = 'Cyclopentane'
                    )

HX_test.set_thermo_BC(su_S = su_S, ex_S = ex_S, su_T = su_T, ex_T = ex_T)

"""
Parameters Setting
"""

HX_test.set_parameters(
                        Layout_angle = 45, # [°]
                        Tube_pass = 2, # [-]
                        n_series = 1, # [-]
                        Baffle_cut = 25, # [%]
                        foul_t = 0,
                        foul_s = 0,
                        tube_cond = 50, # W/(m*K)
                        
                        Shell_Side = 'C',

                        Flow_Type = 'Shell&Tube', 
                        H_DP_ON = True, 
                        C_DP_ON = True, 
                        n_disc = 50
                      )

"""
Geometry Computation Test
"""

HX_test.set_opt_vars_values(
                            {'D_o_inch' : 0.75,
                            'L_shell' : 5,
                            'Shell_ID_inch' : 90,
                            'Central_spac' : 0.5}
                            )

HX_test.compute_geom()

print(vars(HX_test))

# #%%

# """
# Optimization variables : guesses
# """

# D_o_inch_vect = 

# D_o_inch = 3/4 #3/4 # in
# D_o = D_o_inch*25.4*1e-3 # m

# L_shell = 5 # m

# Shell_ID_inch = 90 # in
# Shell_ID = Shell_ID_inch*25.4*1e-3 # m

# Central_spac_low_lim = Shell_ID/5
# Central_spac_high_lim = (74*D_o_inch**(0.75))*25.4*1e-3 # m

# # Random choice of central spacing givingan integer number of cross passes between the lower and upper bounds
# # Central_spac = random_divisor_between_bounds(L_shell, Central_spac_low_lim, Central_spac_high_lim)

# Central_spac = 0.5 # m 

# #%%

# """
# Preliminary Sizing form Kakac Book
# """

# Q_dot_req = 13300*1e3 # W

# parameters = {
#                 'Tube_pass' : 2,
#                 'tube_layout' : 45,
#                 'n_series' : 1,
#                 'pitch_ratio' : 1.25,
#                 'L_tube' : 10,
#                 'Tube_OD' : 3/4*25.4*1e-3,
#                 'Flow_Type' : 'Shell&Tube'
# }

# h_coefs = {
#             "Shell":6000, # W/(m*K)
#             "Tube":2000, # W/(m*K)
#            }

# A_req, A_req_OD, D_s, N_T = ShellAndTube_PrelimSizing(Q_dot = Q_dot_req, T_c_i = su_S_T, T_c_o = ex_S_T, T_h_i = su_T_T, T_h_o = ex_T_T, params = parameters, h_coefs = h_coefs)

# print(f"A_req : {A_req}")
# print(f"A_req_OD : {A_req_OD}")
# print(f"D_s : {D_s}")
# print(f"N_T : {N_T}")

# #%%

# """
# HTX Model Running
# """

# import __init__

# from library.component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.geometry_shell_and_tube_hx import ShellAndTubeGeom
# from library.component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.simulation_model import HeatExchangerMB

# HX = HeatExchangerMB('Shell&Tube')

# HX.set_inputs(
#     # First fluid
#     Hsu_fluid = su_T_fluid,
#     Hsu_T = su_T_T, # K
#     Hsu_p = su_T_p, # Pa
#     Hsu_m_dot = su_T_m_dot, # kg/s

#     # Second fluid
#     Csu_fluid = su_S_fluid,
#     Csu_T = su_S_T, # K
#     Csu_p = su_S_p, # Pa
#     Csu_m_dot = su_S_m_dot, # kg/s  # Make sure to include fluid information
# )

# "Correlation Loading And Setting"

# Corr_H = {"1P" : "Gnielinski", "2P" : "Horizontal_Tube_Internal_Condensation"}
# Corr_C = {"1P" : "shell_htc_DP_kern", "2P" : "shell_htc_DP_kern"}

# HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31

# "Parameters Setting"

# HX.set_parameters(
#     A_eff = A_eff, Baffle_cut = Baffle_cut, S_V_tot = S_V_tot, Shell_ID = Shell_ID, T_V_tot = T_V_tot, 
#     Tube_L = L_tube, Tube_OD = D_o, Tube_pass = n_passes, Tube_t = Tube_t, central_spacing = Central_spac, 
#     cross_passes = Cross_Passes, foul_s = foul_s, foul_t = foul_t, n_series = n_series, n_tubes = n_tubes, 
#     pitch_ratio = pitch_ratio, tube_cond = tube_cond, tube_layout = Layout_angle_deg,

#     Shell_Side = 'C',

#     Flow_Type = 'Shell&Tube', H_DP_ON = True, C_DP_ON = True, n_disc = 50) 

# HX.set_DP()

# # HX.solve()
