from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

# from component.heat_exchanger.moving_boundary.simple_model.modules.U import U_Gnielinski_calibrated, U_DittusBoelter, U_Cooper_calibrater, U_Thonon

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import numpy as np
import math

class HXPinchCst(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su_wf = MassConnector() # Working fluid supply
        self.su_sf = MassConnector() # Secondary fluid supply
        self.ex_wf = MassConnector()
        self.ex_sf = MassConnector()
        self.Q_dot = HeatConnector()
        self.guesses = {}

    def get_required_inputs(self):

        if self.inputs == {}:
            # Hot Fluid
            if self.su_wf.fluid is not None:
                self.inputs['fluid_wf'] = self.su_wf.fluid
            if self.su_wf.h is not None:
                self.inputs['su_wf_h'] = self.su_wf.h
            if self.su_wf.m_dot is not None:
                self.inputs['su_wf_m_dot'] = self.su_wf.m_dot
            if self.su_sf.fluid is not None:
                self.inputs['fluid_sf'] = self.su_sf.fluid
            if self.su_sf.T is not None:
                self.inputs['su_sf_T'] = self.su_sf.T
            if self.su_sf.cp is not None:
                self.inputs['su_sf_cp'] = self.su_sf.cp
            if self.su_sf.m_dot is not None:
                self.inputs['su_sf_m_dot'] = self.su_sf.m_dot
            if self.ex_sf.T is not None:
                self.inputs['ex_sf_T'] = self.ex_sf.T
        
        if self.inputs != {}:
            # Working Fluid
            if 'fluid_wf' in self.inputs:
                self.su_wf.set_fluid(self.inputs['fluid_wf'])
            if 'su_wf_h' in self.inputs:
                self.su_wf.set_h(self.inputs['su_wf_h'])
            if 'su_wf_m_dot' in self.inputs:
                self.su_wf.set_m_dot(self.inputs['su_wf_m_dot'])
            # if 'su_wf_x' in self.inputs:
            #     self.su_wf.set_x(self.inputs['su_wf_x'])
            # Secondary Fluid
            if 'fluid_sf' in self.inputs:
                self.su_sf.set_fluid(self.inputs['fluid_sf'])
            if 'su_sf_T' in self.inputs:
                self.su_sf.set_T(self.inputs['su_sf_T'])
            if 'su_sf_cp' in self.inputs:
                self.su_sf.set_cp(self.inputs['su_sf_cp'])
            if 'su_sf_m_dot' in self.inputs:
                self.su_sf.set_m_dot(self.inputs['su_sf_m_dot'])

        return['fluid_wf', 'su_wf_h', 'su_wf_m_dot', 'fluid_sf', 'su_sf_T', 'su_sf_cp', 'su_sf_m_dot']
    
    def get_required_parameters(self):
        return [
            'Pinch', # pinch point
            'Delta_T_sh_sc', # Superheating or subcooling
            'type_HX' # Evaporator or condenser
        ]
    def get_required_guesses(self):
        return [ 'P_sat' ]
    
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - su_wf: fluid={self.su_wf.fluid}, T={self.su_wf.T}, p={self.su_wf.p}, m_dot={self.su_wf.m_dot}")
        print(f"  - su_sf: fluid={self.su_sf.fluid}, T={self.su_sf.T}, p={self.su_sf.p}, m_dot={self.su_sf.m_dot}")
        print(f"  - ex_wf: fluid={self.ex_wf.fluid}, T={self.ex_wf.T}, p={self.ex_wf.p}, m_dot={self.ex_wf.m_dot}")
        print(f"  - ex_sf: fluid={self.ex_sf.fluid}, T={self.ex_sf.T}, p={self.ex_sf.p}, m_dot={self.ex_sf.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")

        print("\nInputs:")
        for input in self.get_required_inputs():
            if input in self.inputs:
                print(f"  - {input}: {self.inputs[input]}")
            else:
                print(f"  - {input}: Not set")

        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")

        print("======================")

    def system_evap(self, x):
        P_ev = x[0]
        T_ev = PropsSI('T', 'P', P_ev, 'Q', 0.5, self.su_wf.fluid)

        "Refrigerant side"
        "Liquid zone"
        h_ev_su =  self.su_wf.h #PropsSI('H', 'P', P_ev, 'T', self.su_wf.T, self.su_wf.fluid)
        h_ev_l = PropsSI('H', 'P', P_ev, 'Q', 0, self.su_wf.fluid)
        Q_dot_ev_l = self.su_wf.m_dot*(h_ev_l-h_ev_su)
        "Two-phase zone"
        h_ev_v = PropsSI('H', 'P', P_ev, 'Q', 1, self.su_wf.fluid)
        Q_dot_ev_tp = self.su_wf.m_dot*(h_ev_v-h_ev_l)

        "Vapor zone"
        self.T_wf_ex = self.su_wf.T + self.params['Delta_T_sh_sc']
        h_ev_ex = PropsSI('H', 'P', P_ev, 'T', self.T_wf_ex, self.su_wf.fluid)
        Q_dot_ev_v = self.su_wf.m_dot*(h_ev_ex-h_ev_v)

        "Total heat transfer"
        Q_dot_ev = Q_dot_ev_l + Q_dot_ev_tp + Q_dot_ev_v

        "Secondary fluid side"
        "Total heat transfer"
        self.T_sf_ex = self.su_sf.T + Q_dot_ev/(self.su_sf.m_dot*self.su_sf.cp)
        "Vapor zone"
        self.T_sf_v = self.su_sf.T - Q_dot_ev_v/(self.su_sf.m_dot*self.su_sf.cp)
        "Two-phase zone"
        self.T_sf_l = self.T_sf_v - Q_dot_ev_tp/(self.su_sf.m_dot*self.su_sf.cp)
 
        PP = min(self.su_sf.T-self.T_wf_ex, self.T_sf_l-T_ev)
        res = PP - self.params['Pinch']
        self.Q = Q_dot_ev
        self.su_wf.set_p(P_ev)
        return res

    def system_cond(self, x):
        P_cd = x[0]
        print('P_cd', P_cd)
        T_cd = PropsSI('T', 'P', P_cd, 'Q', 0.5, self.su_wf.fluid)
        print('T_cd', T_cd)
        "3. Condenser model"
        "3.1 Refrigerant side"
        h_wf_su_cd = self.su_wf.h # PropsSI('H', 'P', P_cd, 'T', self.su_wf.T, self.su_wf.fluid)
        print('h_wf_su_cd', h_wf_su_cd)
        "3.1.1 Vapor zone"
        h_wf_cd_v = PropsSI('H', 'P', P_cd, 'Q', 1, self.su_wf.fluid)
        Q_dot_cd_v = self.su_wf.m_dot*(h_wf_su_cd-h_wf_cd_v)
 
        "3.1.2 Two_phase zone"
        h_wf_cd_l = PropsSI('H', 'P', P_cd, 'Q', 0, self.su_wf.fluid)
        Q_dot_cd_tp = self.su_wf.m_dot*(h_wf_cd_v-h_wf_cd_l)
 
        "3.1.3 Liquid zone"
        self.T_wf_ex = T_cd-self.params['Delta_T_sh_sc']
        h_wf_ex_cd = PropsSI('H', 'P', P_cd, 'T', self.T_wf_ex, self.su_wf.fluid)
        Q_dot_cd_l = self.su_wf.m_dot*(h_wf_cd_l-h_wf_ex_cd)
 
        "3.1.4 Total power"
        Q_dot_cd = Q_dot_cd_v+Q_dot_cd_tp+Q_dot_cd_l 
 
        "3.2 Water side"
        self.T_sf_ex = self.su_sf.T + Q_dot_cd/(self.su_sf.m_dot*self.su_sf.cp)

        T_sf_cd_l = self.su_sf.T + Q_dot_cd_l/(self.su_sf.m_dot*self.su_sf.cp)
        T_sf_cd_v = T_sf_cd_l + Q_dot_cd_tp/(self.su_sf.m_dot*self.su_sf.cp)
 
        "3.3 Pinch point position"
        Pinch_cd=min(self.T_wf_ex-self.su_sf.T, T_cd-T_sf_cd_v)

        res = Pinch_cd-self.params['Pinch']   
        self.Q = Q_dot_cd
        self.su_wf.set_p(P_cd)
        print('P_cd', P_cd)
        return res


    def solve(self):
        self.check_calculable()
        print('Je rentre dans l échangeur')
        if self.calculable and self.parametrized:
            if self.params['type_HX'] == 'evaporator':
                print('je rentre dans l évapo')
                if self.guesses != {}:
                    P_ev_guess = self.guesses['P_sat']
                if self.guesses == {}:
                    P_ev_guess = PropsSI('P', 'T', self.su_wf.T+20, 'Q', 0.5, self.su_wf.fluid)
                x = [P_ev_guess]
                try:
                    fsolve(self.system_evap, x)
                except:
                    print('Convergence error in evaporator model')
                    return

            elif self.params['type_HX'] == 'condenser':
                print('coucou')
                if self.guesses != {}:
                    P_cd_guess = self.guesses['P_sat']
                if self.guesses == {}:
                    P_cd_guess = PropsSI('P', 'T', self.su_wf.T-20, 'Q', 0.5, self.su_wf.fluid)
                x = [P_cd_guess]
                try:
                    fsolve(self.system_cond, x)
                except:
                    print('Convergence error in condenser model')
                    return

            self.ex_wf.set_fluid(self.su_wf.fluid)
            self.ex_wf.set_T(self.T_wf_ex)
            self.ex_wf.set_p(self.su_wf.p)
            self.ex_wf.set_m_dot(self.su_wf.m_dot)

            self.ex_sf.set_fluid(self.su_sf.fluid)
            self.ex_sf.set_m_dot(self.su_sf.m_dot)
            self.ex_sf.set_T(self.T_sf_ex)
            self.defined = True

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q}")
        print("======================")





