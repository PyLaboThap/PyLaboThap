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
            if self.su_wf.T is not None:
                self.inputs['su_wf_T'] = self.su_wf.T
            if self.su_wf.m_dot is not None:
                self.inputs['su_wf_m_dot'] = self.su_wf.m_dot
            if self.su_sf.fluid is not None:
                self.inputs['fluid_sf'] = self.su_sf.fluid
            if self.su_sf.T is not None:
                self.inputs['su_sf_T'] = self.su_sf.T
            if self.su_sf.p is not None:
                self.inputs['su_sf_cp'] = self.su_sf.cp
            if self.su_sf.m_dot is not None:
                self.inputs['su_sf_m_dot'] = self.su_sf.m_dot
            if self.ex_sf.T is not None:
                self.inputs['ex_sf_T'] = self.ex_sf.T
        
        if self.inputs != {}:
            # Working Fluid
            if 'fluid_wf' in self.inputs:
                self.su_wf.set_fluid(self.inputs['fluid_wf'])
            if 'su_wf_T' in self.inputs:
                self.su_wf.set_T(self.inputs['su_wf_T'])
            if 'su_wf_m_dot' in self.inputs:
                self.su_wf.set_m_dot(self.inputs['su_wf_m_dot'])
            if 'su_wf_x' in self.inputs:
                self.su_wf.set_x(self.inputs['su_wf_x'])
            # Secondary Fluid
            if 'fluid_sf' in self.inputs:
                self.su_sf.set_fluid(self.inputs['fluid_sf'])
            if 'su_sf_T' in self.inputs:
                self.su_sf.set_T(self.inputs['su_sf_T'])
            if 'su_sf_p' in self.inputs:
                self.su_sf.set_cp(self.inputs['su_sf_cp'])
            if 'su_sf_m_dot' in self.inputs:
                self.su_sf.set_m_dot(self.inputs['su_sf_m_dot'])
            if 'ex_sf_T' in self.inputs:
                self.ex_sf.set_T(self.inputs['ex_sf_T'])
            
        return['fluid_wf', 'su_wf_T', 'su_wf_m_dot', 'fluid_sf', 'su_sf_T', 'su_sf_cp', 'su_sf_m_dot', 'ex_sf_T']
    
    def get_required_parameters(self):
        return [
            'Pinch', # Minimum pinch point
            'Delta_T_sh', # Superheating or subcooling
        ]
    def get_required_guesses(self):
        return [ 'P_ev' ]
    
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

    def system(self, P_ev):

        T_ev = PropsSI('Tsat', 'P', P_ev, self.su_wf.fluid)

        "Refrigerant side"
        "Liquid zone"
        h_ev_su = PropsSI('H', 'P', P_ev, 'T', self.su_wf.T, self.su_wf.fluid)
        h_ev_l = PropsSI('H', 'P', P_ev, 'Q', 0, self.su_wf.fluid)
        Q_dot_ev_l = self.su_wf.m_dot*(h_ev_l-h_ev_su)

        "Two-phase zone"
        h_ev_v = PropsSI('H', 'P', P_ev, 'Q', 1, self.su_wf.fluid)
        Q_dot_ev_tp = self.su_wf.m_dot*(h_ev_v-h_ev_l)

        "Vapor zone"
        T_wf_ex = self.su_wf.T + self.params['Delta_T_sh']
        h_ev_ex = PropsSI('H', 'P', P_ev, 'T', T_wf_ex, self.su_wf.fluid)
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
 
        # Q_dot_ev_v = m_dot_oil_ev*c_oil*(T_oil_su_ev-T_oil_ev_v)"guess on DeltaT_sh to remove"
        # Q_dot_ev_tp = m_dot_oil_ev*c_oil*(T_oil_ev_v-T_oil_ev_l)
        # Q_dot_ev_l = m_dot_oil_ev*c_oil*(T_oil_ev_l-T_oil_ex_ev) 

        PP = min(self.su_sf.T-T_wf_ex, self.T_sf_l-T_ev)
        res = PP - self.params['Pinch']

        return res


    def solve(self):
        self.check_calculable()
        print(self.guesses, self.calculable, self.parametrized)
        if self.calculable and self.parametrized:
            if self.guesses != {}:
                P_ev_guess = self.guesses['P_ev']
            if self.guesses == {}:
                P_ev_guess = PropsSI('Psat', 'T', self.su_wf.T, self.su_wf.fluid)
            args = ()
            fsolve(self.system, P_ev_guess, args=args)


        # self.ex_wf.set_fluid(self.su_wf.fluid)
        # self.ex_wf.set_T(self.T_wf_out)
        # self.ex_wf.set_p(P_wf_f)
        # self.ex_wf.set_m_dot(self.su_wf.m_dot)
        # self.ex_wf.set_fluid(self.su_wf.fluid)
        # # self.Q_dot.set_Q_dot(self.Q)

        # self.ex_sf.set_p(self.su_sf.p)
        # self.ex_sf.set_m_dot(self.m_dot_sf)
        # self.ex_sf.set_fluid(self.ex_sf.fluid)

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q}")
        # print(f"Q_sub: {self.Q_sub}")
        # print(f"Q_sat: {self.Q_sat}")
        # print(f"Q_sup: {self.Q_sup}")
        # print(f"m_dot_sf: {self.m_dot_sf}")
        # print(f"T_wf_out: {self.T_wf_out}")
        # print(f"tau: {self.tau}")
        # print(f"A_sub: {self.A_sub}")
        # print(f"A_sat: {self.A_sat}")
        # print(f"A_sup: {self.A_sup}")
        print(f"res: {self.res}")
        # print(f"flag_tau: {self.flag_tau}")
        print("======================")





