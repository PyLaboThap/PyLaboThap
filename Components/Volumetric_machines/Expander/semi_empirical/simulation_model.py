# -*- coding: utf-8 -*-

from components.base_component import BaseComponent
from connectors.mass_connector import MassConnector
from connectors.work_connector import WorkConnector
from connectors.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import numpy as np
import time

class ExpanderSE(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector() # Mass_connector
        self.work_exp = WorkConnector()
        self.heat_amb = HeatConnector()

    def get_required_inputs(self):
        
        if self.inputs == {}:
            if self.su.T is not None:
                self.inputs['su_T'] = self.su.T
            elif self.su.h is not None:
                self.inputs['su_h'] = self.su.h
            if self.su.p is not None:
                self.inputs['su_p'] = self.su.p
            if self.ex.p is not None:
                self.inputs['ex_p'] = self.ex.p
            if self.work_exp.speed is not None:
                self.inputs['N_rot'] = self.work_exp.speed
            if self.heat_amb.temperature_in is not None:
                self.inputs['T_amb'] = self.heat_amb.temperature_in
            if self.su.fluid is not None:
                self.inputs['su_fluid'] = self.su.fluid
        
        if self.inputs != {}:
            self.su.set_fluid(self.inputs['su_fluid'])
            if 'su_T' in self.inputs:
                self.su.set_T(self.inputs['su_T'])
            elif 'su_h' in self.inputs:
                self.su.set_h(self.inputs['su_h'])
            if 'su_p' in self.inputs:
                self.su.set_p(self.inputs['su_p'])
            if 'ex_p' in self.inputs:
                self.ex.set_p(self.inputs['ex_p'])
            if 'N_rot' in self.inputs:
                self.work_exp.set_speed(self.inputs['N_rot'])
            if 'T_amb' in self.inputs:
                self.heat_amb.set_temperature_in(self.inputs['T_amb'])

        return ['su_p', 'su_T', 'ex_p', 'N_rot', 'T_amb', 'su_fluid']

    def get_required_parameters(self):
        return [
            'AU_amb', 'AU_su_n', 'AU_ex_n', 'd_su1', 'm_dot_n', 
            'A_leak', 'W_dot_loss_0', 'alpha', 'C_loss', 'rv_in', 'V_s'
        ]
    
    def print_setup(self):
        print("=== Expander Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")
        print(f"  - W_dot: speed={self.work_exp.speed}")
        print(f"  - Q_dot_amb: temperature_in={self.heat_amb.temperature_in}")

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

    def System(self, x):
        self.m_dot, self.T_w = x
        self.N = self.inputs['N_rot'] /60
            
        #Boundary on the mass flow rate
        self.m_dot = max(self.m_dot, 1e-5)
        
        "1. Supply conditions: su"
        T_su = self.su.T
        P_su = self.su.p
        h_su = self.su.h
        s_su = self.su.s
        rho_su = self.su.D
        Fluid = self.su.fluid
        P_ex = self.ex.p
        self.P_ex = P_ex

        h_ex_is = PropsSI('H', 'P', P_ex, 'S', s_su, Fluid)
        h_max = PropsSI('H', 'P', 4e6, 'T', 500, Fluid)

        T_sat_su = PropsSI('T', 'P', P_su, 'Q', 1, Fluid)
        if T_su<T_sat_su:
            print('----Warning the compressor inlet stream is not in vapor phase---')

        "2. Supply pressure drop: su->su1"
        h_su1 = h_su #Isenthalpic valve
        # Assumption that the density doesn't change too much
        A_su = np.pi*(self.params['d_su1']/2)**2
        V_dot_su = self.m_dot/rho_su
        C_su = V_dot_su/A_su
        h_su_thr1 = h_su-(C_su**2)/2
        h_su_thr = max(h_su_thr1, h_ex_is)
        P_su_thr = PropsSI('P', 'S', s_su, 'H', h_su_thr, Fluid)
        P_su1 = max(P_su_thr, P_ex+1)
        self.DP_su = P_su-P_su1
        T_su1 = PropsSI('T', 'P', P_su1, 'H', h_su1, Fluid)

        "3. Cooling at the entrance: su1->su2"  
        try:
            cp_su1 = PropsSI('CPMASS', 'P', P_su1, 'H', h_su1, Fluid)
        except:
            cp_su1 = PropsSI('CPMASS', 'P', P_su1, 'Q', 0, Fluid)
        AU_su = self.params['AU_su_n']*(self.m_dot/self.params['m_dot_n'])**(0.8)
        C_dot_su = self.m_dot*cp_su1
        NTU_su = AU_su/C_dot_su
        epsilon_su1 = max(0,1-np.exp(-NTU_su))
        Q_dot_su = max(0, epsilon_su1*self.m_dot*cp_su1*(T_su1-self.T_w))
        
        h_su2 = min(h_max, max(max(h_ex_is, PropsSI('H','Q', 0.1, 'P', P_su1, Fluid)), h_su1 - Q_dot_su/self.m_dot))
        P_su2 = P_su1 #No pressure drop just heat transfer
        
        rho_su2 = PropsSI('D', 'P', P_su2, 'H', h_su2, Fluid)
        T_su2 = PropsSI('T', 'P', P_su2, 'H', h_su2, Fluid)
        s_su2 = PropsSI('S', 'P', P_su2, 'H', h_su2, Fluid)
        
        "4. Leakage"
        try:
            cv_su1 = PropsSI('CVMASS', 'P', P_su1, 'H', h_su1, Fluid)
        except:
            cv_su1 = PropsSI('CVMASS', 'P', P_su1, 'Q', 0, Fluid)
        gamma = max(1e-2, cp_su1/cv_su1)
        P_crit = P_su2*(2/(gamma+1))**(gamma/(gamma-1))
        P_thr_leak = max(P_ex, P_crit)
        rho_thr_leak = PropsSI('D', 'P', P_thr_leak, 'S', s_su2, Fluid)
        h_thr_leak = PropsSI('H', 'P', P_thr_leak, 'S', s_su2, Fluid)
        C_thr_leak = np.sqrt(2*(h_su2-h_thr_leak))
        m_dot_leak = self.params['A_leak']*C_thr_leak*rho_thr_leak
        if self.su.m_dot == None:
            m_dot_in = self.N*self.params['V_s']*rho_su2
            m_dot_leak_bis = self.m_dot-m_dot_in
        elif self.su.m_dot != None:
            m_dot_in = self.m_dot-m_dot_leak
            if self.params['N_rot'] == None:
                self.N = m_dot_in/(self.params['V_s']*rho_su2)
        

        "5. Internal expansion"
        "Isentropic expansion to the internal pressure: su2->in"
        rho_in = rho_su2/self.params['rv_in']
        #Trick to not have problems with CoolProp
        try:
            P_in = PropsSI('P', 'D', rho_in, 'S', s_su2, Fluid)
        except:
            delta = 0.0001
            P_in = 0.5*PropsSI('P', 'D', rho_in*(1+delta), 'S', s_su2, Fluid)+0.5*PropsSI('P', 'D', rho_in*(1-delta), 'S', s_su2, Fluid)
            
        h_in = PropsSI('H', 'D', rho_in, 'P', P_in, Fluid)
        w_in_s = h_su2-h_in
        "Expansion at constant volume: in->ex2"
        w_in_v = (P_in-P_ex)/rho_in
        h_ex2 = h_in-w_in_v
        "Total work"
        W_dot_in = m_dot_in*(w_in_s+w_in_v)
        
        "6. Adiabatic mixing between supply and leakage flows: ex2->ex1"
        h_ex1 = max(min((m_dot_in*h_ex2 + m_dot_leak*h_su2)/self.m_dot, h_su2), h_ex2)
        P_ex1 = P_ex
        T_ex1 = PropsSI('T', 'P', P_ex1, 'H', h_ex1, Fluid)
        try:
            cp_ex2 = PropsSI('CVMASS', 'P', P_ex1, 'H', h_ex1, Fluid)
        except:
            cp_ex2 = PropsSI('CVMASS', 'P', P_ex1, 'Q', 0, Fluid)
        AU_ex = self.params['AU_ex_n']*(self.m_dot/self.params['m_dot_n'])**(0.8)
        C_dot_ex = self.m_dot*cp_ex2
        NTU_ex=AU_ex/C_dot_ex
        epsilon_ex = max(0, 1-np.exp(-NTU_ex))
        Q_dot_ex = max(0, epsilon_ex*C_dot_ex*(self.T_w-T_ex1))
        self.h_ex = h_ex1+Q_dot_ex/self.m_dot

        "8. Energy balance"
        Q_dot_amb = self.params['AU_amb']*(self.T_w-self.inputs['T_amb'])
        W_dot_loss = self.params['alpha']*W_dot_in + self.params['W_dot_loss_0'] + self.params['C_loss']*self.N*2*np.pi
        self.W_dot_exp = W_dot_in - W_dot_loss

        "9. Performances"
        W_dot_s = self.m_dot*(h_su-h_ex_is)
        self.epsilon_is = self.W_dot_exp/W_dot_s
        self.m_dot_th = self.N*self.params['V_s']*rho_su
        self.epsilon_v = self.m_dot/self.m_dot_th
        
        "10. Residuals"
        self.res_E = abs((Q_dot_su + W_dot_loss - Q_dot_ex - Q_dot_amb)/(Q_dot_su + W_dot_loss))
        self.res = self.res_E
        self.res_m_leak = abs((m_dot_leak_bis-m_dot_leak)/m_dot_leak)
        self.res = [self.res, self.res_m_leak]
        return self.res

    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if self.calculable and self.parametrized:
            start_time = time.time()
            ff_guess = [0.7, 1.2, 0.8, 1.3, 0.4, 1.7, 3]
            x_T_guess = [0.7, 0.95, 0.8, 0.99, 0.9]
            stop = 0
            j = 0
            while not stop and j < len(x_T_guess):
                k = 0
                while not stop and k < len(ff_guess):
                    # Loop to permit multiple attempts to solve the implicit calculation 
                    m_dot_guess = ff_guess[k] * self.params['V_s'] * self.inputs['N_rot'] / 60 * PropsSI('D', 'P', self.inputs['su_p'], 'T', self.inputs['su_T'], self.inputs['su_fluid'])
                    T_w_guess = x_T_guess[j] * self.inputs['su_T'] + (1 - x_T_guess[j]) * self.inputs['T_amb']
                    #---------------------------------------------------------------------
                    args = ()
                    x = [m_dot_guess, T_w_guess]
                    #--------------------------------------------------------------------------
                    try:
                        fsolve(self.System, x, args=args)
                        res_norm = np.linalg.norm(self.res)
                    except:
                        res_norm = 1
                    if res_norm < 1e-4:
                        stop = 1
                    k += 1
                j += 1
            self.convergence = stop

            elapsed_time = time.time() - start_time
            if self.convergence == 0:
                print("The system did not converge")
            if self.convergence == 1:
                self.ex.set_fluid(self.su.fluid)
                self.ex.set_m_dot(self.m_dot)
                self.ex.set_h(self.h_ex)
                self.ex.set_p(self.P_ex)
                self.defined = True
        else:
            if not self.calculable:
                print("Input of the component not completely known. Required inputs:")
                for input in self.get_required_inputs():
                    if input not in self.inputs:
                        print(f"  - {input}")
            if not self.parametrized:
                print("Parameters of the component not completely known. Required parameters:")
                for param in self.get_required_parameters():
                    if param not in self.params:
                        print(f"  - {param}")
    
    def print_results(self):
        if self.defined:
            print("=== Expander Results ===")
            print(f"T_ex: {self.ex.T}")
            print(f"W_dot_exp: {self.W_dot_exp}")
            print(f"eta_is: {self.epsilon_is}")
            print(f"eta_v: {self.epsilon_v}")

        else:
            print("Expander component is not defined. Ensure it is solved first.")


