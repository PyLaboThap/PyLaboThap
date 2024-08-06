# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:38:08 2023

@author: Elise
"""

from components.base_component import BaseComponent
from connectors.mass_connector import MassConnector
from connectors.work_connector import WorkConnector
from connectors.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import numpy as np
import time

# TO DO (for the model in general not for the purpose of the off-design model of the CB):
#Find a way so that it can take as input either m_dot or N_rot

class CompressorSE(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector() # Mass_connector
        self.work_cp = WorkConnector()
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
            if self.work_cp.speed is not None:
                self.inputs['N_rot'] = self.work_cp.speed
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
                self.work_cp.set_speed(self.inputs['N_rot'])
            if 'T_amb' in self.inputs:
                self.heat_amb.set_temperature_in(self.inputs['T_amb'])

        return ['su_p', 'su_T', 'ex_p', 'N_rot', 'T_amb', 'su_fluid']
    
    def get_required_parameters(self):
        return [
            'AU_amb', 'AU_su_n', 'AU_ex_n', 'd_ex', 'm_dot_n', 
            'A_leak', 'W_dot_loss_0', 'alpha', 'C_loss', 'rv_in', 'V_s'
        ]

    def print_setup(self):
        print("=== Compressor Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")
        print(f"  - W_dot: speed={self.work_cp.speed}")
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

    #------------------------------------------------------------------------
    def System(self, x):
            
        "Modelling section of the code"
        
        # If the rotational speed is an input of the model while the mass flow rate is unknown
        self.T_ex2, self.T_w, self.P_ex2, self.m_dot = x
        # print(x)
        self.N = self.inputs['N_rot']/60
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

        T_sat_su = PropsSI('T', 'P', P_su, 'Q', 1, Fluid)
        if T_su<T_sat_su:
            print('----Warning the compressor inlet stream is not in vapor phase---')

        "2. Supply heating: su->su1"
        try:
            cp_su = PropsSI('C', 'P', P_su, 'H', h_su, Fluid)
        except:
            cp_su = PropsSI('C', 'P', P_su, 'Q', 1, Fluid)
            
        C_dot_su = self.m_dot*cp_su
        AU_su = self.params['AU_su_n']*(self.m_dot/self.params['m_dot_n'])**0.8
        NTU_su = AU_su/C_dot_su
        epsilon_su = 1-np.exp(-NTU_su)
        Q_dot_su = epsilon_su*C_dot_su*(self.T_w-T_su)
 
        h_su1 = h_su + Q_dot_su/self.m_dot
        P_su1 = P_su
        T_su1 = PropsSI('T', 'H', h_su1, 'P', P_su1, Fluid)
        s_su1 = PropsSI('S', 'H', h_su1, 'P', P_su1, Fluid)
       
        #------------------------------------------------------------------
        
        "3. Internal leakage: ex2->su1"
        s_ex2_bis = PropsSI('S', 'T', self.T_ex2, 'P', self.P_ex2, Fluid) # C ICI LE PB??
        h_ex2_bis = PropsSI('H', 'T', self.T_ex2, 'P', self.P_ex2, Fluid)
        
        try: 
            cv_leak = PropsSI('CVMASS','P', self.P_ex2, 'H', h_ex2_bis, Fluid)
            cp_leak = PropsSI('CPMASS','P', self.P_ex2, 'H', h_ex2_bis, Fluid)
        except:
            cv_leak = PropsSI('CVMASS', 'P', self.P_ex2, 'Q', 0, Fluid)
            cp_leak = PropsSI('CPMASS', 'P', self.P_ex2, 'Q', 0, Fluid)
        
        
        gamma_ex = cp_leak/cv_leak
        P_thr_crit = P_ex*(2/(gamma_ex+1))**(gamma_ex/(gamma_ex-1))
        P_thr = max(P_thr_crit, P_su1)
        # print(P_thr, P_thr_crit)
        s_thr = s_ex2_bis # Isentropic until the throat
        rho_thr = PropsSI('D', 'P', P_thr, 'S', s_thr, Fluid)
        h_thr = PropsSI('H', 'P', P_thr, 'S', s_thr, Fluid)
        C_thr = min(300, np.sqrt(2*(h_ex2_bis-h_thr)))
        V_dot_leak = self.params['A_leak']*C_thr
        m_dot_leak = V_dot_leak*rho_thr

        "4. Flow rate calculation: su2"
        P_su2 = P_su1
        m_dot_in = m_dot_leak + self.m_dot
        h_su2 = max(min((self.m_dot*h_su1 + m_dot_leak*h_ex2_bis)/m_dot_in, h_ex2_bis), h_su1)
        rho_su2 = PropsSI('D', 'H', h_su2, 'P', P_su2, Fluid)
        s_su2 = PropsSI('S', 'H', h_su2, 'P', P_su2, Fluid)
        self.N_rot_bis = m_dot_in/self.params['V_s']/rho_su2*60
        

        
        "5. Internal compression: su2->ex2"
        "Isentropic compression: su2->in"
        s_in = s_su2
        rho_in = rho_su2*self.params['rv_in']
        P_in = PropsSI('P', 'D', rho_in, 'S', s_in, Fluid)
        h_in = PropsSI('H', 'D', rho_in, 'P', P_in, Fluid)
        v_in = 1/rho_in
        rp_in = P_in/P_su2

        w_in_is = h_in-h_su2
        
        "Isochoric compression: in->ex2"

        w_in_v = (self.P_ex2-P_in)/rho_in#*1000
        
        "Total internal work"
        w_in = w_in_is + w_in_v
        h_ex2 = h_su2 + w_in
        s_ex2 = PropsSI('S', 'P', self.P_ex2, 'H', h_ex2, Fluid)
        T_ex2 = PropsSI('T', 'P', self.P_ex2, 'H', h_ex2, Fluid)

        
        "6. Pressure drops: ex2->ex1"
        A_ex = np.pi*(self.params['d_ex']/2)**2
        
        h_ex1 = h_ex2 #Isenthalpic valve
        
        P_crit_ex = self.P_ex2*(2/(gamma_ex+1))**(gamma_ex/(gamma_ex-1))
        P_thr_ex = max(P_crit_ex, P_ex)
        h_thr_ex = PropsSI('H', 'P', P_thr_ex, 'S', s_ex2, Fluid)
        v_thr_ex = 1./PropsSI('D', 'P', P_thr_ex, 'S', s_ex2, Fluid)
        
        V_dot_ex = self.m_dot*v_thr_ex
        C_ex = V_dot_ex/A_ex
        h_ex1_bis = h_thr_ex + (C_ex**2)/2
        
        T_ex1 = PropsSI('T', 'H', h_ex1, 'P', P_ex, Fluid)
        
        "7. Cooling at exit: ex1 -> ex "
        AU_ex = self.params['AU_ex_n']*(self.m_dot/self.params['m_dot_n'])**0.8
        
        try:
            cp_ex = PropsSI('C', 'P', P_ex, 'H', h_ex1, Fluid)
        except:
            cp_su = PropsSI('C', 'P', P_ex, 'Q', 1, Fluid)
        
        C_dot_ex = cp_ex*self.m_dot
        NTU_ex = AU_ex/C_dot_ex
        epsilon_ex = 1-np.exp(-NTU_ex)
        Q_dot_ex = epsilon_ex*C_dot_ex*(self.T_w-T_ex1)
        
        h_ex = h_ex1+(Q_dot_ex/self.m_dot)
        
        "Fictious enveloppe heat balance"
        Q_dot_amb = self.params['AU_amb']*(self.T_w-self.inputs['T_amb'])
        
        
        "Compression work and power"
        W_dot_in = m_dot_in*w_in
        W_dot_loss = self.params['alpha']*W_dot_in + self.params['W_dot_loss_0'] + self.params['C_loss']*self.N*2*np.pi
        self.W_dot = W_dot_in + W_dot_loss
        
        "Exit data"
        self.h_ex = h_ex
        self.P_ex = P_ex
        self.T_ex = PropsSI('T', 'H', h_ex, 'P', P_ex, Fluid)
        
        "Isentropic efficiency"
        h_ex_is = PropsSI('H', 'S', s_su, 'P', P_ex, Fluid)
        w_s = h_ex_is-h_su
        W_dot_s = self.m_dot*w_s
        self.epsilon_is = W_dot_s/self.W_dot
        
        "Volumetric efficiency"
        #Theoretical flowrate
        V_s_dot = self.params['V_s']*self.N
        m_dot_th = V_s_dot*rho_su
        
        m_dot_in_bis = V_s_dot*rho_su2
        
        #Volumetric efficiencies definitions
        self.epsilon_v = self.m_dot/m_dot_th
        self.epsilon_v_l = self.m_dot/m_dot_in
        self.epsilon_v_PT = m_dot_in/m_dot_th
        
        "Residue"
        
        self.res_h_ex1 = abs(h_ex1_bis-h_ex1)/h_ex1
        self.resE = abs((W_dot_loss - Q_dot_ex - Q_dot_su - Q_dot_amb)/(W_dot_loss))
        self.res_h_ex2 = abs(h_ex2_bis-h_ex2)/h_ex2
        self.res_m_dot_in = abs(m_dot_in-m_dot_in_bis)/m_dot_in
        self.res = [self.res_h_ex1, self.resE, self.res_h_ex2, self.res_m_dot_in]
        # print(W_dot_loss, Q_dot_ex, Q_dot_su, Q_dot_amb)
        "problem with entropy"
        self.h_su = h_su
        self.h_su1 = h_su1
        self.h_su2 = h_su2
        self.h_in = h_in
        self.h_ex2 = h_ex2
        self.h_ex1 = h_ex1
        self.h_ex = h_ex

        self.s_su = s_su
        self.s_su1 = s_su1
        self.s_su2 = s_su2
        self.s_in = s_in
        self.s_ex2 = s_ex2
        # self.s_ex1 = s_ex2
        self.s_ex = PropsSI('S', 'T', self.T_ex, 'P', P_ex, Fluid)
        # print(self.s_su, self.s_su1, self.s_su2, self.s_in, self.s_ex2, self.s_ex)

        # print(self.res)
        return self.res
    
                
    #------------------------------------------------------------------------
    def solve(self):
        self.check_calculable()
        self.check_parametrized()
        
        if self.calculable and self.parametrized:
            
            start_time = time.time()
            
            x_m_guess = [1.1, 0.99, 1.3, 1, 1.15] #guesses on the filling factor to provide suitable initial point for the iteration
            x_T_guess = [0.9, 1.01, 0.7, 1.1, 0.2] #For the iteration on the T_w
            stop = 0
            
            #If N_rot is known
            j = 0
            while not stop and j < len(x_T_guess):
                k = 0
                while not stop and k < len(x_m_guess):
                        
                    # Loop to permit multiple attempts to solve the implicit calculation
                    T_w_guess = x_T_guess[j]*PropsSI('T', 'P', self.ex.p, 'S', self.su.s, self.su.fluid)+5 #x_T_guess[j]*self.su.T+(1-x_T_guess[j])*self.T_amb
                    m_dot_guess = x_m_guess[k]*self.params['V_s']*self.inputs['N_rot']/60*self.su.D #ff_guess[k]*self.V_s*self.N_rot/60*PropsSI('D', 'P', self.su.p, 'H', self.su.h, self.su.fluid) #initial value for M_dot
                    T_ex2_guess = PropsSI('T','P', self.ex.p,'S', self.su.s, self.su.fluid)+5 #PropsSI('T', 'P', self.su.p*self.rp,'S', self.su.s, self.su.fluid)
                    P_ex2_guess = 0.9*self.ex.p
                    #---------------------------------------------------------------------
                    args = ()
                    x = [T_ex2_guess, T_w_guess, P_ex2_guess, m_dot_guess]
                    #--------------------------------------------------------------------------
                    # ub = 2*x # upper bound for fsolve
                    try:
                        fsolve(self.System, x, args = args)
                        res_norm = np.linalg.norm(self.res)
                        print(res_norm)
                    except:
                        res_norm = 1
                        pass
                    # print(res_norm, self.res)
                    if res_norm < 1e-3:
                        stop = 1
                    k = k + 1
                j = j + 1

                self.convergence = stop
                
            self.ex.set_fluid(self.su.fluid)
            self.ex.set_m_dot(self.m_dot)
            self.ex.set_h(self.h_ex)
            self.ex.set_p(self.P_ex)

            self.defined = True
            
            elapsed_time = time.time() - start_time
            #print("Optimization time:", elapsed_time, "seconds")
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")

    def print_results(self):
        if self.defined:
            print("=== Compressor Results ===")
            print(f"T_ex: {self.ex.T}")
            print(f"W_dot_cp: {self.W_dot}")
            print(f"eta_is: {self.epsilon_is}")
            print(f"eta_v: {self.epsilon_v}")

        else:
            print("Expander component is not defined. Ensure it is solved first.")