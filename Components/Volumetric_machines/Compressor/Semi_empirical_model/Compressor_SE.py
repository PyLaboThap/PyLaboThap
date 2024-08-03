# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:38:08 2023

@author: Elise
"""


from CoolProp.CoolProp import PropsSI
import numpy as np

from scipy.optimize import fsolve
import time

# TO DO (for the model in general not for the purpose of the off-design model of the CB):
#Find a way so that it can take as input either m_dot or N_rot
# I did something in the code I sent to Olivier but needs to be worked on again


class Compressor_Class:
    def __init__(self):

        self.calculable = False
        self.parametrized = False
        self.defined = False
        
        "Inputs"
        self.su = None
        self.ex = None

        self.N_rot = None
        self.T_amb = None
        
        "Design parameters"
        self.AU_amb = None # global heat transfer coefficient for ambient heat losses     [W/K]
        self.AU_su_n = None # global heat transfer coefficient for supply heat transfer   [W/K]
        self.AU_ex_n = None # global heat transfer coefficient for exhaust heat transfer  	[W/K]
        self.d_ex = None # nozzle diameter for the exhaust pressure losses                [W/K]
        self.m_dot_n = None # global heat transfer coefficient for exhaust heat transfer  	[W/K]
        self.A_leak = None # nozzle cross section area for the leakage                    [m^2]
        self.W_dot_loss_0 = None # constant losses term                                    [W]
        self.alpha = None # proportional losses coefficient                                [-]
        self.C_loss = None # Couple losses                                                  [Nm]
        self.rv_in = None # built-in volumetric ratio                                      [-]
        self.V_s = None # machine displacement volume                                     [m^3]
        
    def inputs(self, su, ex, N_exp, T_amb):

        self.su = su
        self.ex = ex

        self.rp = self.ex.p/self.su.p

        self.N_rot = N_exp
        self.T_amb = T_amb

        self.check_calculable()
    
    def check_calculable(self):
        self.Required_inputs = [self.su.p, self.su.h, self.ex.p, self.N_rot, self.T_amb]
        if all(Input is not None for Input in self.Required_inputs):
                self.calculable = True

    def check_parametrized(self):
        if self.AU_amb != None and self.AU_su_n != None and self.AU_ex_n != None and self.d_ex !=None and self.A_leak != None and self.W_dot_loss_0 != None and self.alpha != None and self.C_loss != None and self.rv_in != None and self.V_s != None:
            self.parametrized = True

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
                    print(f"Warning: Parameter '{key}' not found in the parameters.")

            self.check_parametrized()

    #------------------------------------------------------------------------
    def System(self, x):
            
        "Modelling section of the code"
        
        # If the rotational speed is an input of the model while the mass flow rate is unknown
        self.T_ex2, self.T_w, self.P_ex2, self.m_dot = x
        # print(x)
        self.N = self.N_rot/60
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
        AU_su = self.AU_su_n*(self.m_dot/self.m_dot_n)**0.8
        NTU_su = AU_su/C_dot_su
        epsilon_su = 1-np.exp(-NTU_su)
        Q_dot_su = epsilon_su*C_dot_su*(self.T_w-T_su)
 
        h_su1 = h_su + Q_dot_su/self.m_dot
        P_su1 = P_su
        T_su1 = PropsSI('T', 'H', h_su1, 'P', P_su1, Fluid)
        s_su1 = PropsSI('S', 'H', h_su1, 'P', P_su1, Fluid)
        # print('Q_dot_su', Q_dot_su, 'h_su1', h_su1, 'P_su1', P_su1, 'T_su1', T_su1, 's_su1', s_su1)
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
        
        # print(cv_leak, cp_leak)
        gamma_ex = cp_leak/cv_leak
        P_thr_crit = P_ex*(2/(gamma_ex+1))**(gamma_ex/(gamma_ex-1))
        P_thr = max(P_thr_crit, P_su1)
        # print(P_thr, P_thr_crit)
        s_thr = s_ex2_bis # Isentropic until the throat
        rho_thr = PropsSI('D', 'P', P_thr, 'S', s_thr, Fluid)
        h_thr = PropsSI('H', 'P', P_thr, 'S', s_thr, Fluid)
        C_thr = min(300, np.sqrt(2*(h_ex2_bis-h_thr)))
        V_dot_leak = self.A_leak*C_thr
        m_dot_leak = V_dot_leak*rho_thr
        # print('m_dot_leak', m_dot_leak, 'V_dot_leak', V_dot_leak, 'rho_thr', rho_thr, 'h_thr', h_thr, 'P_thr', P_thr, 's_thr', s_thr)

        "4. Flow rate calculation: su2"
        P_su2 = P_su1
        m_dot_in = m_dot_leak + self.m_dot
        h_su2 = max(min((self.m_dot*h_su1 + m_dot_leak*h_ex2_bis)/m_dot_in, h_ex2_bis), h_su1)
        rho_su2 = PropsSI('D', 'H', h_su2, 'P', P_su2, Fluid)
        s_su2 = PropsSI('S', 'H', h_su2, 'P', P_su2, Fluid)
        self.N_rot_bis = m_dot_in/self.V_s/rho_su2*60
        # print(h_su2, h_ex2_bis, h_su1)

        
        "5. Internal compression: su2->ex2"
        "Isentropic compression: su2->in"
        s_in = s_su2
        rho_in = rho_su2*self.rv_in
        P_in = PropsSI('P', 'D', rho_in, 'S', s_in, Fluid)
        h_in = PropsSI('H', 'D', rho_in, 'P', P_in, Fluid)
        v_in = 1/rho_in
        rp_in = P_in/P_su2

        w_in_is = h_in-h_su2
        
        "Isochoric compression: in->ex2"
        # rho_ex2 = rho_in
        # T_ex2 = PropsSI('T', 'P', self.P_ex2, 'D', rho_ex2, Fluid)
        # h_ex2 = PropsSI('H', 'P', self.P_ex2, 'D', rho_ex2, Fluid)
        # s_ex2 = PropsSI('S', 'P', self.P_ex2, 'D', rho_ex2, Fluid)
        # w_in_v = v_in*(self.P_ex2-P_in)
        w_in_v = (self.P_ex2-P_in)/rho_in#*1000
        
        "Total internal work"
        w_in = w_in_is + w_in_v
        h_ex2 = h_su2 + w_in
        s_ex2 = PropsSI('S', 'P', self.P_ex2, 'H', h_ex2, Fluid)
        T_ex2 = PropsSI('T', 'P', self.P_ex2, 'H', h_ex2, Fluid)


        "6. Pressure drops: ex2->ex1"
        A_ex = np.pi*(self.d_ex/2)**2
        
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
        AU_ex = self.AU_ex_n*(self.m_dot/self.m_dot_n)**0.8
        
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
        Q_dot_amb = self.AU_amb*(self.T_w-self.T_amb)
        
        
        "Compression work and power"
        W_dot_in = m_dot_in*w_in
        W_dot_loss = self.alpha*W_dot_in + self.W_dot_loss_0 + self.C_loss*self.N*2*np.pi
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
        V_s_dot = self.V_s*self.N
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
                    m_dot_guess = x_m_guess[k]*self.V_s*self.N_rot/60*self.su.D #ff_guess[k]*self.V_s*self.N_rot/60*PropsSI('D', 'P', self.su.p, 'H', self.su.h, self.su.fluid) #initial value for M_dot
                    T_ex2_guess = PropsSI('T','P', self.ex.p,'S', self.su.s, self.su.fluid)+5 #PropsSI('T', 'P', self.su.p*self.rp,'S', self.su.s, self.su.fluid)
                    P_ex2_guess = 0.9*self.su.p*self.rp
                    #---------------------------------------------------------------------
                    args = ()
                    x = [T_ex2_guess, T_w_guess, P_ex2_guess, m_dot_guess]
                    #--------------------------------------------------------------------------
                    # ub = 2*x # upper bound for fsolve
                    try:
                        fsolve(self.System, x, args = args)
                        res_norm = np.linalg.norm(self.res)
                    except:
                        res_norm = 1
                        pass
                    # print(res_norm, self.res)
                    if res_norm < 1e-3:
                        stop = 1
                    k = k + 1
                j = j + 1

                self.convergence = stop
                
            
            self.N_rot = self.N*60
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