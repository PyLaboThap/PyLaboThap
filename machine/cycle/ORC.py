from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.heat_exchanger.pinch_cst.simulation_model import HXPinchCst
from component.heat_exchanger.moving_boundary.charge_sensitive.simulation_model import HeatExchangerMB
from component.heat_exchanger.moving_boundary.charge_sensitive.modules.geometry_plate_hx_swep import PlateGeomSWEP
from component.volumetric_machine.expander.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from component.pump.constant_efficiency.simulation_model import PumpCstEff
from machine.cycle.cycle import Cycle



from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize
from scipy.optimize import fsolve


from scipy.optimize import minimize
import numpy as np


# class RC(Cycle):
#     def __init__(self, fluid=None):
#         super().__init__(fluid)
#         self.tolerance = 1e-6  # Convergence tolerance for residuals
#         self.prev_pressures = {}  # Store pressures between iterations

#     def solve(self):
#         self.set_cycle_guesses()

#         iteration = 0
#         max_iterations = 100
#         converged = False

#         while not converged and iteration < max_iterations:
#             print(f"Iteration {iteration}: Solving the cycle.")
            
#             # Calculate initial residuals based on guesses or previous iteration
#             initial_residuals = self.get_residuals()

#             # Recursively solve the cycle starting from the pump
#             visited = set()  # Reset visited components for each iteration
#             self.recursive_solve(self.get_component("Pump"), visited)

#             # Calculate final residuals after solving
#             final_residuals = self.get_residuals()

#             # Check if the residuals are within the tolerance
#             converged = self.check_residuals(initial_residuals, final_residuals)

#             # Store calculated pressures for the next iteration
#             self.store_pressures()

#             # Update guesses for the next iteration
#             self.set_cycle_guesses()

#             iteration += 1

#         if converged:
#             print("Cycle solved successfully.")
#         else:
#             print("Cycle did not converge within the maximum number of iterations.")

#     def set_cycle_guesses(self):
#         # Update pressures with previous iteration values if available, otherwise use initial guesses
#         P_ev_guess = self.prev_pressures.get('P_ev', PropsSI('P', 'T', 90 + 273.15, 'Q', 0.5, self.fluid))
#         P_cd_guess = self.prev_pressures.get('P_cd', PropsSI('P', 'T', 25 + 273.15, 'Q', 0.5, self.fluid))

#         subcooling = self.parameters.get('SC_cd', 5)
#         superheating = self.parameters.get('SH_ev', 5)

#         T_su_pp = PropsSI('T', 'P', P_cd_guess, 'Q', 0, self.fluid) - subcooling
#         self.set_cycle_properties(target="Pump:su", T=T_su_pp, P=P_cd_guess)

#         self.set_cycle_properties(target="Pump:ex", P=P_ev_guess)
#         T_ex_ev = PropsSI('T', 'P', P_ev_guess, 'Q', 1, self.fluid) + superheating
#         self.set_cycle_properties(target="Evaporator:ex_wf", P=P_ev_guess, T=T_ex_ev)
#         self.set_cycle_properties(target="Expander:ex", P=P_cd_guess)

#     def store_pressures(self):
#         # Store pressures from the evaporator and condenser models
#         evaporator = self.get_component("Evaporator")
#         condenser = self.get_component("Condenser")

#         self.prev_pressures['P_ev'] = evaporator.model.ex_wf.p
#         self.prev_pressures['P_cd'] = condenser.model.ex_wf.p


#     def recursive_solve(self, component, visited):
#         # Prevent infinite loops by skipping already visited components
#         if component in visited:
#             return

#         # Mark the component as visited and solve it
#         visited.add(component)
#         print(f"Solving {component.name}")
#         component.solve()

#         # Recursively solve connected components
#         for next_component in component.next.values():
#             self.recursive_solve(next_component, visited)

#     def get_residuals(self):
#         # Retrieve components that define the residuals
#         evaporator = self.get_component("Evaporator")
#         condenser = self.get_component("Condenser")
#         pump = self.get_component("Pump")


#         # Calculate residuals based on enthalpies
#         h_ex_ev = evaporator.model.ex_wf.h
#         h_ex_cd = condenser.model.ex_wf.h

#         return h_ex_ev, h_ex_cd

#     def check_residuals(self, residuals_i, residuals_f):
#         # Calculate the difference between initial and final residuals
#         residual_h_ex_ev = abs(residuals_i[0] - residuals_f[0])
#         residual_h_ex_cd = abs(residuals_i[1] - residuals_f[1])

#         # Output residuals for debugging
#         print(f"Residual Evaporator: {residual_h_ex_ev}, Residual Condenser: {residual_h_ex_cd}")

#         # Check if residuals are within the tolerance
#         return residual_h_ex_ev < self.tolerance and residual_h_ex_cd < self.tolerance
    
#     # def set_cycle_guesses_residuals(self, guesses, residuals_var):
#     #     """
#     #     Set the initial guesses and calculate related properties.
#     #     Also, set the target properties used to compute the residuals
#     #     """
#     #     for guess, value in guesses.items():
#     #         comp_port, prop = guess.split('-')
#     #         comp, port = comp_port.split(':')
#     #         fluid = self.components[comp].fluid
#     #         if prop == 'SC':
#     #             P_sat_cd = self.components['Condenser'].model.ex_H.p
#     #             if P_sat_cd is None:
#     #                 P_sat_cd = self.components['Condenser'].model.su_H.p
#     #             T_cd = PropsSI('T', 'P', P_sat_cd, 'Q', 0, fluid) - abs(value)
#     #             kwargs = {'target' : comp_port,
#     #                     'T' : T_cd,
#     #                     'fluid' : fluid}

#     #             self.set_guess(**kwargs)    
#     #         elif prop == 'SH':
#     #             P_sat_ev = self.components['Evaporator'].model.ex_C.p
#     #             if P_sat_ev is None:
#     #                 P_sat_ev = self.components['Evaporator'].model.su_C.p
#     #             T_ev = PropsSI('T', 'P', P_sat_ev, 'Q', 1, fluid) + abs(value)
#     #             kwargs = {'target' : comp_port,
#     #                     'T' : T_ev,
#     #                     'fluid' : fluid}

#     #             self.set_guess(**kwargs)    

#     #         else:
#     #             kwargs = {'target' : comp_port,
#     #                     prop : value,
#     #                     'fluid' : fluid}
#     #             self.set_guess(**kwargs)
#     #     for residual_target in residuals_var:
#     #         self.residuals_var[residual_target] = None
#     #         return


# # Example usage of the Rankine Cycle (RC)
# if __name__ == "__main__":
#     orc_cycle = RC(fluid='R245fa')

#     # Add components
#     PUMP = PumpCstEff()
#     EVAP = HXPinchCst()
#     EXP = ExpanderCstEff()
#     COND = HXPinchCst()

#     # Set component parameters
#     PUMP.set_parameters(eta_is=0.6)
#     EVAP.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='evaporator')
#     COND.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='condenser')
#     EXP.set_parameters(eta_is=0.8)

#     # Add components to the cycle
#     orc_cycle.add_component(PUMP, "Pump")
#     orc_cycle.add_component(EVAP, "Evaporator")
#     orc_cycle.add_component(EXP, "Expander")
#     orc_cycle.add_component(COND, "Condenser")

#     # Link components
#     orc_cycle.link_components("Pump", "m-ex", "Evaporator", "m-su_wf")
#     orc_cycle.link_components("Evaporator", "m-ex_wf", "Expander", "m-su")
#     orc_cycle.link_components("Expander", "m-ex", "Condenser", "m-su_wf")
#     orc_cycle.link_components("Condenser", "m-ex_wf", "Pump", "m-su")

#     # Set the cycle properties
#     orc_cycle.set_cycle_properties(m_dot=0.06, target='Pump:su')

#     orc_cycle.set_cycle_properties(T=30 + 273.15, fluid='Water', m_dot=0.4, target='Condenser:su_sf')
#     orc_cycle.set_cycle_properties(cp=4186, target='Condenser:su_sf')
#     orc_cycle.set_cycle_properties(fluid='Water', target='Condenser:ex_sf')

#     orc_cycle.set_cycle_properties(T=150 + 273.15, fluid='Water', m_dot=0.4, target='Evaporator:su_sf')
#     orc_cycle.set_cycle_properties(cp=4186, target='Evaporator:su_sf')
#     orc_cycle.set_cycle_properties(fluid='Water', target='Evaporator:ex_sf')

#     # Set parameters for the cycle
#     SC_cd = 5
#     SH_ev = 5
#     orc_cycle.set_cycle_parameters(SC_cd=SC_cd, SH_ev=SH_ev)

#     # Initial guesses for pressures
#     T_ev_guess = 120 + 273.15
#     P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R245fa')

#     T_cd_guess = 3 + 273.15
#     P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, 'R245fa')

#     # Define pressure guesses
#     pressure_guesses = {
#         "P_ev": P_ev_guess,
#         "P_cd": P_cd_guess
#     }

#     guesses = {
#     "Evaporator:su_C-P": P_ev_guess,
#     "Condenser:ex_H-P": P_cd_guess,
#     "Condenser:ex_H-SC": SC_cd
#     }

#     # Define the residuals to converge on
#     residuals_var = [
#         ("Evaporator:ex_C-h"),
#         ("Condenser:ex_H-h")
#     ]
 
#     orc_cycle.set_cycle_guesses_residuals(guesses,residuals_var)

#     # Solve the cycle
#     orc_cycle.solve()

    # Optional: Print a summary of the cycle
    # orc_cycle.print_summary()

class RC(Cycle):
    def __init__(self, fluid=None):
        super().__init__(fluid)
        self.tolerance = 1e-6  # Convergence tolerance for residuals
        self.prev_pressures = {}  # Store pressures between iterations
        self.prev_residuals = None  # Store residuals from the previous iteration
        self.residuals_var = []  # Store the residuals to check

    def solve(self, start_key="Pump"):
        self.set_cycle_guesses()

        iteration = 0
        max_iterations = 100
        converged = False

        while not converged and iteration < max_iterations:
            print(f"Iteration {iteration}: Solving the cycle.")

            # Calculate residuals before solving
            self.prev_residuals = self.get_residuals()

            # Recursively solve the cycle starting from the pump
            visited = set()  # Reset visited components for each iteration
            self.recursive_solve(self.get_component(start_key), visited)

            # Calculate residuals after solving
            final_residuals = self.get_residuals()

            # Update guesses for the next iteration
            self.update_guesses()

            # Check if the residuals are within the tolerance
            converged = self.check_residuals(final_residuals)

            iteration += 1

        if converged:
            print("Cycle solved successfully.")
        else:
            print("Cycle did not converge within the maximum number of iterations.")

    def set_cycle_guesses_residuals(self, guesses, residuals_var):
        """
        Set the initial guesses and define the variables used to compute the residuals.
        """
        self.residuals_var = residuals_var
        print(f"Residuals variables: {self.residuals_var}")

        # Set initial guesses for the cycle based on provided values
        for guess, value in guesses.items():
            component_name, prop = guess.split(':')
            connector_name, prop_name = prop.split('-')
            self.set_cycle_guess(target=f"{component_name}:{connector_name}", **{prop_name: value})

    def update_guesses(self):
        """
        Update the guesses for the next iteration based on the pressures from the current iteration.
        """
        evaporator = self.get_component("Evaporator")
        condenser = self.get_component("Condenser")

        self.prev_pressures['P_ev'] = evaporator.model.ex_C.p
        self.prev_pressures['P_cd'] = condenser.model.ex_H.p

        # Update the guesses with the latest pressure values
        self.guesses.update({
            "Evaporator:su_C-P": self.prev_pressures['P_ev'],
            "Condenser:ex_H-P": self.prev_pressures['P_cd']
        })

    def check_residuals(self, final_residuals):
        """
        Check if the difference between the previous and current residuals is within the tolerance.
        """
        if not self.prev_residuals:
            return False  # No previous residuals to compare to

        residual_diff = [abs(f - p) for f, p in zip(final_residuals, self.prev_residuals)]
        
        # Output residuals for debugging
        for i, diff in enumerate(residual_diff):
            print(f"Residual {self.residuals_var[i]}: {diff}")

        return all(diff < self.tolerance for diff in residual_diff)

    def get_residuals(self):
        """
        Calculate the residuals based on the specified variables.
        """
        residuals = [] # Store the calculated residuals
        for residual_target in self.residuals_var: # Iterate over the residuals to calculate
            component_name, connector_prop = residual_target.split(':') # Split the target into component and connector
            connector_name, prop_name = connector_prop.split('-') # Split the connector into name and property
            component = self.get_component(component_name) # Get the component object
            residual_value = getattr(getattr(component.model, connector_name), prop_name) # Get the value of the property from the connector
            residuals.append(residual_value) # Append the calculated residual to the list
            print(f"Residual {residual_target}: {residual_value}") # Output the calculated residual

        return residuals

    def set_cycle_guesses(self):
        # This method is now mainly for setting up the initial conditions, with updates handled in `update_guesses`.
        P_ev_guess = self.prev_pressures.get('P_ev', PropsSI('P', 'T', 90 + 273.15, 'Q', 0.5, self.fluid))
        P_cd_guess = self.prev_pressures.get('P_cd', PropsSI('P', 'T', 25 + 273.15, 'Q', 0.5, self.fluid))

        subcooling = self.parameters['SC_cd']
        superheating = self.parameters['SH_ev']

        T_su_pp = PropsSI('T', 'P', P_cd_guess, 'Q', 0, self.fluid) - subcooling
        self.set_cycle_properties(target="Pump:su", T=T_su_pp, P=P_cd_guess)

        T_ex_ev = PropsSI('T', 'P', P_ev_guess, 'Q', 1, self.fluid) + superheating
        self.set_cycle_properties(target="Evaporator:ex_C", P=P_ev_guess, T=T_ex_ev)
        self.set_cycle_properties(target="Expander:ex", P=P_cd_guess)

    def recursive_solve(self, component, visited):
        # Prevent infinite loops by skipping already visited components
        if component in visited:
            return

        # Mark the component as visited and solve it
        visited.add(component)
        print(f"Solving {component.name}")
        if isinstance(component, Cycle.Component):
            component.solve()

        # Recursively solve connected components
        for next_component in component.next.values():
            self.recursive_solve(next_component, visited)


# Example usage of the Rankine Cycle (RC)
# if __name__ == "__main__":
#     orc_cycle = RC(fluid='R245fa')

#     # Add components
#     PUMP = PumpCstEff()
#     EVAP = HXPinchCst()
#     EXP = ExpanderCstEff()
#     COND = HXPinchCst()

#     # Set component parameters
#     PUMP.set_parameters(eta_is=0.6)
#     EVAP.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='evaporator')
#     COND.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='condenser')
#     EXP.set_parameters(eta_is=0.8)

#     # Add components to the cycle
#     orc_cycle.add_component(PUMP, "Pump")
#     orc_cycle.add_component(EVAP, "Evaporator")
#     orc_cycle.add_component(EXP, "Expander")
#     orc_cycle.add_component(COND, "Condenser")

#     # Link components
#     orc_cycle.link_components("Pump", "m-ex", "Evaporator", "m-su_wf")
#     orc_cycle.link_components("Evaporator", "m-ex_wf", "Expander", "m-su")
#     orc_cycle.link_components("Expander", "m-ex", "Condenser", "m-su_wf")
#     orc_cycle.link_components("Condenser", "m-ex_wf", "Pump", "m-su")

#     # Set the cycle properties
#     orc_cycle.set_cycle_properties(m_dot=0.06, target='Pump:su')

#     orc_cycle.set_cycle_properties(T=30 + 273.15, fluid='Water', m_dot=0.4, target='Condenser:su_sf')
#     orc_cycle.set_cycle_properties(cp=4186, target='Condenser:su_sf')
#     orc_cycle.set_cycle_properties(fluid='Water', target='Condenser:ex_sf')

#     orc_cycle.set_cycle_properties(T=150 + 273.15, fluid='Water', m_dot=0.4, target='Evaporator:su_sf')
#     orc_cycle.set_cycle_properties(cp=4186, target='Evaporator:su_sf')
#     orc_cycle.set_cycle_properties(fluid='Water', target='Evaporator:ex_sf')

#     # Set parameters for the cycle
#     SC_cd = 5
#     SH_ev = 5
#     orc_cycle.set_cycle_parameters(SC_cd=SC_cd, SH_ev=SH_ev)

#     # Initial guesses for pressures
#     T_ev_guess = 120 + 273.15
#     P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R245fa')

#     T_cd_guess = 3 + 273.15
#     P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, 'R245fa')

#     # Define guesses and residuals
#     guesses = {
#         "Evaporator:su_wf-P": P_ev_guess,
#         "Condenser:ex_wf-P": P_cd_guess,
#     }

#     residuals_var = [
#         "Evaporator:ex_wf-h",
#         "Condenser:ex_wf-h"
#     ]

#     orc_cycle.set_cycle_guesses_residuals(guesses, residuals_var)

#     # Solve the cycle
#     orc_cycle.solve()

    # Example usage of ORC
if __name__ == "__main__":
    orc_cycle = RC(fluid='Cyclopentane')
 
    #%%
    "EVAPORATOR PARAMETERS"
    EVAP = HeatExchangerMB('Plate')
    "Geometry Loading"
    EVAP_geom = PlateGeomSWEP()
    EVAP_geom.set_parameters("B20Hx24/1P") 
    Corr_H = {"1P" : "Gnielinski", "2P" : "Han_cond_BPHEX"}
    Corr_C = {"1P" : "Gnielinski", "2P" : "Han_Boiling_BPHEX_HTC"}
 
    EVAP.set_parameters(
        A_c = EVAP_geom.A_c, A_h = EVAP_geom.A_h, h = EVAP_geom.h, l = EVAP_geom.l, l_v = EVAP_geom.l_v, # 5
        C_CS = EVAP_geom.C_CS, C_Dh = EVAP_geom.C_Dh, C_V_tot = EVAP_geom.C_V_tot, C_canal_t = EVAP_geom.C_canal_t, C_n_canals = EVAP_geom.C_n_canals, # 10
        H_CS = EVAP_geom.H_CS, H_Dh = EVAP_geom.H_Dh, H_V_tot = EVAP_geom.H_V_tot, H_canal_t = EVAP_geom.H_canal_t, H_n_canals = EVAP_geom.H_n_canals, # 15
        casing_t = EVAP_geom.casing_t, chevron_angle = EVAP_geom.chevron_angle, fooling = EVAP_geom.fooling, # 18
        n_plates = EVAP_geom.n_plates, plate_cond = EVAP_geom.plate_cond, plate_pitch_co = EVAP_geom.plate_pitch_co, t_plates = EVAP_geom.t_plates, w = EVAP_geom.w, # 23
        Flow_Type = 'CounterFlow', H_DP_ON = True, C_DP_ON = True, n_disc = 50) # 27
    EVAP.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 28
 
    orc_cycle.add_component(EVAP, "Evaporator")
 
    #%%
    "PUMP PARAMETERS"
 
    PUMP = PumpCstEff()
 
    # Set parameters
    PUMP.set_parameters(eta_is=0.6)
    orc_cycle.add_component(PUMP, "Pump")
    #%%
    "EXPANDER PARAMETERS"
    EXP = ExpanderCstEff()
    EXP.set_parameters(eta_is=0.8)
 
    orc_cycle.add_component(EXP, "Expander")
    #%%
    "CONDENSER PARAMETERS"
    COND = HeatExchangerMB('Plate')
    "Geometry Loading"
    COND_geom = PlateGeomSWEP()
    COND_geom.set_parameters("B20Hx24/1P")
    Corr_H = {"1P" : "Gnielinski", "2P" : "Han_cond_BPHEX"}
    Corr_C = {"1P" : "Gnielinski", "2P" : "Han_Boiling_BPHEX_HTC"}
 
    COND.set_parameters(
        A_c = EVAP_geom.A_c, A_h = EVAP_geom.A_h, h = EVAP_geom.h, l = EVAP_geom.l, l_v = EVAP_geom.l_v, # 5
        C_CS = EVAP_geom.C_CS, C_Dh = EVAP_geom.C_Dh, C_V_tot = EVAP_geom.C_V_tot, C_canal_t = EVAP_geom.C_canal_t, C_n_canals = EVAP_geom.C_n_canals, # 10
        H_CS = EVAP_geom.H_CS, H_Dh = EVAP_geom.H_Dh, H_V_tot = EVAP_geom.H_V_tot, H_canal_t = EVAP_geom.H_canal_t, H_n_canals = EVAP_geom.H_n_canals, # 15
        casing_t = EVAP_geom.casing_t, chevron_angle = EVAP_geom.chevron_angle, fooling = EVAP_geom.fooling, # 18
        n_plates = EVAP_geom.n_plates, plate_cond = EVAP_geom.plate_cond, plate_pitch_co = EVAP_geom.plate_pitch_co, t_plates = EVAP_geom.t_plates, w = EVAP_geom.w, # 23
        Flow_Type = 'CounterFlow', H_DP_ON = True, C_DP_ON = True, n_disc = 50) # 27
    COND.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 28
 
    orc_cycle.add_component(COND, "Condenser")
 
    #%% Sources and sinks
    Source_T66 = orc_cycle.Source('T66', orc_cycle.get_component('Evaporator'), 'm-su_H')
    Sink_T66 = orc_cycle.Sink('T66', orc_cycle.get_component('Evaporator'), 'm-ex_H')
    # orc_cycle.set_solve_start('T66', Source_T66)
 
    Source_Water = orc_cycle.Source('Water', orc_cycle.get_component('Condenser'), 'm-su_C')    
    Sink_Water = orc_cycle.Sink('Water', orc_cycle.get_component('Condenser'), 'm-ex_C')    
    # orc_cycle.set_solve_start('Water', Source_Water)
 
    #%%    
 
    # Link components
    orc_cycle.link_components("Pump", "m-ex", "Evaporator", "m-su_C")
    orc_cycle.link_components("Evaporator", "m-ex_C", "Expander", "m-su")
    orc_cycle.link_components("Expander", "m-ex", "Condenser", "m-su_H")
    orc_cycle.link_components("Condenser", "m-ex_H", "Pump", "m-su")
 
    # Set the inputs using set_properties
    orc_cycle.set_cycle_properties(m_dot = 0.014, target='Pump:su')
 
    orc_cycle.set_cycle_properties(T=12 + 273.15, P = 5*1e5, fluid='Water', m_dot=0.2, target="Condenser:su_C")
    orc_cycle.set_cycle_properties(T=243 + 273.15, P = 5*1e5, fluid='INCOMP::T66', m_dot=0.4, target='Evaporator:su_H')
 
    # Define guesses for pressures
    P_ev_guess = 31.5*1e5
    P_cd_guess = 1e5

    orc_cycle.set_cycle_parameters(SH_ev = 5, SC_cd = 5)
    # SC_cd = 5
 
    guesses = {
        "Evaporator:su_C-P": P_ev_guess,
        "Condenser:ex_H-P": P_cd_guess,
        "Condenser:su_H-P": P_cd_guess
        }
 
    # Define the residuals to converge on
    residuals_var = [
        ("Evaporator:ex_C-h"),
        ("Condenser:ex_H-h")
    ]
 
    orc_cycle.set_cycle_guesses_residuals(guesses,residuals_var)
 
    # Solve the cycle with the defined guesses and residuals
    orc_cycle.solve()

