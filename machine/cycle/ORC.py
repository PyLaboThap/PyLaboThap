from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.heat_exchanger.pinch_cst.simulation_model import HXPinchCst
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
#         super().__init__(fluid)  # Ensure Cycle's initialization is done

#     def solve(self, pressure_guesses, residuals, tolerance=1e-3, max_iterations=1000):
#         def objective(x):
#             P_ev, rp = x
#             P_cd = rp * P_ev
#             print(x, 'x')
#             # Step 1: Set the initial guesses and calculate related properties
#             self.set_cycle_guesses(P_ev, P_cd)

#             # Access the evaporator and condenser components directly
#             evaporator = self.get_component("Evaporator")
#             condenser = self.get_component("Condenser")
#             pump = self.get_component("Pump")
#             expander = self.get_component("Expander")

#             # Step 2: Evaluate initial enthalpies from the component ports
#             print(evaporator)
#             # print(pump.model.su.p)
#             # print(pump.model.su.h)
#             print(condenser.model.ex_wf.p)
#             print(evaporator.model.ex_wf.h)
#             h_ex_ev_init = evaporator.model.ex_wf.h
#             h_ex_cd_init = condenser.model.ex_wf.h

#             print(self.components.values())
#             # Step 3: Solve the cycle with the guesses
#             for component in self.components.values():
#                 component.solve()
#                 print(component)

#             print(condenser.model.ex_wf.p)
#             print(pump.model.su.p)
#             print(pump.model.ex.p)
#             # Step 4: Recalculate the enthalpies after solving
#             h_ex_ev_final = evaporator.model.ex_wf.h
#             h_ex_cd_final = condenser.model.ex_wf.h

#             # Step 5: Compute the residuals as differences in enthalpies
#             residual_h_ex_ev = h_ex_ev_init - h_ex_ev_final
#             residual_h_ex_cd = h_ex_cd_init - h_ex_cd_final
#             print(residual_h_ex_ev, 'residual_h_ex_ev')
#             print(residual_h_ex_cd, 'residual_h_ex_cd')
#             res = np.sqrt(residual_h_ex_ev ** 2 + residual_h_ex_cd ** 2)
#             print(res, 'res')
#             # res = [residual_h_ex_ev, residual_h_ex_cd]
#             return res

#         # Set up initial guesses for P_ev and P_cd
#         x0 = [pressure_guesses["P_ev"], 3]
#         bounds = [(1e4, 5e6), (1, 10)]  # Example: (10 kPa to 5 MPa)

#         # Optimize using scipy's minimize function
#         result = minimize(objective, x0, bounds=bounds, tol=tolerance, options={'maxiter': max_iterations})
#         # result = fsolve(objective, x0)
#         print(result, 'result')

#         if result.success:
#             optimized_pressures = result.x
#             print(f"Optimization succeeded: P_ev = {optimized_pressures[0]}, P_cd = {optimized_pressures[1]}")
#         else:
#             print("Optimization failed.")

#     def set_cycle_guesses(self, P_ev_guess, P_cd_guess):
#         """
#         Set the initial guesses for the evaporator and condenser pressures and calculate related properties.
#         """
#         subcooling = 5  # Adjust as needed
#         superheating = 5  # Adjust as needed

#         # Step 1: Calculate temperature at the pump inlet based on P_cd_guess and subcooling
#         T_su_pp = PropsSI('T', 'P', P_cd_guess, 'Q', 0, self.fluid) - subcooling
#         self.set_guess(target="Pump:su", T=T_su_pp, P=P_cd_guess, fluid=self.fluid)
#         # self.set_guess(target="Pump:su", P=P_cd_guess)

#         # Step 2: Set guessed pressure at the expander inlet (evaporation pressure)
#         self.set_guess(target="Pump:ex", P=P_ev_guess)
#         T_ex_ev = PropsSI('T', 'P', P_ev_guess, 'Q', 1, self.fluid) + superheating
#         self.set_guess(target="Evaporator:ex_wf", P=P_ev_guess, T=T_ex_ev)
#         self.set_guess(target="Expander:ex", P=P_cd_guess)

class RC(Cycle):
    def __init__(self, fluid=None):
        super().__init__(fluid)
        self.tolerance = 1e-6  # Convergence tolerance for residuals
        self.prev_pressures = {}  # To store pressures between iterations

    def solve(self):
        # Set initial guesses for pressures and temperatures
        self.set_cycle_guesses()

        iteration = 0
        max_iterations = 100
        converged = False

        while not converged and iteration < max_iterations:
            print(f"Iteration {iteration}: Solving the cycle.")
            
            # En ce moment mon problème, c'est que à la deuxième itération, les guess reviennent à leur valeur initiale et on refait une boucle avec les mêmes valeurs
            # Calculate initial residuals (based on guesses or previous iteration)
            residuals_i = self.get_residuals()
            pump = self.get_component("Pump")
            print(pump.model.su.p, 'pump.model.su.p')
            print(pump.model.ex.p, 'pump.model.ex.p')
            # Solve the cycle starting from the Pump recursively
            visited = set()  # Reset visited components for each iteration
            self.recursive_solve(self.get_component("Pump"), visited)

            # Calculate final residuals after solving
            residuals_f = self.get_residuals()

            # Check if the residuals are within the tolerance
            converged = self.check_residuals(residuals_i, residuals_f)
            # Store the calculated pressures to use as inputs for the next iteration
            self.store_pressures()

            # Clear intermediate states before the next iteration
            self.clear_intermediate_states()

            self.set_cycle_properties(m_dot = 0.06, target='Pump:su')

            self.set_cycle_properties(T=30 + 273.15, fluid='Water', m_dot=0.4, target='Condenser:su_sf')
            self.set_cycle_properties(cp=4186, target='Condenser:su_sf')
            self.set_cycle_properties(fluid='Water', target='Condenser:ex_sf')

            self.set_cycle_properties(T=150 + 273.15, fluid='Water', m_dot=0.4, target='Evaporator:su_sf')
            self.set_cycle_properties(cp=4186, target='Evaporator:su_sf')
            self.set_cycle_properties(fluid='Water', target='Evaporator:ex_sf')

            # Re-set cycle guesses with the updated pressures
            self.set_cycle_guesses()

            iteration += 1

        if converged:
            print("Cycle solved successfully.")
        else:
            print("Cycle did not converge within the maximum number of iterations.")

    def set_cycle_guesses(self):
        # Use the updated pressures if available, otherwise use the initial guess
        if 'P_ev' in self.prev_pressures:
            P_ev_guess = self.prev_pressures['P_ev']
        else:
            T_ev_guess = 90 + 273.15
            P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, self.fluid)

        if 'P_cd' in self.prev_pressures:
            P_cd_guess = self.prev_pressures['P_cd']
        else:
            T_cd_guess = 25 + 273.15
            P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, self.fluid)

        subcooling = 5
        superheating = 5
        print(P_cd_guess, 'P_cd_guess')
        print(P_ev_guess, 'P_ev_guess')
        T_su_pp = PropsSI('T', 'P', P_cd_guess, 'Q', 0, self.fluid) - subcooling
        self.set_cycle_properties(target="Pump:su", T=T_su_pp, P=P_cd_guess)

        self.set_cycle_properties(target="Pump:ex", P=P_ev_guess)
        T_ex_ev = PropsSI('T', 'P', P_ev_guess, 'Q', 1, self.fluid) + superheating
        self.set_cycle_properties(target="Evaporator:ex_wf", P=P_ev_guess, T=T_ex_ev)
        self.set_cycle_properties(target="Expander:ex", P=P_cd_guess)

    def store_pressures(self):
        # Store the pressures calculated from the evaporator and condenser
        evaporator = self.get_component("Evaporator")
        condenser = self.get_component("Condenser")

        self.prev_pressures['P_ev'] = evaporator.model.ex_wf.p
        self.prev_pressures['P_cd'] = condenser.model.ex_wf.p

    def clear_intermediate_states(self):

        for component_name, component in self.components.items():
            print(f"Clearing states for component: {component_name}")
            component.clear_intermediate_states()

    def recursive_solve(self, component, visited):
        # If the component has already been visited, return to avoid infinite loop
        if component in visited:
            return

        # Mark the component as visited
        visited.add(component)

        # Check if the component is calculable and hasn't been solved yet
        print(f"Solving {component.name}")
        component.solve()

        # Recursively solve the next connected components
        for next_component in component.next.values():
            self.recursive_solve(next_component, visited)

    def get_residuals(self):
        # Get the components that define the residuals (Evaporator and Condenser)
        evaporator = self.get_component("Evaporator")
        condenser = self.get_component("Condenser")
        pump = self.get_component("Pump")
        print(evaporator.model.ex_wf.p, 'evaporator.model.ex_wf.p')
        print(condenser.model.ex_wf.p, 'condenser.model.ex_wf.p')
        print(pump.model.ex.p, 'pump.model.ex.p')
        print(pump.model.su.p, 'pump.model.su.p')
        

        # Calculate the residuals
        h_ex_ev = evaporator.model.ex_wf.h #Maybe find a way to say those come from the guesses and the other ones from the model
        h_ex_cd = condenser.model.ex_wf.h
        print(h_ex_ev, 'h_ex_ev')
        print(h_ex_cd, 'h_ex_cd')
        return h_ex_ev, h_ex_cd

    def check_residuals(self, residuals_i, residuals_f):
        # For the first iteration, store the initial enthalpies as the reference
        # if not hasattr(self, 'h_ex_ev_initial'):
        #     self.h_ex_ev_initial = h_ex_ev
        #     self.h_ex_cd_initial = h_ex_cd
        # print(residuals_i, residuals_f, 'residuals')
        residual_h_ex_ev = abs(residuals_i[0] - residuals_f[0])
        residual_h_ex_cd = abs(residuals_i[1] - residuals_f[1])

        print(F"Residual Evaporator: {residual_h_ex_ev}, Residual Condenser: {residual_h_ex_cd}")
        # print(f"Residual Evaporator: {residual_h_ex_ev}, Residual Condenser: {residual_h_ex_cd}")
        # Check if the residuals are within the tolerance
        # print(residual_h_ex_ev < self.tolerance)
        # print(residual_h_ex_cd < self.tolerance)
        # print(self.tolerance, 'tolerance')
        return residual_h_ex_ev < self.tolerance and residual_h_ex_cd < self.tolerance


# class RC(Cycle):
#     def __init__(self, fluid=None):
#         super().__init__(fluid)
#         self.tolerance = 1e-6  # Convergence tolerance for residuals

#     def solve(self):
#         # Set initial guesses for pressures and temperatures
#         self.set_cycle_guesses() # Là on pourrait donner les guesses
#         # print(self.get_component("Pump").model.calculable)
#         # Start solving the cycle from the first component (Pump)
#         # self.recursive_solve(self.get_component("Pump"))
#         # print('should solve a first time then?')
#         visited = set()  # Use a set to track visited components
#         # Check the residuals (enthalpies at key points) and decide if another pass is needed
#         # converged = self.check_residuals()
#         iteration = 0
#         max_iterations = 100
#         converged = False
#         while not converged and iteration < max_iterations:
#             # If not converged, update guesses and solve again
#             print(f"Iteration {iteration}: Not converged, updating guesses and solving again.")
#             # print('P_su_pp = ', self.get_component("Pump").model.su.p)
#             # print('P_ex_pp = ', self.get_component("Pump").model.ex.p)
#             residuals_i = self.get_residuals()
#             visited = set()  # Reset visited components for the new iteration
#             self.recursive_solve(self.get_component("Pump"), visited) # ici get_component il faut mettre strat key
#             residuals_f = self.get_residuals()
#             # print('P_su_pp = ', self.get_component("Pump").model.su.p)
#             # print('P_ex_pp = ', self.get_component("Pump").model.ex.p)
#             # print(residuals_i, residuals_f)
#             converged = self.check_residuals(residuals_i, residuals_f)
#             iteration += 1

#         if converged:
#             print("Cycle solved successfully.")
#         else:
#             print("Cycle did not converge within the maximum number of iterations.")

#     def set_cycle_guesses(self):
#         # Define guesses for pressures
#         T_ev_guess = 120 + 273.15
#         P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R245fa')

#         T_cd_guess = 35 + 273.15
#         P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, 'R245fa')

#         subcooling = 5  # Adjust as needed
#         superheating = 5

#         T_su_pp = PropsSI('T', 'P', P_cd_guess, 'Q', 0, self.fluid) - subcooling
#         self.set_guess(target="Pump:su", T=T_su_pp, P=P_cd_guess)

#         self.set_guess(target="Pump:ex", P=P_ev_guess)
#         T_ex_ev = PropsSI('T', 'P', P_ev_guess, 'Q', 1, self.fluid) + superheating
#         self.set_guess(target="Evaporator:ex_wf", P=P_ev_guess, T=T_ex_ev)
#         self.set_guess(target="Expander:ex", P=P_cd_guess)

#     def recursive_solve(self, component, visited):
#         # If the component has already been visited, return to avoid infinite loop
#         if component in visited:
#             return

#         # Mark the component as visited
#         visited.add(component)

#         # Check if the component is calculable and hasn't been solved yet
#         print(f"Solving {component.name}")
#         component.solve()

#         # Recursively solve the next connected components
#         for next_component in component.next.values():
#             self.recursive_solve(next_component, visited)

#     def get_residuals(self):
#         # Get the components that define the residuals (Evaporator and Condenser)
#         evaporator = self.get_component("Evaporator")
#         condenser = self.get_component("Condenser")

#         # Calculate the residuals
#         h_ex_ev = evaporator.model.ex_wf.h #Maybe find a way to say those come from the guesses and the other ones from the model
#         h_ex_cd = condenser.model.ex_wf.h

#         return h_ex_ev, h_ex_cd

#     def check_residuals(self, residuals_i, residuals_f):
#         # For the first iteration, store the initial enthalpies as the reference
#         # if not hasattr(self, 'h_ex_ev_initial'):
#         #     self.h_ex_ev_initial = h_ex_ev
#         #     self.h_ex_cd_initial = h_ex_cd
#         # print(residuals_i, residuals_f, 'residuals')
#         residual_h_ex_ev = abs(residuals_i[0] - residuals_f[0])
#         residual_h_ex_cd = abs(residuals_i[1] - residuals_f[1])

#         print(F"Residual Evaporator: {residual_h_ex_ev}, Residual Condenser: {residual_h_ex_cd}")
#         # print(f"Residual Evaporator: {residual_h_ex_ev}, Residual Condenser: {residual_h_ex_cd}")
#         # Check if the residuals are within the tolerance
#         # print(residual_h_ex_ev < self.tolerance)
#         # print(residual_h_ex_cd < self.tolerance)
#         # print(self.tolerance, 'tolerance')
#         return residual_h_ex_ev < self.tolerance and residual_h_ex_cd < self.tolerance
#     # def check_residuals(self):


#         # For the first iteration, store the initial enthalpies as the reference
#         if not hasattr(self, 'h_ex_ev_initial'):
#             self.h_ex_ev_initial = h_ex_ev
#             self.h_ex_cd_initial = h_ex_cd

#         residual_h_ex_ev = abs(self.h_ex_ev_initial - h_ex_ev)
#         residual_h_ex_cd = abs(self.h_ex_cd_initial - h_ex_cd)

#         print(F"Residual Evaporator: {residual_h_ex_ev}, Residual Condenser: {residual_h_ex_cd}")
#         # print(f"Residual Evaporator: {residual_h_ex_ev}, Residual Condenser: {residual_h_ex_cd}")
#         # Check if the residuals are within the tolerance
#         return residual_h_ex_ev < self.tolerance and residual_h_ex_cd < self.tolerance


# Example usage of ORC
if __name__ == "__main__":
    orc_cycle = RC(fluid='R245fa')

    # Add components 
    PUMP = PumpCstEff()
    EVAP = HXPinchCst()
    EXP = ExpanderCstEff()
    COND = HXPinchCst()

    # Set parameters
    PUMP.set_parameters(eta_is=0.6)
    EVAP.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='evaporator')
    COND.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='condenser')
    EXP.set_parameters(eta_is=0.8)

    orc_cycle.add_component(PUMP, "Pump")
    orc_cycle.add_component(EVAP, "Evaporator")
    orc_cycle.add_component(EXP, "Expander")
    orc_cycle.add_component(COND, "Condenser")

    # Link components
    orc_cycle.link_components("Pump", "m-ex", "Evaporator", "m-su_wf")
    orc_cycle.link_components("Evaporator", "m-ex_wf", "Expander", "m-su")
    orc_cycle.link_components("Expander", "m-ex", "Condenser", "m-su_wf")
    orc_cycle.link_components("Condenser", "m-ex_wf", "Pump", "m-su")

    # Set the inputs using set_properties
    orc_cycle.set_cycle_properties(m_dot = 0.06, target='Pump:su')

    orc_cycle.set_cycle_properties(T=30 + 273.15, fluid='Water', m_dot=0.4, target='Condenser:su_sf')
    orc_cycle.set_cycle_properties(cp=4186, target='Condenser:su_sf')
    orc_cycle.set_cycle_properties(fluid='Water', target='Condenser:ex_sf')

    orc_cycle.set_cycle_properties(T=150 + 273.15, fluid='Water', m_dot=0.4, target='Evaporator:su_sf')
    orc_cycle.set_cycle_properties(cp=4186, target='Evaporator:su_sf')
    orc_cycle.set_cycle_properties(fluid='Water', target='Evaporator:ex_sf')

    # Define guesses for pressures
    T_ev_guess = 120 + 273.15
    P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R245fa')

    T_cd_guess = 3 + 273.15
    P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, 'R245fa')

    pressure_guesses = {
        "P_ev": P_ev_guess,
        "P_cd": P_cd_guess
    }

    # Define the residuals to converge on
    residuals = [
        ("Evaporator:ex_wf", "h"),
        ("Condenser:ex_wf", "h")
    ]

    # Solve the cycle with the defined guesses and residuals
    # orc_cycle.solve(pressure_guesses, residuals)
    orc_cycle.solve()
    # Print a summary
    # orc_cycle.print_summary()
