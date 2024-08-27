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


class RC(Cycle):
    def __init__(self, fluid=None):
        super().__init__(fluid)
        self.tolerance = 1e-6  # Convergence tolerance for residuals
        self.prev_pressures = {}  # Store pressures between iterations

    def solve(self):
        self.set_cycle_guesses()

        iteration = 0
        max_iterations = 100
        converged = False

        while not converged and iteration < max_iterations:
            print(f"Iteration {iteration}: Solving the cycle.")
            
            # Calculate initial residuals based on guesses or previous iteration
            initial_residuals = self.get_residuals()

            # Recursively solve the cycle starting from the pump
            visited = set()  # Reset visited components for each iteration
            self.recursive_solve(self.get_component("Pump"), visited)

            # Calculate final residuals after solving
            final_residuals = self.get_residuals()

            # Check if the residuals are within the tolerance
            converged = self.check_residuals(initial_residuals, final_residuals)

            # Store calculated pressures for the next iteration
            self.store_pressures()

            # Update guesses for the next iteration
            self.set_cycle_guesses()

            iteration += 1

        if converged:
            print("Cycle solved successfully.")
        else:
            print("Cycle did not converge within the maximum number of iterations.")

    def set_cycle_guesses(self):
        # Update pressures with previous iteration values if available, otherwise use initial guesses
        P_ev_guess = self.prev_pressures.get('P_ev', PropsSI('P', 'T', 90 + 273.15, 'Q', 0.5, self.fluid))
        P_cd_guess = self.prev_pressures.get('P_cd', PropsSI('P', 'T', 25 + 273.15, 'Q', 0.5, self.fluid))

        subcooling = 5
        superheating = 5

        T_su_pp = PropsSI('T', 'P', P_cd_guess, 'Q', 0, self.fluid) - subcooling
        self.set_cycle_properties(target="Pump:su", T=T_su_pp, P=P_cd_guess)

        self.set_cycle_properties(target="Pump:ex", P=P_ev_guess)
        T_ex_ev = PropsSI('T', 'P', P_ev_guess, 'Q', 1, self.fluid) + superheating
        self.set_cycle_properties(target="Evaporator:ex_wf", P=P_ev_guess, T=T_ex_ev)
        self.set_cycle_properties(target="Expander:ex", P=P_cd_guess)

    def store_pressures(self):
        # Store pressures from the evaporator and condenser models
        evaporator = self.get_component("Evaporator")
        condenser = self.get_component("Condenser")

        self.prev_pressures['P_ev'] = evaporator.model.ex_wf.p
        self.prev_pressures['P_cd'] = condenser.model.ex_wf.p


    def recursive_solve(self, component, visited):
        # Prevent infinite loops by skipping already visited components
        if component in visited:
            return

        # Mark the component as visited and solve it
        visited.add(component)
        print(f"Solving {component.name}")
        component.solve()

        # Recursively solve connected components
        for next_component in component.next.values():
            self.recursive_solve(next_component, visited)

    def get_residuals(self):
        # Retrieve components that define the residuals
        evaporator = self.get_component("Evaporator")
        condenser = self.get_component("Condenser")
        pump = self.get_component("Pump")


        # Calculate residuals based on enthalpies
        h_ex_ev = evaporator.model.ex_wf.h
        h_ex_cd = condenser.model.ex_wf.h

        return h_ex_ev, h_ex_cd

    def check_residuals(self, residuals_i, residuals_f):
        # Calculate the difference between initial and final residuals
        residual_h_ex_ev = abs(residuals_i[0] - residuals_f[0])
        residual_h_ex_cd = abs(residuals_i[1] - residuals_f[1])

        # Output residuals for debugging
        print(f"Residual Evaporator: {residual_h_ex_ev}, Residual Condenser: {residual_h_ex_cd}")

        # Check if residuals are within the tolerance
        return residual_h_ex_ev < self.tolerance and residual_h_ex_cd < self.tolerance


# Example usage of the Rankine Cycle (RC)
if __name__ == "__main__":
    orc_cycle = RC(fluid='R245fa')

    # Add components
    PUMP = PumpCstEff()
    EVAP = HXPinchCst()
    EXP = ExpanderCstEff()
    COND = HXPinchCst()

    # Set component parameters
    PUMP.set_parameters(eta_is=0.6)
    EVAP.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='evaporator')
    COND.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='condenser')
    EXP.set_parameters(eta_is=0.8)

    # Add components to the cycle
    orc_cycle.add_component(PUMP, "Pump")
    orc_cycle.add_component(EVAP, "Evaporator")
    orc_cycle.add_component(EXP, "Expander")
    orc_cycle.add_component(COND, "Condenser")

    # Link components
    orc_cycle.link_components("Pump", "m-ex", "Evaporator", "m-su_wf")
    orc_cycle.link_components("Evaporator", "m-ex_wf", "Expander", "m-su")
    orc_cycle.link_components("Expander", "m-ex", "Condenser", "m-su_wf")
    orc_cycle.link_components("Condenser", "m-ex_wf", "Pump", "m-su")

    # Set the cycle properties
    orc_cycle.set_cycle_properties(m_dot=0.06, target='Pump:su')

    orc_cycle.set_cycle_properties(T=30 + 273.15, fluid='Water', m_dot=0.4, target='Condenser:su_sf')
    orc_cycle.set_cycle_properties(cp=4186, target='Condenser:su_sf')
    orc_cycle.set_cycle_properties(fluid='Water', target='Condenser:ex_sf')

    orc_cycle.set_cycle_properties(T=150 + 273.15, fluid='Water', m_dot=0.4, target='Evaporator:su_sf')
    orc_cycle.set_cycle_properties(cp=4186, target='Evaporator:su_sf')
    orc_cycle.set_cycle_properties(fluid='Water', target='Evaporator:ex_sf')

    # Initial guesses for pressures
    T_ev_guess = 120 + 273.15
    P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R245fa')

    T_cd_guess = 3 + 273.15
    P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, 'R245fa')

    # Define pressure guesses
    pressure_guesses = {
        "P_ev": P_ev_guess,
        "P_cd": P_cd_guess
    }

    # Solve the cycle
    orc_cycle.solve()

    # Optional: Print a summary of the cycle
    # orc_cycle.print_summary()

