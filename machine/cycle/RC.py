from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.heat_exchanger.pinch_cst.simulation_model import HXPinchCst
from component.volumetric_machine.expander.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from component.pump.constant_efficiency.simulation_model import PumpCstEff
from machine.cycle.circuit import Circuit


from CoolProp.CoolProp import PropsSI
import numpy as np



class RC(Circuit):
    def __init__(self, fluid=None):
        super().__init__(fluid)
        self.tolerance = 1e-6  # Convergence tolerance for residuals
        self.prev_pressures = {}  # Store pressures between iterations
        self.prev_residuals = None  # Store residuals from the previous iteration
        self.residuals_var = []  # Store the residuals to check

    def solve(self, start_key="Pump"):
        self.set_cycle_guesses()

        iteration = 0
        max_iterations = 10
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

        self.prev_pressures['P_ev'] = evaporator.model.ex_wf.p
        self.prev_pressures['P_cd'] = condenser.model.ex_wf.p

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
        self.set_cycle_properties(target="Evaporator:ex_wf", P=P_ev_guess, T=T_ex_ev)
        self.set_cycle_properties(target="Expander:ex", P=P_cd_guess)

    def recursive_solve(self, component, visited):
        # Prevent infinite loops by skipping already visited components
        if component in visited:
            return

        # Mark the component as visited and solve it
        visited.add(component)
        print(f"Solving {component.name}")
        if isinstance(component, Circuit.Component):
            component.solve()

        # Recursively solve connected components
        for next_component in component.next.values():
            self.recursive_solve(next_component, visited)


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

    # Set parameters for the cycle
    SC_cd = 5
    SH_ev = 5
    orc_cycle.set_cycle_parameters(SC_cd=SC_cd, SH_ev=SH_ev)

    # Initial guesses for pressures
    T_ev_guess = 120 + 273.15
    P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R245fa')

    T_cd_guess = 3 + 273.15
    P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, 'R245fa')

    # Define guesses and residuals
    guesses = {
        "Evaporator:su_wf-P": P_ev_guess,
        "Condenser:ex_wf-P": P_cd_guess,
    }

    residuals_var = [
        "Evaporator:ex_wf-h",
        "Condenser:ex_wf-h"
    ]

    orc_cycle.set_cycle_guesses_residuals(guesses, residuals_var)

    # Solve the cycle
    orc_cycle.solve()

