from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.heat_exchanger.pinch_cst.simulation_model import HXPinchCst
from component.volumetric_machine.expander.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from component.pump.constant_efficiency.simulation_model import PumpCstEff
from machine.cycle.cycle import Cycle

from CoolProp.CoolProp import PropsSI


class ORC(Cycle):
    def __init__(self, fluid=None):
        super().__init__(fluid)  # Initialize the base class (Cycle)
        self.additional_params = {}

    def add_component(self, model, name):
        # Optionally, add ORC-specific logic here
        super().add_component(model, name)  # Call the base class method

    def set_additional_params(self, **kwargs):
        """Set additional parameters specific to ORC."""
        self.additional_params.update(kwargs)

    def solve(self, pressure_guesses, residuals, tolerance=1e-5, max_iterations=100):
        """
        Override solve method to adjust connector values based on P_ev and P_cd guesses,
        then call the base Cycle's solve method.
        
        Args:
            pressure_guesses (dict): Dictionary of guessed pressures (e.g., {"P_ev": P_ev_guess, "P_cd": P_cd_guess}).
            residuals (list): List of tuples specifying the residuals to be minimized.
            tolerance (float): Convergence tolerance.
            max_iterations (int): Maximum number of iterations.
        """

        # Extract guessed pressures
        P_ev_guess = pressure_guesses.get("P_ev")
        P_cd_guess = pressure_guesses.get("P_cd")

        # Step 1: Calculate initial temperature at the pump inlet based on P_cd_guess and subcooling
        # Assume a certain subcooling value (for example 5 degrees)
        subcooling = 5  # Adjust as needed
        T_su_pp = PropsSI('T', 'P', P_cd_guess, 'Q', 0, self.fluid) - subcooling
        h_su_pp = PropsSI('H', 'P', P_cd_guess, 'T', T_su_pp, self.fluid)

        # Set properties at the pump inlet
        self.set_properties(target="Pump:su", T=T_su_pp, P=P_cd_guess, fluid=self.fluid) #Condenser working fluid outlet state FIRST GUESS

        # Step 2: Set guessed pressure values at relevant ports
        self.set_guess(target="Pump:ex", P=P_ev_guess)  # Expander inlet pressure (evaporation pressure)

        superheating = 5  # Adjust as needed
        T_ex_ev = PropsSI('T', 'P', P_ev_guess, 'Q', 1, self.fluid) + superheating
        self.set_properties(target='Evaporator:su_wf', P=P_ev_guess, T=T_ex_ev)  # Evaporator working fluid inlet state FIRST GUESS

        # Step 3: Call the base class's solve method, focusing on the specified residuals
        super().solve(residuals, tolerance=tolerance, max_iterations=max_iterations)

    def print_summary(self):
        """Print summary of the ORC cycle solution."""
        print("ORC Cycle Summary:")
        for component in self.components:
            model = component["model"]
            print(f"Component: {component['name']}, Model: {model}")

# Example usage of ORC
if __name__ == "__main__":
    orc_cycle = ORC(fluid='R245fa')

    # Add components 
    PUMP = PumpCstEff()
    EVAP = HXPinchCst()
    COND = HXPinchCst()
    EXP = ExpanderCstEff()

    # Set parameters
    PUMP.set_parameters(eta_is=0.6)
    EVAP.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='evaporator')
    COND.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='condenser')
    EXP.set_parameters(eta_is=0.8)

    orc_cycle.add_component(PUMP, "Pump")
    orc_cycle.add_component(EVAP, "Evaporator")
    orc_cycle.add_component(COND, "Condenser")
    orc_cycle.add_component(EXP, "Expander")

    # Link components
    orc_cycle.link_components("Pump", "m-ex", "Evaporator", "m-su_wf")
    orc_cycle.link_components("Evaporator", "m-ex_wf", "Expander", "m-su")
    orc_cycle.link_components("Expander", "m-ex", "Condenser", "m-su_wf")
    orc_cycle.link_components("Condenser", "m-ex_wf", "Pump", "m-su")

    # Set the inputs using set_properties
    orc_cycle.set_properties(T=30 + 273.15, fluid='Water', m_dot=0.4, target='Condenser:su_sf')
    orc_cycle.set_properties(cp=4186, target='Condenser:su_sf')
    orc_cycle.set_properties(fluid='Water', target='Condenser:ex_sf')

    # Define guesses for pressures
    T_ev_guess = 120 + 273.15
    P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R245fa')

    T_cd_guess = 30 + 273.15
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
    orc_cycle.solve(pressure_guesses, residuals)

    # Print a summary
    orc_cycle.print_summary()
