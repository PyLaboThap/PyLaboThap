from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI

class PumpCstEff(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W_cp = WorkConnector()
        # self.input = {}  # Store actual input values
        self.inputs_provided = False  # Track if inputs have been set manually

    def clear_intermediate_states(self):
        """Clear all intermediate states by resetting the connectors."""
        self.su.reset()
        self.ex.reset()

    def get_required_inputs(self):

        self.sync_inputs()
        # # If inputs have been provided manually, use those values
        # if self.inputs_provided:
        #     self.su.set_fluid(self.input_values.get('su_fluid'))
        #     if 'su_T' in self.input_values:
        #         self.su.set_T(self.input_values['su_T'])
        #     elif 'su_h' in self.input_values:
        #         self.su.set_h(self.input_values['su_h'])
        #     self.su.set_p(self.input_values.get('su_p'))
        #     self.ex.set_p(self.input_values.get('ex_p'))
        # else:
        #     # Populate input_values from connector properties
        #     if self.su.fluid is not None:
        #         self.input_values['su_fluid'] = self.su.fluid
        #     if self.su.T is not None:
        #         self.input_values['su_T'] = self.su.T
        #     elif self.su.h is not None:
        #         self.input_values['su_h'] = self.su.h
        #     if self.su.p is not None:
        #         self.input_values['su_p'] = self.su.p
        #     if self.ex.p is not None:
        #         self.input_values['ex_p'] = self.ex.p

        # Return a list of required inputs
        return ['su_p', 'su_T', 'ex_p', 'su_fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        elif self.su.h is not None:
            self.inputs['su_h'] = self.su.h
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.ex.p is not None:
            self.inputs['ex_p'] = self.ex.p

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        # Update the connectors based on the new inputs
        if 'su_fluid' in self.inputs:
            self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        elif 'su_h' in self.inputs:
            self.su.set_h(self.inputs['su_h'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'ex_p' in self.inputs:
            self.ex.set_p(self.inputs['ex_p'])

    def get_required_parameters(self):
        return ['eta_is']

    def print_setup(self):
        print("=== Pump Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")

        print("\nInputs:")
        for input_name in self.get_required_inputs():
            if input_name in self.input_values:
                print(f"  - {input_name}: {self.input_values[input_name]}")
            else:
                print(f"  - {input_name}: Not set")

        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")
        print("======================")

    def solve(self):
        # print(self.input_values)
        self.check_calculable()
        self.check_parametrized()
        print(self.calculable)
        if self.calculable and self.parametrized:
            try:
                print('Inlet pump:')
                print(self.su.print_resume())
                
                # Calculate the outlet enthalpy based on efficiency
                h_ex_is = PropsSI('H', 'P', self.ex.p, 'S', self.su.s, self.su.fluid)
                h_ex = self.su.h + (h_ex_is - self.su.h) * self.params['eta_is']
                self.ex.set_h(h_ex)
                self.ex.set_fluid(self.su.fluid)
                self.ex.set_m_dot(self.su.m_dot)
                
                self.defined = True
            except Exception as e:
                self.defined = False
                print(f"Convergence problem in pump model: {e}")

            print('Outlet pump:')
            print(self.ex.print_resume())

    def print_results(self):
        print("=== Pump Results ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")

        print("\nResults:")
        print(f"  - h_ex: {self.ex.h}")
        print(f"  - T_ex: {self.ex.T}")
        print("=========================")
