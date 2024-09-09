from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

class CompressorCstEff(BaseComponent):
    """
    Component: Compressor
    Model: Constant isentropic efficiency
    -----------------------------------------------------------
    Connectors:
        su (MassConnector): Mass connector for the suction side.
        ex (MassConnector): Mass connector for the exhaust side.
        W_cp (WorkConnector): Work connector.
    -----------------------------------------------------------
    Parameters:
        eta_is: Isentropic efficiency.
    -----------------------------------------------------------
    Inputs:
        su_p: Suction side pressure.
        su_T: Suction side temperature.
        ex_p: Exhaust side pressure.
        su_fluid: Suction side fluid.
    -----------------------------------------------------------
    Ouputs:
        h_ex: Exhaust side specific enthalpy.
        T_ex: Exhaust side temperature.
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector() # Mass_connector
        self.W_cp = WorkConnector()

    def get_required_inputs(self):
        if self.inputs == {}:
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

        return ['su_p', 'su_T', 'ex_p', 'su_fluid']
    
    def get_required_parameters(self):
        return [
            'eta_is',
        ]
    
    def print_setup(self):
        print("=== Compressor Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")


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

    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if self.calculable and self.parametrized:

            h_ex_is = PropsSI('H', 'P', self.ex.p, 'S', self.su.s, self.su.fluid)
            h_ex = self.su.h + (h_ex_is - self.su.h) / self.params['eta_is']
            self.ex.set_h(h_ex)

    def print_results(self):
        print("=== Expander Results ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")

        print("\nResults:")
        print(f"  - h_ex: {self.ex.h}")
        print(f"  - T_ex: {self.ex.T}")
        print("=========================")


    def plot_component_comp_cst_eff(self):

        inputs_names = ['P_{su}', 'T_{su}', 'P_{ex}']
        outputs_names = ['h_{ex}', 'T_{ex}']
        parameters_names = ['\eta_{is}']

        self.plot_component(inputs_names, outputs_names, parameters_names)
        
        # Save the plot as a PNG file
        output_dir = '../docs/figures/component'  # Path to folder one level up
        os.makedirs(output_dir, exist_ok=True)  # Create the folder if it doesn't exist
        file_path = os.path.join(output_dir, 'constant_isentropic_efficiency_compressor_in_out.png')
        plt.savefig(file_path, format='png', bbox_inches='tight')  # Save the figure

        plt.show()

    def plot_connectors_comp_cst_eff(self):

        supply_connectors_names = ['\dot{m}_{su} \cdot h_{su}', '\dot{W}_{cp}']
        exhaust_connectors_names = ['\dot{m}_{ex} \cdot h_{ex}']

        self.plot_connectors(supply_connectors_names, exhaust_connectors_names)
        
        # Save the plot as a PNG file
        output_dir = '../docs/figures/component'  # Path to folder one level up
        os.makedirs(output_dir, exist_ok=True)  # Create the folder if it doesn't exist
        file_path = os.path.join(output_dir, 'constant_isentropic_efficiency_compressor_connectors.png')
        plt.savefig(file_path, format='png', bbox_inches='tight')  # Save the figure

        plt.show()
