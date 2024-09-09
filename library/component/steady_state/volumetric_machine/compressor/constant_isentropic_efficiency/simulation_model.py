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


    def plot_component(self):
        """
        Plot a visual representation of the component with connectors as inputs/outputs and parameters.
        """
        # Create figure and axis
        fig, ax = plt.subplots(figsize=(8, 8))

        # Draw the block (component)
        block = patches.FancyBboxPatch((0.4, 0.5), 0.3, 0.3, boxstyle="round,pad=0.1", 
                                        edgecolor="black", facecolor="#D4E9C7", zorder=2)
        ax.add_patch(block)
        ax.text(0.55, 0.65, self.__class__.__name__, horizontalalignment='center', 
                verticalalignment='center', fontsize=14, fontweight='bold')

        # Define positions for inputs, outputs, and parameters
        input_pos = [0.2, 0.65]  # Position for input
        output_pos = [0.8, 0.65]  # Position for output
        param_pos = [0.55, 0.25]  # Position for parameters

        # Inputs
        inputs = self.get_required_inputs()
        if inputs:
            ax.text(input_pos[0] - 0.1, input_pos[1] + 0.2, 'Inputs', fontsize=12, fontweight='bold', color='black')
            ax.annotate('', xy=(input_pos[0] + 0.1, input_pos[1]), xytext=(input_pos[0] - 0.05, input_pos[1]),
                        fontsize=12, color='black',
                        arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

            # List all inputs next to the arrow
            for i, input_name in enumerate(inputs):
                y = input_pos[1] + (i - len(inputs) / 2) * 0.1
                ax.text(input_pos[0] - 0.15, y, input_name, fontsize=12, color='black', verticalalignment='center')

        # Outputs
        outputs = ['h_ex', 'T_ex']  # Modify this depending on your outputs
        if outputs:
            ax.text(output_pos[0] + 0.05, output_pos[1] + 0.2, 'Outputs', fontsize=12, fontweight='bold', color='black')
            ax.annotate('', xy=(output_pos[0] + 0.15, output_pos[1]), xytext=(output_pos[0], output_pos[1]),
                        fontsize=12, color='black',
                        arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

            # List all outputs next to the arrow
            for i, output_name in enumerate(outputs):
                y = output_pos[1] + (i - len(outputs) / 2) * 0.1
                ax.text(output_pos[0] + 0.15, y, output_name, fontsize=12, color='black', verticalalignment='center')

        # Parameters
        params = self.get_required_parameters()
        if params:
            ax.text(param_pos[0], param_pos[1], 'Parameters', fontsize=12, fontweight='bold', color='black', horizontalalignment='center')
            ax.annotate('', xy=(0.55, param_pos[1] + 0.1), xytext=(0.55, param_pos[1] + 0.05),
                        fontsize=12, color='black',
                        arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

            # List all parameters below "Parameters"
            for i, param_name in enumerate(params):
                y = param_pos[1] - 0.1 * (i + 1)
                ax.text(param_pos[0], y, param_name, fontsize=12, color='black', verticalalignment='center', horizontalalignment='center')

        # Plot formatting
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')  # Hide axis

        # Save the plot as a PNG file
        output_dir = '../docs/figures/component'  # Update this to your desired folder path
        os.makedirs(output_dir, exist_ok=True)  # Create the folder if it doesn't exist
        file_path = os.path.join(output_dir, 'constant_isentropic_efficiency_compressor_in_out.png')
        plt.savefig(file_path, format='png', bbox_inches='tight')  # Save the figure

        plt.show()
