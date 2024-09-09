
"""
09-09-2024
Elise Neven
elise.neven@uliege.be
"""


import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class BaseComponent:
    def __init__(self):
        self.calculable = False
        self.parametrized = False
        self.solved = False
        self.inputs = {}
        self.params = {}
        self.guesses = {}

    def set_inputs(self, **inputs):
        for key, value in inputs.items():
            self.inputs[key] = value

    def set_parameters(self, **parameters):
        for key, value in parameters.items():
            self.params[key] = value
        self.check_parametrized()

    def set_guesses(self, **guesses):
        for key, value in guesses.items():
            self.guesses[key] = value

    def check_calculable(self):
        required_inputs = self.get_required_inputs() 
        self.calculable = all(self.inputs.get(inp) is not None for inp in required_inputs) # check if all required inputs are set
        return self.calculable

    def check_parametrized(self):
        required_params = self.get_required_parameters()
        self.parametrized = all(self.params.get(param) is not None for param in required_params) # check if all required parameters are set
        return self.parametrized

    def get_required_inputs(self):
        # This method should be overridden in derived classes
        return []

    def get_required_parameters(self):
        # This method should be overridden in derived classes
        return []
    
    def get_required_guesses(self):

        return []

    def solve(self):
        # This method should be overridden in derived classes
        raise NotImplementedError("The 'solve' method should be implemented in derived classes.")
    

    def plot_component(self, inputs_names=None, outputs_names=None, parameters_names=None):
        """
        Plot a visual representation of the component with inputs/outputs and parameters.
        """

        # Enable LaTeX rendering in Matplotlib
        plt.rcParams['text.usetex'] = True

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(8, 8))

        # Draw the block (component)
        block = patches.FancyBboxPatch((0.4, 0.6), 0.3, 0.1, boxstyle="round,pad=0.1", 
                                        edgecolor="black", facecolor="#D4E9C7", zorder=2)
        ax.add_patch(block)
        ax.text(0.55, 0.65, self.__class__.__name__, horizontalalignment='center', 
                verticalalignment='center', fontsize=16, fontweight='bold')

        # Define positions for inputs, outputs, and parameters
        input_pos = [0.2, 0.65]  # Position for input
        output_pos = [0.8, 0.65]  # Position for output
        param_pos = [0.55, 0.4]  # Position for parameters

        # Inputs
        if inputs_names:
            ax.text(input_pos[0] - 0.15, input_pos[1] + 0.2, r'\textbf{Inputs}', fontsize=15, fontweight='bold', color='black')
            ax.annotate('', xy=(input_pos[0] + 0.1, input_pos[1]), xytext=(input_pos[0] - 0.05, input_pos[1]),
                        fontsize=12, color='black',
                        arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

             # Calculate spacing for inputs
            num_inputs = len(inputs_names)
            vertical_spacing = 0.1
            start_y = input_pos[1] - (num_inputs - 1) * vertical_spacing / 2

            # List all inputs next to the arrow
            for i, input_name in enumerate(inputs_names):
                y = start_y + i * vertical_spacing
                ax.text(input_pos[0] - 0.15, y, f'${input_name}$', fontsize=15, color='black', verticalalignment='center')


        # Outputs
        if outputs_names:
            ax.text(output_pos[0] + 0.2, output_pos[1] + 0.2, r'\textbf{Outputs}', fontsize=15, fontweight='bold', color='black')
            ax.annotate('', xy=(output_pos[0] + 0.15, output_pos[1]), xytext=(output_pos[0], output_pos[1]),
                        fontsize=12, color='black',
                        arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

            # Calculate spacing for outputs
            num_outputs = len(outputs_names)
            vertical_spacing = 0.1
            start_y = output_pos[1] - (num_outputs - 1) * vertical_spacing / 2

            # List all outputs next to the arrow
            for i, output_name in enumerate(outputs_names):
                y = start_y + i * vertical_spacing
                ax.text(output_pos[0] + 0.2, y, f'${output_name}$', fontsize=15, color='black', verticalalignment='center')


        # Parameters
        if parameters_names:
            ax.text(param_pos[0], param_pos[1], r'\textbf{Parameters}', fontsize=15, fontweight='bold', color='black', horizontalalignment='center')
            ax.annotate('', xy=(param_pos[0], param_pos[1] + 0.1), xytext=(param_pos[0], param_pos[1] + 0.05),
                        fontsize=12, color='black',
                        arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

            # Calculate spacing for parameters
            num_params = len(parameters_names)
            vertical_spacing = 0.1
            start_y = param_pos[1] - (num_params - 1) * vertical_spacing / 2

            # List all parameters below "Parameters"
            for i, param_name in enumerate(parameters_names):
                y = start_y - i * vertical_spacing
                ax.text(param_pos[0], y-0.05, f'${param_name}$', fontsize=15, color='black', verticalalignment='center', horizontalalignment='center')

        # Plot formatting
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')  # Hide axis


    def plot_connectors(self, supply_connectors_names=None, exhaust_connectors_names=None):

        """
        Plot a visual representation of the component with inputs/outputs and parameters.
        """

        # Enable LaTeX rendering in Matplotlib
        plt.rcParams['text.usetex'] = True

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(8, 8))

        # Draw the block (component)
        block = patches.FancyBboxPatch((0.4, 0.6), 0.3, 0.1, boxstyle="round,pad=0.1", 
                                        edgecolor="black", facecolor="skyblue", zorder=2)
        ax.add_patch(block)
        ax.text(0.55, 0.65, self.__class__.__name__, horizontalalignment='center', 
                verticalalignment='center', fontsize=16, fontweight='bold')

        # Define positions for inputs, outputs, and parameters
        input_pos = [0.2, 0.65]  # Position for input
        output_pos = [0.8, 0.65]  # Position for output

        # Inputs
        if supply_connectors_names:
            ax.annotate('', xy=(input_pos[0] + 0.1, input_pos[1]), xytext=(input_pos[0] - 0.05, input_pos[1]),
                        fontsize=12, color='black',
                        arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

             # Calculate spacing for inputs
            num_inputs = len(supply_connectors_names)
            vertical_spacing = 0.1
            start_y = input_pos[1] - (num_inputs - 1) * vertical_spacing / 2

            # List all inputs next to the arrow
            for i, input_name in enumerate(supply_connectors_names):
                y = start_y + i * vertical_spacing
                ax.text(input_pos[0] - 0.15, y, f'${input_name}$', fontsize=15, color='black', verticalalignment='center')


        # Outputs
        if exhaust_connectors_names:
            ax.annotate('', xy=(output_pos[0] + 0.15, output_pos[1]), xytext=(output_pos[0], output_pos[1]),
                        fontsize=12, color='black',
                        arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

            # Calculate spacing for outputs
            num_outputs = len(exhaust_connectors_names)
            vertical_spacing = 0.1
            start_y = output_pos[1] - (num_outputs - 1) * vertical_spacing / 2

            # List all outputs next to the arrow
            for i, output_name in enumerate(exhaust_connectors_names):
                y = start_y + i * vertical_spacing
                ax.text(output_pos[0] + 0.2, y, f'${output_name}$', fontsize=15, color='black', verticalalignment='center')


        # Plot formatting
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')  # Hide axis




