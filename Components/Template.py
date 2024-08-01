from Base_Component import BaseComponent
from Connectors.Mass_connector import Mass_connector
from Connectors.Work_connector import Work_connector
from Connectors.Heat_connector import Heat_connector


class Template(BaseComponent):
    def __init__(self):
        super().__init__()
        # Define connectors here
        self.su = Mass_connector()
        self.ex = Mass_connector() # Mass_connector
        self.work_su = Work_connector()
        self.heat_ex = Heat_connector()

    def get_required_inputs(self):
        
        if self.inputs == {}:
            "If the inputs are not set directly BUT throught the connectors"
            # Example where the supply temperature is needed:
            if self.su.T is not None:
                self.inputs['su_T'] = self.su.T
        
        if self.inputs != {}:
            "If the inputs are set directly, we put them in the connectors"
            # Example where the supply temperature is needed:
            if 'su_T' in self.inputs:
                self.su.set_T(self.inputs['su_T'])
        
        # Return the required inputs
        return ['su_T']

    def get_required_parameters(self):
        # Return the required parameters
        return [
            'param1', 'param2', 'param3'
        ]
    
    def print_setup(self):
        # Print the setup of the component (the connectors and the inputs)

        print("=== Expander Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}")
        print(f"  - ex: ...")
        print(f"  - W_dot: ...")
        print(f"  - Q_dot_amb: ...")

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

    def System(self):
        # Define the system of equations to solve
        # This function should return the residuals of the equations
        # The residuals should be stored in self.res
        pass

    def solve(self):
        # Check if the component is completely defined
        self.check_calculable()
        self.check_parametrized()

        if self.calculable and self.parametrized:
            # If the component is completely defined, solve it
            self.System()
        else:
            if not self.calculable:
                print("Input of the component not completely known. Required inputs:")
                for input in self.get_required_inputs():
                    if input not in self.inputs:
                        print(f"  - {input}")
            if not self.parametrized:
                print("Parameters of the component not completely known. Required parameters:")
                for param in self.get_required_parameters():
                    if param not in self.params:
                        print(f"  - {param}")
    
    def print_results(self):
        # Print the results of the component
        if self.defined:
            print("=== Expander Results ===")
            print(f"T_ex: {self.ex.T}")

        else:
            print("Expander component is not defined. Ensure it is solved first.")


#-------------------------------------------------------------
if __name__ == '__main__':
    # Example usage
    template = Template()

    #Print the set up
    template.print_setup()

    "If the inputs are not set directly BUT throught the connectors"
    # Set fluid states for connectors
    template.su.set_fluid('R134a')
    template.ex.set_fluid('R134a')

    # Set properties for su connector
    template.su.set_T(300)

    # Set properties for ex connector
    
    # Set properties for work connector

    # Set properties for heat connector

    "If the inputs are set directlys"
    # Setting inputs
    template.set_inputs(
        su_T=300,
        su_fluid='R134a'  # Make sure to include fluid information
    )
    
    # Setting parameters
    template.set_parameters(
        param1=1, param2=2, param3=3
    )

    # Solve the component
    template.solve()
    template.print_results()

