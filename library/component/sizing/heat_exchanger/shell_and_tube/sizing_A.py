
import __init__

from library.connector.mass_connector import MassConnector
from library.component.base_component import BaseComponent

class SizingAHTX(BaseComponent):

    def __init__(self):

        super().__init__()

        self.su_H = MassConnector()
        self.su_C = MassConnector()

        self.ex_H = MassConnector()
        self.ex_C = MassConnector()

        self.fouling = None
        self.HTX_Type = None
        self.geom_input = None

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return['fouling', 'HTX_Type', 'geom_input']
    
    def sync_inputs(self):
        
        if self.fouling is not None:
            self.inputs['fouling'] = self.fouling
        if self.HTX_Type is not None:
            self.inputs['HTX_Type'] = self.HTX_Type
        if self.fouling is not None:
            self.inputs['geom_input'] = self.geom_input

        return

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'fouling' in self.inputs:
            self.fouling = self.inputs['fouling']
        if 'HTX_Type' in self.inputs:
            self.HTX_Type = self.inputs['HTX_Type']
        if 'geom_input' in self.inputs:
            self.geom_input = self.inputs['geom_input']

        return['su_C_fluid', 'su_C_h', 'su_C_m_dot', 'su_H_fluid', 'su_H_T', 'su_H_cp', 'su_H_m_dot']
    
    def get_required_parameters(self):
        return [
            'su_C_fluid', 'su_C_T', 'su_C_m_dot', 'su_C_p', 'su_H_fluid', 'su_H_T', 'su_H_p', 'su_H_m_dot'
        ]
    
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")

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

    def system(self):
        
 

        return 
    
    def opt_size():
        
        return

#%%

HTX = SizingAHTX()

HTX.set_parameters(
                    su_C_fluid = 'Water',
                    su_C_T = 273.15 + 24, # K
                    su_C_p = 101.3*1e3, # Pa
                    su_C_m_dot = 1000, # kg/s

                    su_H_fluid = 'Cyclopentane',
                    su_H_T = 273.15 + 38.43, # K
                    su_H_p = 51.75*1e3, # Pa
                    su_H_m_dot = 32.85 # kg/s                    
                    )



HTX.set_inputs(
                fouling = 0, HTX_Type = 'Shell&Tube', geom_input = ['Tube_OD','L_shell','Central_Spacing','Shell_ID']
                )

HTX.opt_size()
