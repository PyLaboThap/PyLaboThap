# -*- coding: utf-8 -*-
"""
Created on Wed Jul 07 11:47:52 2024
    
@author: elise neven
"""

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.heat_exchanger.pinch_cst.simulation_model import HXPinchCst
from component.volumetric_machine.expander.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from component.pump.constant_efficiency.simulation_model import PumpCstEff

from CoolProp.CoolProp import PropsSI

# Second option:
class Cycle:
    class Component:
        def __init__(self, name, model, fluid=None):
            self.name = name
            self.model = model
            self.previous = {}
            self.next = {}
            self.fluid = fluid

        def add_previous(self, port, component):
            self.previous[port] = component

        def add_next(self, port, component):
            self.next[port] = component

        def link(self, output_port, target_component, input_port):
            connector_type = output_port.split('-')[0]

            if connector_type == "m":  # Mass connector
                connector = MassConnector(fluid=self.fluid)
            elif connector_type == "q":  # Heat connector
                connector = HeatConnector()
            elif connector_type == "w":  # Work connector
                connector = WorkConnector()
            else:
                raise ValueError(f"Unknown connector type: {connector_type}")

            setattr(self.model, output_port.split('-')[1], connector)
            setattr(target_component.model, input_port.split('-')[1], connector) # Voir si ça fait juste référence ou si ça crée un nouvel objet

            self.add_next(output_port, target_component)
            target_component.add_previous(input_port, self)

            print(f"Linked {self.name}.{output_port} to {target_component.name}.{input_port}")

        def set_properties(self, port_name, **kwargs): #Je pense que ça vient d'cic!!!
            port = getattr(self.model, port_name)
            # print("Update values for", port_name, ":", kwargs) # Non en fait on dirait c'est okay???
            port.set_properties(**kwargs)

        def solve(self):
            self.model.check_calculable()
            if self.model.calculable:
                self.model.solve()

        def clear_intermediate_states(self):
            print(f"Clearing intermediate states for {self.name}")
            self.model.clear_intermediate_states()

    def __init__(self, fluid=None):
        self.components = {}  # Store components using a dictionary for direct access
        self.fixed_properties = {}
        self.guesses = {}
        self.fluid = fluid

    def add_component(self, model, name):
        component = Cycle.Component(name, model, self.fluid)
        self.components[name] = component

    def get_component(self, name):
        if name in self.components:
            return self.components[name]
        raise ValueError(f"Component {name} not found")

    def link_components(self, component1_name, output_port, component2_name, input_port):
        component1 = self.get_component(component1_name)
        component2 = self.get_component(component2_name)
        component1.link(output_port, component2, input_port)

    def set_cycle_properties(self, **kwargs):
        target = kwargs.pop('target')
        component_name, port_name = target.split(':')
        component = self.get_component(component_name)
        component.set_properties(port_name, **kwargs)

        if target not in self.fixed_properties:
            self.fixed_properties[target] = {}
        self.fixed_properties[target].update(kwargs)

    def set_guess(self, **kwargs):
        target = kwargs.pop('target')
        component_name, port_name = target.split(':')
        component = self.get_component(component_name)
        component.set_properties(port_name, **kwargs)

        if target not in self.guesses:
            self.guesses[target] = {}
        self.guesses[target].update(kwargs)

    


    def solve(self, variables_to_converge, tolerance=1e-5, max_iterations=100):
        print('variable to converge', variables_to_converge)
        for i in range(max_iterations):
            self.solve_loop()

            converged = True
            for target, variable_name in variables_to_converge:
                component_name, port_name = target.split(':')
                component = self.get_component(component_name)
                port = getattr(component.model, port_name)

                value = getattr(port, variable_name)
                guess = self.guesses[target][variable_name]

                if abs(value - guess) > tolerance:
                    converged = False
                    self.guesses[target][variable_name] = value
                    port.set_properties(**{variable_name: value})

            if converged:
                print(f"Converged after {i + 1} iterations")
                break
        else:
            print(f"Warning: Convergence not achieved after {max_iterations} iterations")

    def solve_loop(self): # will have to be changed
        for component in self.components.values():
            component.solve()
