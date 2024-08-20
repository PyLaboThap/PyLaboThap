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



# Second option:
class Cycle:
    class Component:
        def __init__(self, name, model):
            self.name = name
            self.model = model
            self.previous = {}  # Dictionary to store connections to previous components
            self.next = {}  # Dictionary to store connections to next components

        def add_previous(self, port, component):
            self.previous[port] = component

        def add_next(self, port, component):
            self.next[port] = component

    def __init__(self):
        self.component_list = []  # List to store the components in the cycle
        self.connectors = {}
        self.fluid = None

    def add_component(self, model, name):
        component = self.Component(name, model)
        self.component_list.append(component)

    def get_component(self, name):
        # Retrieve a component by name
        for component in self.component_list:
            if component.name == name:
                return component
        raise ValueError(f"Component {name} not found")

    def link_components(self, component1_name, output_port, component2_name, input_port):
        component1 = self.get_component(component1_name)
        component2 = self.get_component(component2_name)

        # Determine connector type based on the port name prefix
        connector_type = output_port.split('-')[0]

        if connector_type == "m":  # Mass connector
            connector = MassConnector()
        elif connector_type == "q":  # Heat connector
            connector = HeatConnector()
        elif connector_type == "w":  # Work connector
            connector = WorkConnector()
        else:
            raise ValueError(f"Unknown connector type: {connector_type}")

        # Dynamically link the connectors
        setattr(component1.model, output_port.split('-')[1], connector)
        setattr(component2.model, input_port.split('-')[1], connector)

        # Update the component connections
        component1.add_next(output_port, component2)
        component2.add_previous(input_port, component1)

        print(f"Linked {component1_name}.{output_port} to {component2_name}.{input_port}")

    def set_properties(self, **kwargs):
        """
        Set properties for a specific component's port.
        Example usage: ORC.set_properties(T=117+273.15, fluid='R245fa', target='Condenser:su_wf')
        """
        target = kwargs.pop('target')
        component_name, port_name = target.split(':')
        component = self.get_component(component_name)

        # Set the properties for the given port
        port = getattr(component.model, port_name)
        port.set_properties(**kwargs)

    def check_parametrized(self):
        """
        Check if all the required parameters have been passed for all components.
        """
        self.flag_not_parametrized = 0  # Set to 1 if at least one component is not fully parametrized
        for component in self.component_list:
            if component.model.parametrized:
                pass
            else:
                self.flag_not_parametrized = 1
                print(f"Error: Parameters of {component.name} are not completely known")

        if self.flag_not_parametrized == 0:
            print("The cycle components are fully parametrized")

    def solve(self):
        # First, check if all the components are parametrized
        self.check_parametrized()
        if self.flag_not_parametrized == 1:
            return

        # Solve the components in the cycle
        for component in self.component_list:
            component.model.check_calculable()
            if component.model.calculable:
                component.model.solve()

    def set_guess(self):
        pass


# Example usage:
if __name__ == "__main__":
    # Create a cycle
    ORC = Cycle()

    # Create components
    PUMP = PumpCstEff()
    EVAP = HXPinchCst()
    COND = HXPinchCst()
    EXP = ExpanderCstEff()

    # Set parameters
    PUMP.set_parameters(eta_is=0.6)

    EVAP.set_parameters(**{
        'Pinch': 5,
        'Delta_T_sh_sc': 5,
        'type_HX': 'evaporator'
    })

    COND.set_parameters(**{
        'Pinch': 5,
        'Delta_T_sh_sc': 5,
        'type_HX': 'condenser'
    })

    EXP.set_parameters(eta_is=0.8)

    # Add components to the cycle
    ORC.add_component(EXP, "Expander")
    ORC.add_component(COND, "Condenser")
    ORC.add_component(PUMP, "Pump")
    ORC.add_component(EVAP, "Evaporator")

    # Link components using specified ports
    ORC.link_components("Pump", "m-ex", "Evaporator", "m-su_wf")
    ORC.link_components("Evaporator", "m-ex_wf", "Expander", "m-su")
    ORC.link_components("Expander", "m-ex", "Condenser", "m-su_wf")
    ORC.link_components("Condenser", "m-ex_wf", "Pump", "m-su")

    # Set the inputs using the new set_properties method
    ORC.set_properties(T=117 + 273.15, fluid='R245fa', m_dot=0.06, target='Condenser:su_wf')
    ORC.set_properties(T=30 + 273.15, fluid='Water', m_dot=0.4, target='Condenser:su_sf')
    ORC.set_properties(cp=4186, target='Condenser:su_sf')
    ORC.set_properties(fluid='Water', target='Condenser:ex_sf')

    ORC.set_properties(fluid='R245fa', target='Pump:su')
    ORC.set_properties(P=880273, fluid='R245fa', target='Pump:ex')

    ORC.set_properties(T=150 + 273.15, fluid='Water', m_dot=0.4, target='Evaporator:su_sf')
    ORC.set_properties(cp=4186, target='Evaporator:su_sf')
    ORC.set_properties(fluid='Water', target='Evaporator:ex_sf')

    # Solve the cycle
    ORC.solve()

    print(ORC.get_component("Expander").model.ex.p)



# class System:
#     def __init__(self):
#         self.components = {}
#         self.connectors = {}

#     def add_component(self, component, name):
#         self.components[name] = component
    
#     def add_connector(self, connector, name):
#         self.connectors[name] = connector
    
#     def set_connector(self, name, **kwargs):
#         self.connectors[name].set_properties(**kwargs)



# class Thermal_System:
#     #-------------------------------------------------------------------------
#     class Cycles():
#         pass
    
#     class Components():
#         pass
    
#     class Components_chain_element():
#         def __init__(self, name, nb_input, nb_output, cycle_names):
#             self.name = name
#             self.previous = [None]*nb_input # will be a thermodynamical point
#             self.next = [None]*nb_output # will be a thermodynamical point
#             self.cycles = cycle_names # Cycle(s) that the component is a part of (s in case of an HTX) 
            
#     class Sources():
#         pass
    
#     class Source_chain_element():
#         def __init__(self, name):
#             self.name = name
#             self.next = [None] # will be a thermodynamical point
            
#     class Sinks():
#         pass
    
#     class Sink_chain_element():
#         def __init__(self, name):
#             self.name = name
#             self.previous = [None] # will be a thermodynamical point
            
#     class Params():
#         pass
#     #-------------------------------------------------------------------------
    
#     def __init__(self, fluid, h_fluid, c_fluid):
#         """
#         Inputs:
#         -------
#         fluid : ORC working fluid
#         h_fluid : Hot source fluid/mixture
#         c_fluid : Cold source fluid/mixture
#         --------------------------------------------------
        
#         An upgrade would be to enter as input an array of fluids in order to have something else than 3 subcycles in the system
#         """
        
#         self.parametrized = False
        
#         self.finished_it = False
#         self.start_key = None
        
#         #self.start_key = {fluid : None, h_fluid : {}, c_fluid : {}}
        
#         self.Params.fluid = fluid # Working Fluid
#         self.Params.h_fluid = h_fluid # Hot source fluid
#         self.Params.c_fluid = c_fluid # Cold source fluid
        
#         self.Components = self.Components()
#         self.Components.dict = {} # Creates a dictionnary that has as keys the name of the components and as values their models
#         self.Components.chains = {fluid : {}, h_fluid : {}, c_fluid : {}} # Create chains that will contain for each cycles the components and thermodynamical points -> Will be iterated on
        
#         self.Sources = self.Sources()
#         self.Sources.dict = {} # Creates a dictionnary to contain all sources (no model or inputs will be passed to them, sources and sinks just ensure that something is there without having a loop)

#         self.Sinks = self.Sinks()
#         self.Sinks.dict = {} # Creates a dictionnary to contain all sources

#         self.Cycles = {fluid : {}, h_fluid : {}, c_fluid : {}} # Will contain only the thermodynamical points related to the cycles 

#     #-------------------------------------------------------------------------
# #%%

#     def add_component(self, component, name, nb_input, nb_output, cycle_names):
#         if name in self.Components.dict:
#             print(f"Component '{name}' already added.")
#             return
        
#         # Add the component to the dictionnary
#         self.Components.dict[name] = component
        
#         # Create the component chain element
#         new_chain_component = self.Components_chain_element(name, nb_input, nb_output, cycle_names)
        
#         for cycle_name in cycle_names:
#             self.Components.chains[cycle_name][name] = new_chain_component
            
#     def link_components(self, name1, output_1, name2, input_2, cycle_name):
                
#         if self.Components.chains[cycle_name][name1].next[output_1-1] != None or self.Components.chains[cycle_name][name2].previous[input_2-1] != None:
#             print('One of ',name1, 'or', name2, 'is already linked at the specified port' )
#             return
#         else:
#             # Create a thermodynamical point and list it in the cycles and component.chains sub_objects
#             new_point_key = cycle_name[0] + str(len(self.Cycles[cycle_name]) + 1)
            
#             self.Cycles[cycle_name][new_point_key] = Mass_connector()
#             self.Components.chains[cycle_name][new_point_key] = self.Cycles[cycle_name][new_point_key] 
            
#             # Link the component 1 to the thermodynamical point
#             self.Components.chains[cycle_name][name1].next[output_1-1] = self.Cycles[cycle_name][new_point_key]
#             self.Cycles[cycle_name][new_point_key].previous = self.Components.chains[cycle_name][name1]
#             self.Components.dict[name1].point_ex[output_1-1] = self.Cycles[cycle_name][new_point_key]
            
#             # Link the thermodynamical point to the component 2
#             self.Components.chains[cycle_name][name2].previous[input_2-1] = self.Cycles[cycle_name][new_point_key]
#             self.Cycles[cycle_name][new_point_key].next = self.Components.chains[cycle_name][name2]
#             self.Components.dict[name2].point_su[input_2-1] = self.Cycles[cycle_name][new_point_key]
            
#             # Set the fluid of the connector
#             self.Cycles[cycle_name][new_point_key].set_fluid(cycle_name)        
    
#                     # oil_source,"Oil_source_1",h_fluid,'Evaporator',2
#     def add_source(self, source, name, cycle_name, name_next_comp, input_port):
#         if name in self.Sources.dict:
#             print(f"Source '{name}' already added.")
#             return
        
#         # Put the source in the source dictionnary
#         self.Sources.dict[name] = source
#         component = self.Components.chains[cycle_name][name_next_comp]

#         # Create the source chain element
#         Source_chain = self.Source_chain_element(name)
        
#         # Create a thermodynamical point and add it to the subcycle
#         new_point_key = cycle_name[0] + str(len(self.Cycles[cycle_name]) + 1)                    
#         self.Cycles[cycle_name][new_point_key] = Mass_connector()
        
#         # Link the source to the thermodynamical point
#         self.Components.chains[cycle_name][new_point_key] = self.Cycles[cycle_name][new_point_key]
#         Source_chain.next = self.Cycles[cycle_name][new_point_key]
#         self.Cycles[cycle_name][new_point_key].previous = Source_chain
        
#         # Link the thermodynamical point to the component after the source
#         component.previous[input_port-1] = self.Cycles[cycle_name][new_point_key]
#         self.Cycles[cycle_name][new_point_key].next = component
#         self.Components.dict[name_next_comp].point_su[input_port-1] = self.Cycles[cycle_name][new_point_key]

#         # Set the fluid of the connector
#         self.Cycles[cycle_name][new_point_key].set_fluid(cycle_name)
        
#     def add_sink(self, sink, name, cycle_name, name_previous_comp,output_port):
#         if name in self.Sinks.dict:
#             print(f"Source '{name}' already added.")
#             return
        
#         # Put the sink in the source dictionnary
#         self.Sinks.dict[name] = sink
#         component = self.Components.chains[cycle_name][name_previous_comp]
        
#         # Create the source chain element   
#         Sink_chain = self.Sink_chain_element(name)
        
#         # Create a thermodynamical point and add it to the subcycle
#         new_point_key = cycle_name[0] + str(len(self.Cycles[cycle_name]) + 1)    
#         self.Cycles[cycle_name][new_point_key] = Mass_connector()
        
#         # Link the thermodynamical point to the sink
#         self.Components.chains[cycle_name][new_point_key] = self.Cycles[cycle_name][new_point_key] 
#         Sink_chain.previous = self.Cycles[cycle_name][new_point_key]
#         self.Cycles[cycle_name][new_point_key].next = Sink_chain
        
#         # Link the the component before the sink to the thermodynamical point
#         component.next[output_port-1] = self.Cycles[cycle_name][new_point_key]
#         self.Cycles[cycle_name][new_point_key].previous = component
#         self.Components.dict[name_previous_comp].point_ex[output_port-1] = self.Cycles[cycle_name][new_point_key]
        
#         # Set the fluid of the connector
#         self.Cycles[cycle_name][new_point_key].set_fluid(cycle_name)        

# #%%

#     def check_parametrized(self):
#         """
#         Check for all components models if all the required parameters have been passed
#         """        
#         flag_not_parametrized = 0 # =1 if at least one component is not fully parametrized
        
#         for comp in self.Components.dict:
#             if self.Components.dict[comp].parametrized:
#                 pass
#             else:
#                 self.defined = False
#                 flag_not_parametrized = 1
#                 print(f"Error: Parameters of {comp} not completely known")
        
#         if flag_not_parametrized == 0:
#             self.parametrized = True
#             print("The cycle components are fully parametrized")

#     def check_connections(self):
#         """
#         Just checks if all components have something in their inputs / outputs
#         Possible upgrade : also check that all components are part of a cycle or are between a source and a sink

#         """
        
#         for key_chain in self.Components.chains:
#             iter_chain = self.Components.chains[key_chain]
            
#             for key_iter_component in iter_chain:
#                 iter_component = iter_chain[key_iter_component]
            
#                 if (not (isinstance(iter_component, self.Sink_chain_element) or isinstance(iter_component, Mass_connector))):
                    
#                     for next_comp in iter_component.next:
#                         if next_comp == None:
#                             print("There is a missing 'next' connection for the component : ", iter_component.name)
#                             self.Params.is_cycle = False
#                             return False
                
#                 if (not (isinstance(iter_component, self.Source_chain_element) or isinstance(iter_component, Mass_connector))):
#                     for prev_comp in iter_component.previous:
#                         if prev_comp == None:
#                             print("There is a missing 'previous' connection for the component : ", iter_component.name)
#                             self.Params.is_cycle = False
#                             return False
        
#         self.Params.is_cycle = True
#         return True

# #%% 

#     def set_su(self, chain_name, comp_name, port, variable, value):
        
#         component_chain_elem = self.Components.chains[chain_name][comp_name]
#         connector = component_chain_elem.previous[port]
        
#         if variable == 'T':
#             connector.set_T(value)            
#         elif variable == 'P':
#             connector.set_p(value)
#         elif variable == 'M_dot':
#             for connector in self.Cycles[chain_name]:
#                 self.Cycles[chain_name][connector].set_m_dot(value)
            
#         return
    
#     def set_ex(self, chain_name, comp_name, port, variable, value):
        
#         component_chain_elem = self.Components.chains[chain_name][comp_name]
#         connector = component_chain_elem.next[port]
        
#         if variable == 'T':
#             connector.set_T(value)            
#         elif variable == 'P':
#             connector.set_p(value)
#         elif variable == 'M_dot':
#             for key in self.Cycles[chain_name]:
#                 connector = self.Cycles[chain_name][key]
#                 connector.set_m_dot(value)
            
#         return

# #%% 
#     def solve_chain_element(self, key_chain_elem, cycle, start_key):
#         print("START KEY : ", start_key)
        
#         if self.finished_it == False: # Check if all components were already computed for this iteration
#             print(key_chain_elem)
#             chain_element = self.Components.chains[cycle][key_chain_elem]
            
#             if isinstance(chain_element, Mass_connector): # We are at a Mass connector
#                 print("Mass connector")
                
#                 if isinstance(chain_element.next, self.Sink_chain_element):
#                     print("Next is a sink")
#                     print("----------------")
#                     return
                
#                 if chain_element.state_known: # State is known, try to solve the next component
#                     print("State known !!!!")
                    
#                     if self.start_key == None: # Real start of the cycle iteration -> saved
#                         self.start_key = key_chain_elem
#                         print("Set of the start key")

#                     print("----------------")

#                     next_comp = chain_element.next # Shall be a component
                    
#                     # Right now I am using the keys to know beter what is happening, it might be possible to only pass the chain_element instead of the key to the method (that would be more efficient)
#                     # Iterate through the items and find the key of the next component
#                     for key, value in self.Components.chains[cycle].items():
#                         if value == next_comp:
#                             next_key = key
#                             break
                        
#                     if self.start_key == None:
#                         self.solve_chain_element(next_key,cycle,start_key)
#                     else:
#                         self.solve_chain_element(next_key,cycle,self.start_key)
                    
#                 else: # State is not known, try to look directly for the next mass connector
#                     print("State not known")
                    
#                     next_key = None
#                     for next_point in chain_element.next.next:
#                         if next_key != None:
#                             break
#                         # Iterate through the items and find the key
#                         for key, value in self.Components.chains[cycle].items():
#                             if value == next_point:
#                                 next_key = key
#                                 break
                            
#                     print("----------------")
                    
#                     if self.start_key == None:
#                         if next_key != start_key:
#                             self.solve_chain_element(next_key,cycle,start_key)
#                         else:
#                             print("Back to start key")
#                             print("!!!!!!!!!!!!!!!!!")
#                             return
#                     else:
#                         if next_key != self.start_key:
#                             self.solve_chain_element(next_key,cycle,self.start_key)   
#                         else:
#                             print("Back to start key")
#                             print("!!!!!!!!!!!!!!!!!")
#                             return
                        
#             else: # We are at a Component
            
#                 next_key = None
#                 for next_point in chain_element.next:
#                     if next_key != None:
#                         break
#                     # Iterate through the items and find the key
#                     for key, value in self.Components.chains[cycle].items():
#                         if value == next_point:
#                             next_key = key
#                             break  
                
#                 if next_key != self.start_key:
#                     comp_model = self.Components.dict[key_chain_elem]
                    
#                     if comp_model.defined == False:
#                         comp_model.solve()
#                         print("Component Solved")
#                         print("????????????????")
                        
#                         for i in range(len(comp_model.point_ex)):
                            
#                             chain_element.next[i].set_p(comp_model.point_ex[i].p)                        
#                             chain_element.next[i].set_h(comp_model.point_ex[i].h)  
                            
#                             # print("- - - - - - - - - - ")
        
#                             # chain_element.next[i].print_resume()
                            
#                             # print("- - - - - - - - - - ")
                        
#                         self.solve_chain_element(next_key,cycle,start_key)
#                     else: 
#                         print("Component already known")
#                         print("----------")
                        
#                 else: # Dernier point à calculer
#                     comp_model = self.Components.dict[key_chain_elem]
#                     comp_model.solve()
                    
#                     print("Component Solved")
#                     print("????????????????")
                        
#                     for i in range(len(comp_model.point_ex)):
                        
#                         chain_element.next[i].set_p(comp_model.point_ex[i].p)                        
#                         chain_element.next[i].set_h(comp_model.point_ex[i].h)  
                        
#                     print("Back to start key")
#                     print("----------")
                    
#                     self.finished_it = True
        
#         else:
#             print("All components are computed")
        
#         return

#     def solve(self):
#         if self.parametrized == True:
#             res = 1e10
#             n_it_max = 1
            
#             tol = 1
#             n_it = 0
            
#             while res > tol and n_it < n_it_max:
#                 print("--------------------------------------\nIteration :",n_it)
#                 for cycle in self.Cycles:
#                     start_key = next(iter(self.Cycles[cycle]))
#                     self.finished_it = False
#                     self.solve_chain_element(start_key, cycle, start_key)
                    
#                     #Search the next component calculable
#                     # if not self.compo_list[i].calculable:
#                     #     self.compo_list[i].check_calculable()

#                     # if self.compo_list[i].calculable and not self.compo_list[i].defined: #If the component is calculable and not defined then we solve it
#                     #     self.compo_list[i].solve()
#                     #     n_solved +=1

#                     # else:
#                     #     pass
                
#                 n_it +=1
#             if res < tol and n_it < n_it_max:
#                 print("Cycle solved")
                
#             else:
#                 print("Error: Cycle couldn't be solved")

#         else:
#             print("Error: Parameters of all components not completely known")

# #%%
    
#     def plot_cycle(self, fluid_name):
#         fig, ax = plt.subplots()
#         ax.axis('off')  # Turn off the axis
    
#         components = self.Components.chains[fluid_name]
#         coordinates = {}
    
#         for i, (comp_name, component) in enumerate(components.items()):
#             x = i * 100
#             y = 0
#             coordinates[comp_name] = (x, y)
#             if isinstance(component, Mass_connector):
#                 ax.plot(x, y, 'o', markersize=12, color='black', fillstyle='none', markeredgecolor='white')  # Represent Mass_connector with circles (black border)
#             else:
#                 ax.text(x, y, comp_name, ha='center', va='center', fontsize=8, bbox=dict(facecolor='white', edgecolor='black'))
    
#         for i, (point_name, point) in enumerate(self.Cycles[fluid_name].items()):
#             x, y = coordinates[point_name]
#             ax.plot(x, y, 'o', markersize=12, color='white', markeredgecolor='black')  # Use markeredgecolor to represent the white border
#             ax.text(x, y, point_name, ha='center', va='center', fontsize=8)
    
#         plt.show()
        
#     def show_cycle_info(self,wf):
        
#         print("-----------------------")
#         print("-----------------------")
#         print("Cycle Results")
#         print("-----------------------")
                
#         for key in self.Cycles[wf]:
#             print(key)
            
#             point = self.Cycles[wf][key]
            
#             print("T :", point.T)
#             print("P :", point.p)
#             print("H :", point.h)
#             print("---------")
            
#             plt.plot(point.h, point.p, marker='o', label=key, color = 'green')
        
#         W_dot_pump = (self.Components.chains[wf]["Pump"].next[0].h - self.Components.chains[wf]["Pump"].previous[0].h)*self.Components.chains[wf]["Pump"].previous[0].m_dot
#         W_dot_exp = (self.Components.chains[wf]["Expander"].previous[0].h - self.Components.chains[wf]["Expander"].next[0].h)*self.Components.chains[wf]["Expander"].previous[0].m_dot
#         Q_dot_evap = (self.Components.chains[wf]["Evaporator"].next[0].h - self.Components.chains[wf]["Evaporator"].previous[0].h)*self.Components.chains[wf]["Evaporator"].previous[0].m_dot
#         Q_dot_cond = (self.Components.chains[wf]["Condenser"].previous[0].h - self.Components.chains[wf]["Condenser"].next[0].h)*self.Components.chains[wf]["Condenser"].previous[0].m_dot
        
#         print("W_dot_pump :",W_dot_pump," [W]")
#         print("W_dot_exp :",W_dot_exp," [W]")
#         print("Q_dot_evap :",Q_dot_evap," [W]")
#         print("Q_dot_cond :",Q_dot_cond," [W]")

#         print("-----------------------")

#         # Saturation curve
#         P_sat = np.linspace(1e5,45e5,1001)
        
#         h_sat1 = np.zeros(len(P_sat))
#         h_sat2 = np.zeros(len(P_sat))
        
#         for i in range(len(h_sat1)):
#             h_sat1[i] = PropsSI('H','P',P_sat[i],'Q',0,'Cyclopentane')
        
#         for i in range(len(h_sat2)):
#             h_sat2[i] = PropsSI('H','P',P_sat[i],'Q',1,'Cyclopentane')
        
#         plt.plot(h_sat1, P_sat, label=key, color = 'black')
#         plt.plot(h_sat2, P_sat, label=key, color = 'black')

#         # Add labels and legend
#         plt.xlabel('Enthalpy (h)')
#         plt.ylabel('Pressure (p)')
#         plt.title('P-h Diagram')
#         plt.grid()
#         plt.show()
        
#         return
        
# #%%

# """
# See for add source_sink and the connection
# """