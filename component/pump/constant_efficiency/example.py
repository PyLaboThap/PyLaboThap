"""
Created on Aug 03 21:31:37 2023

@author: Elise
"""

from component.pump.constant_efficiency.simulation_model import PumpCstEff

import numpy as np

# Example usage
PP = PumpCstEff()

# "If the inputs are not set directly BUT throught the connectors"
PP.su.set_properties(P=319296.5575177148, T=331.033964665788, fluid='R1233ZDE')
PP.ex.set_properties(P=606240.1433176235, fluid='R1233ZDE')
PP.set_parameters(eta_is=0.9)
# PP.print_setup()

PP.solve()
PP.print_results()

# PP.clear_intermediate_states
PP.su.set_properties(P=400000)
PP.ex.set_properties(P=800000, fluid='R1233ZDE')
PP.solve()
PP.print_results()

PP.set_inputs(su_T=400)
PP.solve()
PP.print_results()

# Problème peut être résolu ici! C'est le clear_intermediate_states qui est problématique
