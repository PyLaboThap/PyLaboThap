"""
Created on Aug 03 21:31:37 2023

@author: Elise
"""

from Components.Base_Component import BaseComponent
from Connectors.Mass_connector import Mass_connector
from Connectors.Work_connector import Work_connector
from Connectors.Heat_connector import Heat_connector

from Components.Volumetric_machines.Compressor.Constant_isentropic_efficiency.Compressor_model import Compressor_CstEff

import numpy as np

# Example usage
CP = Compressor_CstEff()
CP.print_setup()

# "If the inputs are not set directly BUT throught the connectors"
CP.su.set_properties(P=319296.5575177148, T=331.033964665788, fluid='R1233ZDE')
CP.ex.set_properties(P=606240.1433176235, fluid='R1233ZDE')
CP.set_parameters(eta_is=0.8)
CP.print_setup()

CP.solve()
CP.print_results()