# -*- coding: utf-8 -*-
"""
Created on Aug 03 21:31:37 2023

@author: Elise
"""

from Components.Base_Component import BaseComponent
from Connectors.Mass_connector import Mass_connector
from Connectors.Work_connector import Work_connector
from Connectors.Heat_connector import Heat_connector

from Components.Volumetric_machines.Expander.Constant_isentropic_efficiency.Expander_model import Expander_CstEff

import numpy as np

# Example usage
EXP = Expander_CstEff()
EXP.print_setup()

# "If the inputs are not set directly BUT throught the connectors"
EXP.su.set_properties(P=955214.9, T=374.18, fluid='R134a')
EXP.ex.set_properties(P=293940.1, fluid='R1233ZDE')
EXP.set_parameters(eta_is=0.8)
EXP.print_setup()

EXP.solve()
EXP.print_results()
