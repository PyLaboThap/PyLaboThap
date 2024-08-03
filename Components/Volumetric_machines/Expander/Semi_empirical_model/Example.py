# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:39:37 2023

@author: Elise
"""

from Components.Base_Component import BaseComponent
from Connectors.Mass_connector import Mass_connector
from Connectors.Work_connector import Work_connector
from Connectors.Heat_connector import Heat_connector

from Components.Volumetric_machines.Expander.Semi_empirical_model.Expander_model import Expander

import numpy as np

# Example usage
expander = Expander()
expander.print_setup()

"If the inputs are not set directly BUT throught the connectors"
expander.su.set_fluid('R134a')
expander.ex.set_fluid('R134a')

# # Set properties for su connector
expander.su.set_p(955214.9)
expander.su.set_T(374.1)  # You need to set su.h appropriately

# # Set properties for ex connector
expander.ex.set_p(293940.1)

# Set rotational speed
expander.work_exp.set_speed(1500)

# Set ambient temperature
expander.heat_amb.set_temperature_in(293)

# Setting inputs
# expander.set_inputs(
#     N_rot=3000,
#     T_amb=298.15,
#     su_p=1e5,
#     su_T=300,
#     ex_p=1e4,
#     su_fluid='R134a'  # Make sure to include fluid information
# )

# Setting parameters
expander.set_parameters(
    AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
    A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha=7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712
)

# Solve the expander component
expander.solve()
expander.print_results()

# print(expander.defined)  # Should print True if the component was successfully solved



