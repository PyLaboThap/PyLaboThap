from component.heat_exchanger.pinch_cst.simulation_model import HXPinchCst
import numpy as np


"Evaporator test"

#Exo ORC M&S
EVAP = HXPinchCst()

EVAP.set_inputs(
    fluid_wf = 'R245fa',
    su_wf_T = 60+273.15,
    su_wf_m_dot = 0.4621,
    fluid_sf = 'Water',
    su_sf_T = 70+273.15,
    su_sf_cp = 4186,
    su_sf_m_dot = 2.425,
)

EVAP.set_parameters(**{
    'Pinch': 2,
    'Delta_T_sh_sc': 5,
    'type_HX': 'evaporator'
})

EVAP.solve()
EVAP.print_results()

# "Condenser test"

# COND = HXPinchCst()

# COND.set_inputs(
#     fluid_wf = 'R245fa',
#     su_wf_T = 117+273.15,
#     su_wf_m_dot = 0.4,
#     fluid_sf = 'Water',
#     su_sf_T = 30+273.15,
#     su_sf_cp = 4186,
#     su_sf_m_dot = 0.4
# )

# COND.set_parameters(**{
#     'Pinch': 2,
#     'Delta_T_sh_sc': 5,
#     'type_HX': 'condenser'
# })

# COND.solve()
# COND.print_results()
# print(COND.ex_wf.p)