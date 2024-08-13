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
    ex_sf_T = 62+273.15
)

EVAP.set_parameters(**{
    'Pinch': 2,
    'Delta_T_sh': 5
})

EVAP.solve()
EVAP.print_results()