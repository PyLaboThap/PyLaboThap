"""
Created on Fri June 19 2024

@author: Elise Neven
@email: elise.neven@uliege.be

"""

from CoolProp.CoolProp import PropsSI
from connector.mass_connector import MassConnector

# Create an instance of the Mass_connector class
point = MassConnector()
point.set_properties(T=500, m_dot=0.5, fluid = 'INCOMP::DowQ', P=101325)
point.print_resume(unit_T='C', unit_p='bar')

point2 = MassConnector()
point2.set_fluid('R1233zd(E)')
point2.set_T(500)
point2.set_cp(2000)
print(point2.cp)
# point.set_T(70+273.15)
# point.print_resume(unit_T='C', unit_p='bar')
# # Set the pressure at 200000 Pa
# point.set_p(200000)
# point.print_resume(unit_T='C', unit_p='bar')


# print(connector.state_known)
# print(connector.completely_known)
# connector.set_V_dot(400)
# connector.print_resume(unit_T='C', unit_p='bar')
# connector.set_m_dot(1)
# connector.print_resume(unit_T='C', unit_p='bar')
# h = PropsSI('H', 'T', 17.91499302642086+273.15, 'P', 1e5, 'R1233zd(E)')


