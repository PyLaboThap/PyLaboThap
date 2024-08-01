class Heat_connector:
    def __init__(self):
        """
        Parameters
        ----------
        /
        """
        self.variables_input = []
        
        self.heat = None           # Heat [J]
        self.temperature_in = None # Inlet Temperature [K]
        self.temperature_out = None # Outlet Temperature [K]
        self.UA = None             # Heat transfer coefficient [W/K]
        self.m_dot = None          # Mass flow rate [kg/s]
        self.C_p = None            # Specific heat capacity [J/kg/K]

    def calculate_properties(self):
        if self.variables_input:
            try:
                # Example calculations (for demonstration purposes only)
                if self.m_dot is not None and self.C_p is not None and self.temperature_in is not None and self.temperature_out is not None:
                    self.heat = self.m_dot * self.C_p * (self.temperature_out - self.temperature_in)  # Heat in Joules
                elif self.heat is not None and self.UA is not None and self.temperature_in is not None and self.temperature_out is not None:
                    self.m_dot = self.heat / (self.C_p * (self.temperature_out - self.temperature_in))  # Mass flow rate in kg/s
            except:
                print("Error: Calculation failed")

    def set_heat(self, value):
        self.heat = value
        self.variables_input = self.variables_input + [['heat', value]]
        self.calculate_properties()

    def set_temperature_in(self, value):
        self.temperature_in = value
        self.variables_input = self.variables_input + [['temperature_in', value]]
        self.calculate_properties()

    def set_temperature_out(self, value):
        self.temperature_out = value
        self.variables_input = self.variables_input + [['temperature_out', value]]
        self.calculate_properties()

    def set_UA(self, value):
        self.UA = value
        self.variables_input = self.variables_input + [['UA', value]]
        self.calculate_properties()

    def set_m_dot(self, value):
        self.m_dot = value
        self.variables_input = self.variables_input + [['m_dot', value]]
        self.calculate_properties()

    def set_C_p(self, value):
        self.C_p = value
        self.variables_input = self.variables_input + [['C_p', value]]
        self.calculate_properties()

    def print_resume(self):
        """
        Print a summary of the heat connector properties
        """
        print("Heat: " + str(self.heat) + " [J]")
        print("Inlet Temperature: " + str(self.temperature_in) + " [K]")
        print("Outlet Temperature: " + str(self.temperature_out) + " [K]")
        print("Heat Transfer Coefficient (UA): " + str(self.UA) + " [W/K]")
        print("Mass Flow Rate: " + str(self.m_dot) + " [kg/s]")
        print("Specific Heat Capacity: " + str(self.C_p) + " [J/kg/K]")


# # Example usage
# hc = Heat_connector()
# hc.set_temperature_in(300)  # 300 K
# hc.set_temperature_out(350) # 350 K
# hc.set_m_dot(1)             # 1 kg/s
# hc.set_C_p(4186)            # 4186 J/kg/K (specific heat capacity of water)
# hc.print_resume()
