class WorkConnector:
    def __init__(self):
        """
        Parameters
        ----------
        /
        """
        self.variables_input = []
        
        self.w = None              # Specific work [J/kg]
        self.W_dot = None          # Work power [W]
        self.N = None              # Speed [rpm]
        self.T = None            # Torque [Nm]


    def calculate_properties(self):
        pass

    def set_w(self, value):
        self.w = value
        self.variables_input = self.variables_input + [['w', value]]
        self.calculate_properties()

    def set_W_dot(self, value):
        self.W_dot = value
        self.variables_input = self.variables_input + [['W_dot', value]]
        self.calculate_properties()

    def set_N(self, value):
        self.N = value
        self.variables_input = self.variables_input + [['N', value]]
        self.calculate_properties()

    def set_T(self, value):
        self.T = value
        self.variables_input = self.variables_input + [['T', value]]
        self.calculate_properties()


    def print_resume(self):
        """
        Print a summary of the work connector properties
        """
        print("Work: " + str(self.w) + " [J]")
        print("Speed: " + str(self.W_dot) + " [rpm]")
        print("Torque: " + str(self.N) + " [Nm]")
        print("Power: " + str(self.T) + " [W]")
