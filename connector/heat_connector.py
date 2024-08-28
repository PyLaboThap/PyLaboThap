class HeatConnector:
    def __init__(self):
        """
        Parameters
        ----------
        /
        """
        self.variables_input = []
        
        self.q = None                 # specific heat [J/kg]
        self.Q_dot = None             # Heat power [W]

    def calculate_properties(self):
        pass

    def set_q(self, value):
        self.q = value
        self.variables_input = self.variables_input + [['q', value]]
        self.calculate_properties()

    def set_Q_dot(self, value):
        self.Q_dot = value
        self.variables_input = self.variables_input + [['Q_dot', value]]
        self.calculate_properties()


    def print_resume(self):
        """
        Print a summary of the heat connector properties
        """
        print("Heat: " + str(self.q) + " [J]")
        print("Inlet Temperature: " + str(self.Q_dot) + " [W]")

