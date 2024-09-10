class HeatConnector:
    """
    A class to handle transfer of heat power.

    Attributes:
    ^^^^^^^^^^^
    Q_dot : float, optional
        Heat power in W.
    variables_input : list of lists
        A list of the variables used to define the heat connector. Each entry is a list of [variable_name, value].

    Methods:
    ^^^^^^^^
    __init__(self):
        Initializes the HeatConnector object with.

    set_Q_dot(self, value):
        Sets the heat power and updates the list of known variables.

    print_resume(self):
        Print a summary of the heat connector properties.
    """

    def __init__(self):

        self.variables_input = []
        
        self.Q_dot = None             # Heat power [W]

    def calculate_properties(self):
        pass

    def set_Q_dot(self, value):
        self.Q_dot = value
        self.variables_input = self.variables_input + [['Q_dot', value]]
        self.calculate_properties()


    def print_resume(self):
        """
        Print a summary of the heat connector properties
        """
        print("Heat power: " + str(self.Q_dot) + " [W]")

