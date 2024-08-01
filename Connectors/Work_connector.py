class Work_connector:
    def __init__(self):
        """
        Parameters
        ----------
        /
        """
        self.variables_input = []
        
        self.work = None              # Work [J]
        self.speed = None             # Speed [rpm]
        self.torque = None            # Torque [Nm]
        self.power = None             # Power [W]


    def calculate_properties(self):
        if self.variables_input:
            try:
                # Example calculations (for demonstration purposes only)
                if self.torque is not None and self.speed is not None:
                    self.power = self.torque * (self.speed * 2 * 3.141592653589793 / 60)  # Power in Watts
                elif self.power is not None and self.speed is not None:
                    self.torque = self.power / (self.speed * 2 * 3.141592653589793 / 60)  # Torque in Nm
                elif self.speed is not None and self.work is not None:
                    self.power = self.work / (self.speed * 2 * 3.141592653589793 / 60)  # Power in Watts
            except:
                print("Error: Calculation failed")

    def set_work(self, value):
        self.work = value
        self.variables_input = self.variables_input + [['work', value]]
        self.calculate_properties()

    def set_speed(self, value):
        self.speed = value
        self.variables_input = self.variables_input + [['speed', value]]
        self.calculate_properties()

    def set_torque(self, value):
        self.torque = value
        self.variables_input = self.variables_input + [['torque', value]]
        self.calculate_properties()

    def set_power(self, value):
        self.power = value
        self.variables_input = self.variables_input + [['power', value]]
        self.calculate_properties()

    def set_efficiency(self, value):
        self.efficiency = value
        self.variables_input = self.variables_input + [['efficiency', value]]
        self.calculate_properties()

    def print_resume(self):
        """
        Print a summary of the work connector properties
        """
        print("Work: " + str(self.work) + " [J]")
        print("Speed: " + str(self.speed) + " [rpm]")
        print("Torque: " + str(self.torque) + " [Nm]")
        print("Power: " + str(self.power) + " [W]")


# # Example usage
# wc = Work_connector()
# wc.set_speed(3000)   # 3000 rpm
# wc.set_torque(10)    # 10 Nm
# wc.print_resume()
