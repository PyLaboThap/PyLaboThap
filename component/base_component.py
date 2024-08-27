class BaseComponent:
    def __init__(self):
        self.calculable = False
        self.parametrized = False
        self.solved = False
        self.inputs = {}
        self.params = {}
        self.guesses = {}

    def set_inputs(self, **inputs):
        for key, value in inputs.items():
            self.inputs[key] = value

    def set_parameters(self, **parameters):
        for key, value in parameters.items():
            self.params[key] = value
        self.check_parametrized()

    def set_guesses(self, **guesses):
        for key, value in guesses.items():
            self.guesses[key] = value

    def check_calculable(self):
        required_inputs = self.get_required_inputs() 
        self.calculable = all(self.inputs.get(inp) is not None for inp in required_inputs) # check if all required inputs are set
        return self.calculable

    def check_parametrized(self):
        required_params = self.get_required_parameters()
        self.parametrized = all(self.params.get(param) is not None for param in required_params) # check if all required parameters are set
        return self.parametrized

    def get_required_inputs(self):
        # This method should be overridden in derived classes
        return []

    def get_required_parameters(self):
        # This method should be overridden in derived classes
        return []
    
    def get_required_guesses(self):

        return []

    def solve(self):
        # This method should be overridden in derived classes
        raise NotImplementedError("The 'solve' method should be implemented in derived classes.")
