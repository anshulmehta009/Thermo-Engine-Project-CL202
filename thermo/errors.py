CATEGORIES = ["Invalid input", "Insufficient independent properties", "Unsupported fluid-model combination",
              "Missing parameters", "Mathematical failure", "Physically inadmissible state", "Numerical non-convergence"]

class ThermoError(Exception):
    def __init__(self, category, message, field=None):
        super().__init__(message)
        self.category, self.message, self.field = category, message, field
