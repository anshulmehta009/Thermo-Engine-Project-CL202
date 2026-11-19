"""Normalize the application's structured problem request.

This parser handles the explicit form/API payload; it intentionally does not
claim natural-language interpretation. A future NL interface can translate
text into this stable request schema before calling the engine.
"""
from .errors import ThermoError


def parse_request(request):
    if not isinstance(request, dict):
        raise ThermoError("Invalid input", "Problem request must be a dictionary.", "request")
    required = ("fluid", "model", "inputs", "targets", "phase")
    missing = [key for key in required if key not in request]
    if missing:
        raise ThermoError("Invalid input", "Missing request fields: " + ", ".join(missing), "request")
    if not isinstance(request["inputs"], (list, tuple)):
        raise ThermoError("Invalid input", "Inputs must be a list of (property, value, unit) entries.", "inputs")
    normalized = dict(request)
    normalized["fluid"] = str(request["fluid"]).strip()
    normalized["model"] = str(request["model"]).strip()
    normalized["targets"] = list(request["targets"])
    normalized["inputs"] = [tuple(item) for item in request["inputs"]]
    for item in normalized["inputs"]:
        if len(item) != 3:
            raise ThermoError("Invalid input", "Each input must contain property, value and unit.", "inputs")
    return normalized
