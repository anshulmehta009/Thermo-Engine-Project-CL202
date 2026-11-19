"""Unit handling (Pint). Internal base units: K, Pa, m^3/kg, kg/m^3, J/kg, J/(kg K)."""
import pint
ureg = pint.UnitRegistry()
BASE = {"temperature": "K", "pressure": "Pa", "volume": "m^3/kg", "density": "kg/m^3",
        "energy": "J/kg", "entropy": "J/(kg*K)"}
UNITS = {
    "SI": {"temperature": ["K", "degC"], "pressure": ["Pa", "kPa", "MPa", "bar"], "volume": ["m^3/kg"],
           "density": ["kg/m^3"], "energy": ["J/kg", "kJ/kg"], "entropy": ["J/(kg*K)", "kJ/(kg*K)"]},
    "CGS": {"temperature": ["K", "degC"], "pressure": ["dyn/cm^2", "bar"], "volume": ["cm^3/g"],
            "density": ["g/cm^3"], "energy": ["erg/g"], "entropy": ["erg/(g*K)"]}}
PROPS = {"T": ("Temperature", "temperature"), "P": ("Pressure", "pressure"), "v": ("Specific volume", "volume"),
         "rho": ("Density", "density"), "u": ("Internal energy", "energy"), "h": ("Enthalpy", "energy"),
         "s": ("Entropy", "entropy"), "A": ("Helmholtz energy", "energy"), "G": ("Gibbs energy", "energy")}

class UnitError(Exception):
    pass

def convert(value, unit, kind, to_base=True):
    """to_base=True: value in `unit` -> base unit. False: base value -> `unit`."""
    try:
        if to_base:
            return float(ureg.Quantity(float(value), unit).to(BASE[kind]).magnitude)
        return float(ureg.Quantity(float(value), BASE[kind]).to(unit).magnitude)
    except (pint.errors.PintError, ValueError, TypeError) as e:
        raise UnitError(str(e))

UNITS["SI"]["temperature"] += ["degF", "degR"]; UNITS["SI"]["pressure"] += ["atm", "psi", "mmHg"]
UNITS["SI"]["volume"] += ["L/kg"]; UNITS["SI"]["density"] += ["g/L"]
UNITS["US Customary"] = {"temperature": ["degF", "degR", "K", "degC"], "pressure": ["psi", "atm", "bar"], "volume": ["ft^3/lb"],
                         "density": ["lb/ft^3"], "energy": ["Btu/lb"], "entropy": ["Btu/(lb*degR)"]}
UNITS["Engineering (metric)"] = {"temperature": ["degC", "K"], "pressure": ["bar", "MPa", "kPa", "atm"], "volume": ["m^3/kg", "L/kg"],
                                 "density": ["kg/m^3"], "energy": ["kJ/kg", "kcal/kg"], "entropy": ["kJ/(kg*K)", "kcal/(kg*K)"]}
