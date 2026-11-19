"""Bundled pure-fluid data plus lazy access to the installed CoolProp catalogue.

The application can run its bundled 22-fluid subset without CoolProp. When
CoolProp is installed, additional catalogue fluids are exposed and their
critical constants and ideal-gas Cp are loaded on demand.
"""
from . import coolprop_data


def _f(Tc, Pc, w, M, cp, cool):
    return dict(Tc=Tc, Pc=Pc, w=w, M=M, cp=cp, cool=cool, source="bundled project data")

FLUIDS = {
 "Methane": _f(190.564, 4.5992e6, 0.01142, 16.043, (1.702, 9.081e-3, -2.164e-6, 0), "Methane"),
 "Ethane": _f(305.322, 4.8722e6, 0.0995, 30.07, (1.131, 19.225e-3, -5.561e-6, 0), "Ethane"),
 "Propane": _f(369.89, 4.2512e6, 0.1521, 44.097, (1.213, 28.785e-3, -8.824e-6, 0), "Propane"),
 "n-Butane": _f(425.12, 3.796e6, 0.200, 58.123, (1.935, 36.915e-3, -11.402e-6, 0), "n-Butane"),
 "Isobutane": _f(407.8, 3.64e6, 0.186, 58.123, (1.677, 37.853e-3, -11.945e-6, 0), "IsoButane"),
 "n-Pentane": _f(469.7, 3.370e6, 0.252, 72.15, (2.464, 45.351e-3, -14.111e-6, 0), "n-Pentane"),
 "n-Hexane": _f(507.6, 3.025e6, 0.301, 86.177, (3.025, 53.722e-3, -16.791e-6, 0), "n-Hexane"),
 "Ethylene": _f(282.34, 5.041e6, 0.087, 28.054, (1.424, 14.394e-3, -4.392e-6, 0), "Ethylene"),
 "Propylene": _f(365.57, 4.665e6, 0.140, 42.081, (1.637, 22.706e-3, -6.915e-6, 0), "Propylene"),
 "Acetylene": _f(308.3, 6.139e6, 0.187, 26.038, (6.132, 1.952e-3, 0, -1.299e5), "Acetylene"),
 "Benzene": _f(562.2, 4.898e6, 0.210, 78.114, (-0.206, 39.064e-3, -13.301e-6, 0), "Benzene"),
 "Carbon dioxide": _f(304.1282, 7.3773e6, 0.22394, 44.01, (5.457, 1.045e-3, 0, -1.157e5), "CarbonDioxide"),
 "Carbon monoxide": _f(132.85, 3.494e6, 0.048, 28.01, (3.376, 0.557e-3, 0, -0.031e5), "CarbonMonoxide"),
 "Nitrogen": _f(126.192, 3.3958e6, 0.0372, 28.0134, (3.280, 0.593e-3, 0, 0.040e5), "Nitrogen"),
 "Oxygen": _f(154.581, 5.043e6, 0.0222, 31.999, (3.639, 0.506e-3, 0, -0.227e5), "Oxygen"),
 "Argon": _f(150.687, 4.863e6, -0.002, 39.948, (2.5, 0, 0, 0), "Argon"),
 "Water": _f(647.096, 22.064e6, 0.3443, 18.015, (3.470, 1.450e-3, 0, 0.121e5), "Water"),
 "Ammonia": _f(405.4, 11.333e6, 0.2526, 17.031, (3.578, 3.020e-3, 0, -0.186e5), "Ammonia"),
 "Hydrogen sulfide": _f(373.5, 8.963e6, 0.094, 34.08, (3.931, 1.490e-3, 0, -0.232e5), "HydrogenSulfide"),
 "Sulfur dioxide": _f(430.8, 7.884e6, 0.245, 64.065, (5.699, 0.801e-3, 0, -1.015e5), "SulfurDioxide"),
 "Methanol": _f(512.6, 8.084e6, 0.565, 32.042, (2.211, 12.216e-3, -3.450e-6, 0), "Methanol"),
 "Ethanol": _f(513.9, 6.148e6, 0.645, 46.069, (3.518, 20.001e-3, -6.002e-6, 0), "Ethanol"),
}
CP_RANGE = (250.0, 1500.0)


def list_fluids(include_coolprop=True):
    names = set(FLUIDS)
    if include_coolprop and coolprop_data.COOLPROP_AVAILABLE:
        try:
            names.update(coolprop_data.fluid_names())
        except Exception:
            pass
    return sorted(names)


def get_fluid(name):
    if name in FLUIDS:
        return FLUIDS[name]
    d = coolprop_data.fluid_data(name)
    return {"Tc": d["Tc"], "Pc": d["Pc"], "w": d["omega"],
            "M": d["M_kg_per_mol"] * 1000.0, "cp": None,
            "cool": d["coolprop_name"], "source": "CoolProp fluid database",
            "coolprop_metadata": d}
