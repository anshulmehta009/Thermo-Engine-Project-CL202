"""Thermodynamic equation/knowledge registry used by the route planner."""
PROPERTY_DEPENDENCIES = {
    "v": {"minimum_state": ["T", "P"]}, "rho": {"minimum_state": ["T", "P"]},
    "Z": {"minimum_state": ["T", "P"]},
    "h": {"minimum_state": ["T", "P"], "needs": ["Cp ideal-gas data", "EOS residual enthalpy"]},
    "s": {"minimum_state": ["T", "P"], "needs": ["Cp ideal-gas data", "EOS residual entropy"]},
    "u": {"derived_from": ["h", "P", "v"], "relation": "u = h - P v"},
    "A": {"derived_from": ["u", "T", "s"], "relation": "A = u - T s"},
    "G": {"derived_from": ["h", "T", "s"], "relation": "G = h - T s"},
}
THERMODYNAMIC_RELATIONS = {
    "first_law": "du = T ds - P dv",
    "enthalpy": "h = u + P v",
    "helmholtz": "A = u - T s",
    "gibbs": "G = h - T s",
    "maxwell_1": "(∂s/∂v)_T = (∂P/∂T)_v",
    "maxwell_2": "(∂s/∂P)_T = -(∂v/∂T)_P",
    "cp_cv": "Cp - Cv = T (∂P/∂T)_v (∂v/∂T)_P",
    "dh": "dh = Cp dT + [v - T(∂v/∂T)_P] dP",
    "ds": "ds = Cp/T dT - (∂v/∂T)_P dP",
}
EOS_MODELS = ["Ideal Gas", "vdW", "RK", "SRK", "PR", "PR78"]


def dependencies(targets):
    out = {}
    for target in targets:
        out[target] = PROPERTY_DEPENDENCIES.get(target, {"minimum_state": ["T", "P"]})
    return out
