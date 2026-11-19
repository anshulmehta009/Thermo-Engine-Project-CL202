"""SymPy utilities for symbolic thermodynamic manipulation."""
import sympy as sp

T, v, P, R, a, b = sp.symbols("T v P R a b", positive=True)


def differentiate(expression, variable):
    return sp.diff(expression, variable)


def integrate(expression, variable):
    return sp.integrate(expression, variable)


def simplify(expression):
    return sp.simplify(expression)


def maxwell_pressure_derivative(eos_expression):
    """Return (∂P/∂T)_v and (∂P/∂v)_T for a SymPy EOS expression."""
    return sp.diff(eos_expression, T), sp.diff(eos_expression, v)
