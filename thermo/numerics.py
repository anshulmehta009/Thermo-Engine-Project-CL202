"""Numerical methods used by the thermodynamic engine.

Root discovery: NumPy polynomial roots. Refinement: custom bisection/Newton
and SciPy Brent. Derivatives: central/forward finite differences. Integrals:
SciPy adaptive Gauss-Kronrod quadrature or composite Simpson.
"""
import numpy as np
from scipy import integrate, optimize


def bisection(f, a, b, tol=1e-12, maxit=200):
    fa, fb = f(a), f(b)
    if not np.isfinite(fa) or not np.isfinite(fb) or fa * fb > 0:
        raise ValueError("bracket must have finite endpoint values and contain a sign change")
    n, m = 0, a
    while n < maxit:
        m = 0.5 * (a + b); fm = f(m); n += 1
        if not np.isfinite(fm):
            raise ValueError("non-finite function value during bisection")
        if fm == 0 or 0.5 * abs(b - a) < tol:
            return dict(root=m, iters=n, evals=n + 2, converged=True, method="Bisection")
        if fa * fm < 0: b, fb = m, fm
        else: a, fa = m, fm
    return dict(root=m, iters=n, evals=n + 2, converged=False, method="Bisection")


def brent(f, a, b, tol=1e-12, maxit=200):
    r, info = optimize.brentq(f, a, b, xtol=tol, maxiter=maxit, full_output=True, disp=False)
    return dict(root=r, iters=info.iterations, evals=info.function_calls,
                converged=bool(info.converged), method="SciPy Brent")


def newton(f, df, x, tol=1e-12, maxit=100):
    for n in range(1, maxit + 1):
        fx, d = f(x), df(x)
        if not np.isfinite(d) or d == 0: raise ValueError("zero or non-finite derivative in Newton-Raphson")
        dx = fx / d; x -= dx
        if abs(dx) < tol and abs(f(x)) <= max(tol, 1e-14):
            return dict(root=x, iters=n, evals=2*n+1, converged=True, method="Newton-Raphson")
    return dict(root=x, iters=maxit, evals=2*maxit+1, converged=False, method="Newton-Raphson")


def secant(f, x0, x1, tol=1e-12, maxit=100):
    info = optimize.root_scalar(f, x0=x0, x1=x1, method="secant", xtol=tol, maxiter=maxit)
    return dict(root=float(info.root), iters=info.iterations, evals=info.function_calls,
                converged=bool(info.converged), method="SciPy secant")


def fd_central(f, x, h):
    if h <= 0: raise ValueError("finite-difference step must be positive")
    return (f(x + h) - f(x - h)) / (2 * h)


def fd_forward(f, x, h):
    if h <= 0: raise ValueError("finite-difference step must be positive")
    return (f(x + h) - f(x)) / h


def fd_five_point(f, x, h):
    """Fourth-order central finite difference."""
    if h <= 0: raise ValueError("finite-difference step must be positive")
    return (f(x - 2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h)) / (12*h)


def integrate_1d(f, a, b, method="Adaptive quadrature", tol=1e-10, n=200):
    """Return (integral, estimated absolute error, function evaluations)."""
    if not np.isfinite(a) or not np.isfinite(b):
        if method != "Adaptive quadrature":
            raise ValueError("Simpson integration requires finite limits")
    if method == "Adaptive quadrature":
        val, err, info = integrate.quad(f, a, b, epsabs=tol, epsrel=tol, limit=500, full_output=1)[:3]
        return float(val), float(err), int(info["neval"])
    if method != "Simpson's rule":
        raise ValueError(f"Unknown integration method: {method}")
    n = max(2, int(n)); n += n % 2
    x = np.linspace(a, b, n + 1); y = np.array([f(t) for t in x])
    fine = integrate.simpson(y, x=x)
    coarse = integrate.simpson(y[::2], x=x[::2])
    return float(fine), float(abs(fine - coarse) / 15.0), int(n + 1)
