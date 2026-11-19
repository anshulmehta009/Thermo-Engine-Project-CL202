"""Numerical tolerance and first-order input-uncertainty propagation helpers."""
import math


def validate_tolerances(tol, maxit, residual_tol=None):
    if not math.isfinite(float(tol)) or float(tol) <= 0:
        raise ValueError("Numerical tolerance must be a finite positive number.")
    if int(maxit) < 1:
        raise ValueError("Maximum iterations must be at least 1.")
    if residual_tol is not None and (not math.isfinite(float(residual_tol)) or float(residual_tol) <= 0):
        raise ValueError("Residual tolerance must be a finite positive number.")
    return {"step_tolerance": float(tol), "residual_tolerance": float(residual_tol or tol), "max_iterations": int(maxit)}


def propagate_independent_uncertainties(derivatives, standard_uncertainties):
    """RSS propagation: sigma_y = sqrt(sum((dy/dx_i * sigma_x_i)^2))."""
    terms = {k: float(derivatives[k]) * float(standard_uncertainties[k])
             for k in derivatives if k in standard_uncertainties}
    return {"standard_uncertainty": math.sqrt(sum(x*x for x in terms.values())),
            "contributions": terms,
            "assumption": "independent input uncertainties; first-order linear propagation"}
