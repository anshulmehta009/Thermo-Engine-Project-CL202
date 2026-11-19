"""Feasibility-based route planning; never silently substitutes a requested route."""
from . import routes, knowledge_base


def analyse(targets, available_inputs, params, T, P, Z, has_cp=True):
    energy = bool(set(targets) & {"h", "s", "u", "A", "G"})
    feasibility = routes.feasibility(params, T, P, Z, has_cp, energy)
    candidates = []
    for route_id, (ok, reason) in feasibility.items():
        candidates.append({"route": route_id, "name": routes.ROUTES[route_id][0],
                           "feasible": bool(ok), "reason": reason,
                           "cost_rank": {1: 1, 3: 2, 2: 3}[route_id]})
    feasible = [c for c in candidates if c["feasible"]]
    selected = min(feasible, key=lambda c: c["cost_rank"]) if feasible else None
    return {"targets": list(targets), "available_inputs": list(available_inputs),
            "dependencies": knowledge_base.dependencies(targets),
            "routes": candidates,
            "recommended_route": selected["route"] if selected else None,
            "selection_basis": "feasibility first; then analytic/closed-form before numerical quadrature"}
