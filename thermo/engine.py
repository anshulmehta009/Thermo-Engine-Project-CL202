"""Backend entry points: validate -> prepare -> feasibility -> calculate. Structured data only; no UI code."""
from . import units, eos, routes
from .fluids import FLUIDS, CP_RANGE, get_fluid, list_fluids
from .errors import ThermoError
from . import planner, uncertainty, coolprop_data
from .request_parser import parse_request

ENERGY = {"h", "s", "u", "A", "G"}
DEFAULT_NM = dict(root="Brent's method", diff="Central finite difference", integ="Adaptive quadrature",
                  tol=1e-10, residual_tol=1e-8, maxit=100, fd_rel=1e-4, simpson_n=200)
ENGINE_REF = (298.15, 1e5, 0.0, 0.0)   # ideal gas at 298.15 K, 1 bar: h = s = 0

def validate(req):
    req = parse_request(req)
    issues, si = [], {}
    add = lambda c, f, m, sev="error": issues.append(dict(category=c, field=f, message=m, severity=sev))
    import math
    for k, val, unit in req["inputs"]:
        if k not in units.PROPS: add("Invalid input", str(k), "Unknown input property."); continue
        if k in si: add("Invalid input", k, "Property provided more than once."); continue
        try:
            if val is None or not math.isfinite(float(val)):
                add("Invalid input", k, "Value must be a finite number."); continue
        except (TypeError, ValueError):
            add("Invalid input", k, "Value must be a finite number."); continue
        try: x = units.convert(val, unit, units.PROPS[k][1])
        except units.UnitError as e: add("Invalid input", k, f"Unit '{unit}' not usable for {units.PROPS[k][0]}: {e}"); continue
        if k in ("T", "P", "v", "rho") and x <= 0: add("Physically inadmissible state", k, f"{units.PROPS[k][0]} must be > 0 (absolute scale)."); continue
        si[k] = x
    p = None
    try: p = eos.make_params(req["fluid"], req["model"])
    except ThermoError as e: add(e.category, "model/fluid", e.message); p = None
    tg = set(req["targets"])
    if not tg: add("Invalid input", "targets", "Select at least one target property.")
    if tg & ENERGY and p and p.cp_func is None:
        add("Missing parameters", "fluid", f"Ideal-gas Cp data are not available for {req['fluid']}; energy-type targets unsupported.")
    if "T" not in si or not ({"P", "v", "rho"} & si.keys()):
        add("Insufficient independent properties", "inputs", "Need T plus one of P, v, rho to fix the state (single-phase pure fluid).")
    for k in ENERGY & si.keys():
        add("Insufficient independent properties", k, f"Inverse state determination from {units.PROPS[k][0]} is not implemented; use T with P, v or rho.")
    if "rho" in si:
        v = 1 / si["rho"]
        if "v" in si and abs(v - si["v"]) > 1e-6 * v: add("Physically inadmissible state", "v/rho", "v and rho are inconsistent.")
        si["v"] = v
    if "T" in tg and "T" not in si: add("Insufficient independent properties", "targets", "T as a target requires inverse solving; not implemented.")
    if "P" in tg and "P" not in si and "v" not in si: add("Insufficient independent properties", "targets", "P target requires T and v (or rho).")
    if p and not any(i["severity"] == "error" for i in issues):
        T = si["T"]
        if "v" in si and si["v"] <= p.b: add("Physically inadmissible state", "v", f"v <= EOS co-volume b = {p.b:.4e} m^3/kg.")
        elif "P" in si and "v" in si:
            Pe = eos.pressure(p, T, si["v"])
            if abs(Pe - si["P"]) > 1e-6 * si["P"]: add("Physically inadmissible state", "P/v", f"P and v inconsistent with {p.model}: P_EOS(T,v) = {Pe:.6g} Pa.")
        ph, P = req["phase"], si.get("P")
        if P and ph == "Supercritical" and not (T > p.Tc and P > p.Pc): add("Physically inadmissible state", "phase", "Supercritical needs T>Tc and P>Pc.")
        if ph == "Liquid" and T >= p.Tc: add("Physically inadmissible state", "phase", "Liquid requested above Tc.")
        if tg & ENERGY and p.cp is not None and not CP_RANGE[0] <= T <= CP_RANGE[1]:
            add("Physically inadmissible state", "T", f"Bundled ideal-gas Cp polynomial is documented for approximately {CP_RANGE[0]:.0f}-{CP_RANGE[1]:.0f} K; verify its range.", "warning")
    return dict(ok=not any(i["severity"] == "error" for i in issues), issues=issues, si=si)

def prepare(req, si, nm):
    p = eos.make_params(req["fluid"], req["model"]); T = si["T"]; tr = []
    tr.append(("Converted inputs (SI base units)", si))
    tr.append(("EOS parameters", dict(model=p.model, a=p.a, b=p.b, Tc=p.Tc, Pc=p.Pc, omega=p.omega, kappa=p.k, Rs=p.Rs)))
    tr.append(("EOS equations used", {"latex": eos.EOS_LATEX[p.model]}))
    tr.append(("Fluid data provenance", {"source": p.data_source, "CoolProp_available": coolprop_data.COOLPROP_AVAILABLE}))
    if "P" in si:
        P = si["P"]; cands = eos.z_candidates(p, T, P, nm["root"], nm["tol"], nm["maxit"], nm.get("residual_tol"))
        best = eos.select_root(p, T, P, cands, req["phase"])
        tr.append(("Root finding & filtering (cubic in Z)", [{k: c.get(k) for k in ("Z", "status", "reason", "residual", "iters", "evals", "label")} for c in cands]))
        if best is None:
            fails = [c for c in cands if c["status"] == "numerical failure"]
            raise ThermoError("Numerical non-convergence" if fails else "Physically inadmissible state",
                              "No root survives filtering: " + "; ".join(f"Z={c.get('Z')}: {c['reason']}" for c in cands))
        Z, A, B = best["Z"], best["A"], best["B"]; v = Z * p.Rs * T / P; label = best["label"]; res = best["residual"]
    else:
        v = si["v"]; P = eos.pressure(p, T, v); _, A, B = eos.z_poly(p, T, P); Z = P * v / (p.Rs * T); cands = []; res = None
        label = eos.label_state(p, T, P, Z, None)
        tr.append(("Pressure from EOS (T,v given; no root solving needed)", dict(P=P, Z=Z)))
    pr = abs(eos.pressure(p, T, v) - P) / P
    plan = planner.analyse(req["targets"], list(si), p, T, P, Z, p.cp_func is not None)
    tr.append(("Property dependencies & route analysis", plan))
    return dict(p=p, T=T, P=P, v=v, Z=Z, A=A, B=B, label=label, cands=cands, trace=tr, z_residual=res, P_rel_residual=pr, plan=plan)

def calculate(req, st, route, nm, ref=ENGINE_REF, cross_check=False):
    try:
        uncertainty.validate_tolerances(nm.get("tol", 1e-10), nm.get("maxit", 100), nm.get("residual_tol", None))
    except (TypeError, ValueError) as exc:
        raise ThermoError("Invalid input", str(exc), "tolerance")
    p = st["p"]; tg = set(req["targets"]); T, P, v = st["T"], st["P"], st["v"]
    if len(ref) != 4 or not all(__import__("math").isfinite(float(x)) for x in ref) or ref[0] <= 0 or ref[1] <= 0:
        raise ThermoError("Invalid input", "Reference state requires finite values and positive T_ref and P_ref.", "reference state")
    out = dict(T=T, P=P, v=v, rho=1 / v, Z=st["Z"]); trace = list(st["trace"]); info = {}
    if tg & ENERGY:
        feas = routes.feasibility(p, T, P, st["Z"], p.cp_func is not None, True)
        if route not in feas or not feas[route][0]:
            raise ThermoError("Unsupported fluid-model combination", f"Route {route} unsupported here: {feas.get(route, (0, 'unknown'))[1]}. No substitution performed.")
        r = routes.ROUTES[route][1](p, st, ref, nm); h = r["h_ig"] + r["h_res"]; s = r["s_ig"] + r["s_res"]
        u = h - P * v; out.update(h=h, s=s, u=u, A=u - T * s, G=h - T * s)
        info = dict(route=route, route_name=routes.ROUTES[route][0], evals=r["evals"], integration_error=r["error"],
                    parts={k: r[k] for k in ("h_ig", "s_ig", "h_res", "s_res")})
        trace += [("Selected route", info["route_name"]), ("Equations evaluated", {"latex": routes.EQUATIONS[route] + routes.DERIVED}), ("Numerical settings", nm), ("Ideal + residual parts", info["parts"])]
        if cross_check:
            cc = {}
            for q, (ok, _) in feas.items():
                if q != route and ok:
                    r2 = routes.ROUTES[q][1](p, st, ref, nm)
                    cc[f"route {q}"] = dict(dh=(r2["h_ig"] + r2["h_res"]) - h, ds=(r2["s_ig"] + r2["s_res"]) - s)
            info["cross_check"] = cc
            trace.append(("Cross-route comparison (internal consistency, NOT independent verification)", cc))
    trace.append(("Verification", dict(EOS_pressure_relative_residual=st["P_rel_residual"], cubic_residual=st["z_residual"],
                                       integration_error_estimates=info.get("integration_error", "not applicable / closed form"))))
    return dict(values=out, info=info, trace=trace, label=st["label"], nm=nm, ref=ref)

def coolprop_reference(fluid, T, P):
    """Independent CoolProp comparison; never overwrites the selected EOS result."""
    try:
        return coolprop_data.props(fluid, T, P)
    except Exception as e:
        return dict(source="CoolProp", error=str(e))


def fluid_database(fluid=None):
    """Return CoolProp catalogue and, optionally, detailed raw fluid metadata."""
    if fluid is None:
        names = list_fluids()
        return {"coolprop_available": coolprop_data.COOLPROP_AVAILABLE, "count": len(names), "fluids": names}
    try:
        return coolprop_data.fluid_data(fluid)
    except Exception as exc:
        return {"error": str(exc), "bundled_record": FLUIDS.get(fluid)}


def isotherm(st, n=250):
    """P-v curves (state isotherm and critical isotherm) from the EOS, plus the state point."""
    import numpy as np, pandas as pd
    p = st["p"]; vc = p.Zc * p.Rs * p.Tc / p.Pc; rows = []
    for nm_, T in ((f"T = {st['T']:.4g} K", st["T"]), (f"Tc = {p.Tc:.4g} K", p.Tc)):
        for v in np.logspace(np.log10(p.b * 1.03), np.log10(max(st["v"], vc) * 30), n):
            P = eos.pressure(p, T, v)
            if P > 0: rows.append(dict(v=v, P=P, curve=nm_))
    return pd.DataFrame(rows), pd.DataFrame([dict(v=st["v"], P=st["P"])])
