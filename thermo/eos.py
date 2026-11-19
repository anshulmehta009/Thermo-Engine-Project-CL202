"""Own PR / SRK implementation. P = R T/(v-b) - a*alpha(T)/(v^2 + u b v + w b^2), Soave alpha."""
import functools
from dataclasses import dataclass
import numpy as np, sympy as sp
from .fluids import FLUIDS, get_fluid
from .errors import ThermoError
from . import numerics as nx

RU = 8.314462618
def _pr78(w): return 0.37464 + 1.54226*w - 0.26992*w*w if w <= 0.491 else 0.379642 + 1.48503*w - 0.164423*w*w + 0.016666*w**3
MODELS = {
 "vdW": dict(name="van der Waals", u=0, w=0, oa=27/64, ob=1/8, Zc=0.375, alpha="const", kap=lambda w: 0.0, info="Original cubic EOS; qualitative only, large quantitative errors."),
 "RK": dict(name="Redlich–Kwong", u=1, w=0, oa=0.42748, ob=0.08664, Zc=1/3, alpha="rk", kap=lambda w: 0.0, info="alpha = Tr^-1/2; good for gases, weak for liquids."),
 "SRK": dict(name="Soave–Redlich–Kwong", u=1, w=0, oa=0.42748, ob=0.08664, Zc=1/3, alpha="soave", kap=lambda w: 0.480 + 1.574*w - 0.176*w*w, info="Soave alpha with acentric factor; widely used for hydrocarbons."),
 "PR": dict(name="Peng–Robinson (1976)", u=2, w=-1, oa=0.45724, ob=0.07780, Zc=0.3074, alpha="soave", kap=lambda w: 0.37464 + 1.54226*w - 0.26992*w*w, info="Better liquid densities than SRK."),
 "PR78": dict(name="Peng–Robinson (1978)", u=2, w=-1, oa=0.45724, ob=0.07780, Zc=0.3074, alpha="soave", kap=_pr78, info="PR with extended kappa for heavy fluids (omega > 0.491)."),
}
_C = [r"A=\frac{a\alpha P}{R^{2}T^{2}},\quad B=\frac{bP}{RT},\quad Z=\frac{Pv}{RT}", r"(Z-B-1)\,(Z^{2}+uBZ+wB^{2})+A\,(Z-B)=0"]
_G = r"P=\frac{RT}{v-b}-\frac{a\,\alpha(T)}{v^{2}+u\,b\,v+w\,b^{2}}"
EOS_LATEX = {
 "vdW": [_G, r"u=w=0,\; a=\tfrac{27}{64}\frac{R^{2}T_c^{2}}{P_c},\; b=\tfrac18\frac{RT_c}{P_c},\; \alpha=1"] + _C,
 "RK": [_G, r"u=1,w=0,\; a=0.42748\frac{R^{2}T_c^{2}}{P_c},\; b=0.08664\frac{RT_c}{P_c},\; \alpha=T_r^{-1/2}"] + _C,
 "SRK": [_G, r"u=1,w=0,\; \alpha=[1+\kappa(1-\sqrt{T_r})]^{2},\; \kappa=0.480+1.574\omega-0.176\omega^{2}"] + _C,
 "PR": [_G, r"u=2,w=-1,\; a=0.45724\frac{R^{2}T_c^{2}}{P_c},\; b=0.07780\frac{RT_c}{P_c}", r"\alpha=[1+\kappa(1-\sqrt{T_r})]^{2},\; \kappa=0.37464+1.54226\omega-0.26992\omega^{2}"] + _C,
 "PR78": [_G, r"\kappa=0.379642+1.48503\omega-0.164423\omega^{2}+0.016666\omega^{3}\;(\omega>0.491)"] + _C,
}

@functools.lru_cache(None)
def sym(model):
    m = MODELS[model]
    T, v, a, b, Tc, k, R = sp.symbols("T v a b Tc k R", positive=True)
    al = {"soave": (1 + k * (1 - sp.sqrt(T / Tc))) ** 2, "rk": sp.sqrt(Tc / T), "const": sp.Integer(1)}[m["alpha"]]
    P = R * T / (v - b) - a * al / (v**2 + m["u"] * b * v + m["w"] * b**2)
    L = lambda e: sp.lambdify((T, v, a, b, Tc, k, R), e, "numpy")
    aa = a * al
    return dict(P=L(P), dPdT=L(sp.diff(P, T)), dPdv=L(sp.diff(P, v)), aa=L(aa), daa=L(sp.diff(aa, T)),
                eu=L(sp.simplify(P - T * sp.diff(P, T))),   # u_res = int_v^inf eu dv
                es=L(sp.simplify(R / v - sp.diff(P, T))))   # s_res(T,v) = int_v^inf es dv

@dataclass
class Params:
    model: str; fluid: str; u: float; w: float; a: float; b: float; Tc: float; Pc: float
    k: float; omega: float; Rs: float; cp: tuple; M: float; Zc: float; cp_func: object = None; data_source: str = "bundled project data"

def make_params(fluid, model):
    if model not in MODELS:
        raise ThermoError("Unsupported fluid-model combination", f"EOS '{model}' is not implemented.")
    # Resolve both bundled records and fluids supplied by the installed CoolProp
    # catalogue. The UI intentionally exposes the full catalogue when available;
    # rejecting non-bundled names here made those selections impossible to run.
    try:
        f = get_fluid(fluid)
    except Exception as exc:
        raise ThermoError(
            "Unsupported fluid-model combination",
            f"Could not load fluid '{fluid}'. Check that it is a supported pure fluid "
            f"and that CoolProp is installed for non-bundled fluids. Details: {exc}",
        )
    required = ("Tc", "Pc", "w", "M")
    missing = [key for key in required if key not in f or f[key] is None]
    if missing:
        raise ThermoError(
            "Missing parameters",
            f"Fluid '{fluid}' is missing required EOS constants: {', '.join(missing)}.",
        )
    if not all(np.isfinite(float(f[key])) for key in required) or f["Tc"] <= 0 or f["Pc"] <= 0 or f["M"] <= 0:
        raise ThermoError(
            "Missing parameters",
            f"Fluid '{fluid}' has invalid critical constants, acentric factor or molar mass.",
        )
    m = MODELS[model]
    Rs = RU / f["M"] * 1000.0
    cp_tuple = f.get("cp")
    if cp_tuple is not None:
        cp_func = lambda T: Rs * (cp_tuple[0] + cp_tuple[1]*T + cp_tuple[2]*T*T + cp_tuple[3]/T**2)
    elif f.get("source") == "CoolProp fluid database":
        from .coolprop_data import ideal_cp_mass

        def cp_func(T):
            try:
                return ideal_cp_mass(f["cool"], T)
            except Exception as exc:
                raise ThermoError(
                    "Missing parameters",
                    f"CoolProp could not evaluate ideal-gas Cp for '{fluid}' at "
                    f"T={float(T):.6g} K. This fluid may not expose CP0MASS over "
                    f"the requested range. Details: {exc}",
                    "Cp",
                ) from exc
    else:
        cp_func = None
    return Params(model, fluid, m["u"], m["w"], m["oa"] * Rs**2 * f["Tc"]**2 / f["Pc"], m["ob"] * Rs * f["Tc"] / f["Pc"],
                  f["Tc"], f["Pc"], m["kap"](f["w"]), f["w"], Rs, cp_tuple, f["M"], m["Zc"], cp_func, f.get("source", "unknown"))

def _args(p, T, v): return (T, v, p.a, p.b, p.Tc, p.k, p.Rs)
def pressure(p, T, v): return float(sym(p.model)["P"](*_args(p, T, v)))
def dPdv(p, T, v): return float(sym(p.model)["dPdv"](*_args(p, T, v)))

def z_poly(p, T, P):
    aa = float(sym(p.model)["aa"](*_args(p, T, 1.0)))
    A, B = aa * P / (p.Rs * T)**2, p.b * P / (p.Rs * T)
    c = np.polyadd(np.polymul([1, -(B + 1)], [1, p.u * B, p.w * B * B]), A * np.array([1, -B]))
    return c, A, B

def ln_phi(p, Z, A, B):
    sq = np.sqrt(p.u**2 - 4 * p.w)
    if sq == 0: return Z - 1 - np.log(Z - B) - A / Z   # van der Waals
    return Z - 1 - np.log(Z - B) - A / (B * sq) * np.log((2 * Z + B * (p.u + sq)) / (2 * Z + B * (p.u - sq)))

def zmax(p, T, P):
    c, A, B = z_poly(p, T, P); f = np.poly1d(c); df = f.deriv()
    z = max(r.real for r in np.roots(c) if abs(r.imag) < 1e-9 * max(1, abs(r)))
    for _ in range(3): z -= f(z) / df(z)
    return z

def z_candidates(p, T, P, method="Brent's method", tol=1e-12, maxit=100, residual_tol=None):
    """All 3 mathematical roots, each refined with the chosen solver and classified."""
    c, A, B = z_poly(p, T, P); f = np.poly1d(c); df = f.deriv(); out = []
    for r in np.roots(c):
        d = dict(raw=complex(r), A=A, B=B)
        if abs(r.imag) > 1e-9 * max(1, abs(r)):
            out.append({**d, "Z": None, "status": "mathematically invalid", "reason": "complex root"}); continue
        z0 = r.real
        try:
            if method == "Newton-Raphson": info = nx.newton(f, df, z0, tol, maxit)
            elif method == "Bisection": info = nx.bisection(f, z0 * (1 - 1e-4), z0 * (1 + 1e-4), tol, maxit)
            elif method == "Secant method": info = nx.secant(f, z0 * (1 - 1e-4), z0 * (1 + 1e-4), tol, maxit)
            else: info = nx.brent(f, z0 * (1 - 1e-4), z0 * (1 + 1e-4), tol, maxit)
        except ValueError as e:
            out.append({**d, "Z": z0, "status": "numerical failure", "reason": str(e)}); continue
        Z = info["root"]; d.update(Z=Z, residual=abs(f(Z)), iters=info["iters"], evals=info["evals"], converged=info["converged"])
        v = Z * p.Rs * T / P
        rtol = float(residual_tol if residual_tol is not None else max(tol, 1e-10))
        if not info["converged"]: d.update(status="numerical failure", reason="root solver did not converge within max iterations")
        elif abs(f(Z)) > rtol: d.update(status="numerical failure", reason=f"cubic residual {abs(f(Z)):.3e} exceeds residual tolerance {rtol:.3e}")
        elif Z <= B: d.update(status="physically inadmissible", reason="Z <= B, i.e. v <= b (volume below co-volume)")
        elif dPdv(p, T, v) >= 0: d.update(status="physically inadmissible", reason="mechanically unstable: (dP/dv)_T >= 0")
        else: d.update(status="retained", reason="real, v>b, mechanically stable")
        out.append(d)
    return out

def label_state(p, T, P, Z, multi_pos):
    """multi_pos 'low'/'high' when >=2 stable roots exist; single-root labelling is a density heuristic."""
    if multi_pos: return "liquid" if multi_pos == "low" else "vapour"
    if T >= p.Tc: return "supercritical" if P >= p.Pc else "vapour"
    return "liquid" if Z * p.Rs * T / P < p.Zc * p.Rs * p.Tc / p.Pc else "vapour"

def select_root(p, T, P, cands, phase):
    ret = sorted([c for c in cands if c["status"] == "retained"], key=lambda c: c["Z"])
    for i, c in enumerate(ret):
        pos = None if len(ret) < 2 else ("low" if i == 0 else "high" if i == len(ret) - 1 else None)
        c["label"] = label_state(p, T, P, c["Z"], pos)
    want = phase.lower().replace("vapor", "vapour")
    if want != "any phase":
        for c in ret:
            if c["label"] != want: c["status"], c["reason"] = "inconsistent with selected phase", f"root is {c['label']}, requested {want}"
        ret = [c for c in ret if c["status"] == "retained"]
    if not ret: return None
    for c in ret: c["ln_phi"] = float(ln_phi(p, c["Z"], c["A"], c["B"]))
    best = min(ret, key=lambda c: c["ln_phi"])   # lowest Gibbs energy at fixed T,P
    best["status"] = "selected"; best["reason"] = "lowest fugacity coefficient (minimum Gibbs energy) among admissible roots"
    return best
