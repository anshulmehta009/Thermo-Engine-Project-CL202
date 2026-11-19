"""Three genuinely different computational routes for h and s. Reference: ideal-gas state ref=(T0,P0,h0,s0)."""
import numpy as np
from scipy import integrate
from . import eos, numerics as nx

def ideal_closed(p, T, P, ref):
    T0, P0, h0, s0 = ref
    if p.cp is not None:
        A, B, C, D = p.cp; R = p.Rs
        h = h0 + R * (A*(T-T0) + B/2*(T**2-T0**2) + C/3*(T**3-T0**3) - D*(1/T - 1/T0))
        s = s0 + R * (A*np.log(T/T0) + B*(T-T0) + C/2*(T**2-T0**2) - D/2*(1/T**2 - 1/T0**2)) - R*np.log(P/P0)
        return h, s
    from scipy.integrate import quad
    cp = p.cp_func
    hi = quad(cp, T0, T, epsabs=1e-7, epsrel=1e-8, limit=300)[0]
    si = quad(lambda t: cp(t)/t, T0, T, epsabs=1e-9, epsrel=1e-8, limit=300)[0]
    return h0 + hi, s0 + si - p.Rs*np.log(P/P0)

def feasibility(p, T, P, Z, has_cp, tgt_energy):
    if not tgt_energy: return {r: (False, "no energy-type target requested; routes not needed") for r in (1, 2, 3)}
    if not has_cp: return {r: (False, "ideal-gas cp parameters missing for this fluid") for r in (1, 2, 3)}
    out = {1: (p.u**2 - 4*p.w > 0, "closed-form residual properties from the EOS" if p.u**2 - 4*p.w > 0 else "closed form for u^2-4w = 0 (van der Waals) not implemented")}
    if T < p.Tc: out[2] = (False, "isothermal pressure path from P->0 crosses the saturation curve for T<Tc; two-phase path not implemented")
    elif abs(Z - eos.zmax(p, T, P)) > 1e-8: out[2] = (False, "selected root is not the vapour-side root tracked by the pressure path")
    else: out[2] = (True, "T>=Tc: single-branch pressure path available")
    out[3] = (True, "volume integrals of Maxwell-derived EOS derivatives (valid on any branch)")
    return out

def route1(p, st, ref, nm):
    T, P, Z, B = st["T"], st["P"], st["Z"], st["B"]; sy = eos.sym(p.model)
    aa, daa = [float(sy[k](*eos._args(p, T, 1.0))) for k in ("aa", "daa")]
    sq = np.sqrt(p.u**2 - 4*p.w); L = np.log((2*Z + B*(p.u+sq)) / (2*Z + B*(p.u-sq)))
    hres = p.Rs*T*(Z-1) + (T*daa - aa)/(p.b*sq)*L
    sres = p.Rs*np.log(Z-B) + daa/(p.b*sq)*L
    hi, si = ideal_closed(p, T, P, ref)
    return dict(h_ig=hi, s_ig=si, h_res=hres, s_res=sres, evals=0, error=None,
                equations=["h^R = RT(Z-1) + (T da/dT - a)/(b sqrt(u^2-4w)) * ln[(2Z+B(u+sqrt))/(2Z+B(u-sqrt))]",
                           "s^R = R ln(Z-B) + (da/dT)/(b sqrt(u^2-4w)) * ln[...]", "h = h_ig(T) + h^R ; s = s_ig(T,P) + s^R"])

def route2(p, st, ref, nm):
    T, P = st["T"], st["P"]; T0, P0, h0, s0 = ref
    cp = p.cp_func
    tol, n, m = nm["tol"], nm["simpson_n"], nm["integ"]
    hi, e1, n1 = nx.integrate_1d(cp, T0, T, m, tol, n); si, e2, n2 = nx.integrate_1d(lambda t: cp(t)/t, T0, T, m, tol, n)
    h_ig, s_ig = h0 + hi, s0 + si - p.Rs*np.log(P/P0)
    hstep = nm["fd_rel"] * T
    fd = (nx.fd_five_point if nm["diff"].startswith("Five-point") else
          nx.fd_central if nm["diff"].startswith("Central") else nx.fd_forward)
    dZdT = lambda Pp: fd(lambda t: eos.zmax(p, t, Pp), T, hstep)
    # x = ln P':  (v - T dv/dT) P' = -R T^2 dZ/dT ;  (R/P' - dv/dT) P' = -R (Z-1 + T dZ/dT)
    fh = lambda x: -p.Rs*T*T*dZdT(np.exp(x))
    fs = lambda x: -p.Rs*(eos.zmax(p, T, np.exp(x)) - 1 + T*dZdT(np.exp(x)))
    x0, x1 = np.log(P*1e-6), np.log(P)
    hr, e3, n3 = nx.integrate_1d(fh, x0, x1, m, tol, n); sr, e4, n4 = nx.integrate_1d(fs, x0, x1, m, tol, n)
    return dict(h_ig=h_ig, s_ig=s_ig, h_res=hr, s_res=sr, evals=n1+n2+n3+n4, error=dict(h=e1+e3, s=e2+e4),
                equations=["dH = Cp dT + [V - T(dV/dT)_P] dP ; dS = (Cp/T) dT - (dV/dT)_P dP",
                           "h_ig, s_ig: numerical integration of cp(T); (dV/dT)_P by finite differences of the vapour root Z(T,P')",
                           "pressure integral along the isotherm from P'=1e-6 P (ideal limit) to P"])

def route3(p, st, ref, nm):
    T, P, Z, v = st["T"], st["P"], st["Z"], st["v"]; sy = eos.sym(p.model); a = eos._args
    eu = lambda x: float(sy["eu"](*a(p, T, x))); es = lambda x: float(sy["es"](*a(p, T, x)))
    tol, n, m = nm["tol"], nm["simpson_n"], nm["integ"]
    if m == "Adaptive quadrature":
        ur, e1, i1 = integrate.quad(eu, v, np.inf, epsabs=1e-12, epsrel=tol, limit=500, full_output=1)[:3]
        sr, e2, i2 = integrate.quad(es, v, np.inf, epsabs=1e-12, epsrel=tol, limit=500, full_output=1)[:3]
        ev = i1["neval"] + i2["neval"]
    else:   # Simpson on t=v/v' in (0,1], truncated at t=1e-6 (tail ~1e-6 relative)
        g = lambda f: (lambda t: f(v/t) * v / t**2)
        ur, e1, n1 = nx.integrate_1d(g(eu), 1e-6, 1.0, m, tol, n); sr, e2, n2 = nx.integrate_1d(g(es), 1e-6, 1.0, m, tol, n); ev = n1 + n2
    hres = ur + p.Rs*T*(Z-1); sres = sr + p.Rs*np.log(Z)
    hi, si = ideal_closed(p, T, P, ref)
    return dict(h_ig=hi, s_ig=si, h_res=hres, s_res=sres, evals=ev, error=dict(h=e1, s=e2),
                equations=["Maxwell: (ds/dv)_T = (dP/dT)_v ; (du/dv)_T = T(dP/dT)_v - P  [(dP/dT)_v = SymPy derivative of the EOS]",
                           "u^R(T,v) = int_v^inf [P - T(dP/dT)_v] dv ; s^R(T,v) = int_v^inf [R/v - (dP/dT)_v] dv",
                           "h^R = u^R + RT(Z-1) ; s^R(T,P) = s^R(T,v) + R ln Z"])

ROUTES = {1: ("Direct EOS / property", route1), 2: ("Cp integration + pressure correction", route2), 3: ("Maxwell relations + EOS", route3)}

_L = r"\ln\frac{2Z+B(u+\sqrt{u^{2}-4w})}{2Z+B(u-\sqrt{u^{2}-4w})}"
EQUATIONS = {
 1: [r"h^{R}=RT(Z-1)+\frac{T\,\frac{d(a\alpha)}{dT}-a\alpha}{b\sqrt{u^{2}-4w}}" + _L, r"s^{R}=R\ln(Z-B)+\frac{d(a\alpha)/dT}{b\sqrt{u^{2}-4w}}" + _L, r"h=h^{ig}(T)+h^{R},\quad s=s^{ig}(T,P)+s^{R}"],
 2: [r"dH=C_P\,dT+\left[V-T\left(\frac{\partial V}{\partial T}\right)_P\right]dP", r"dS=\frac{C_P}{T}\,dT-\left(\frac{\partial V}{\partial T}\right)_P dP",
     r"h^{R}=-RT^{2}\int_{0}^{P}\left(\frac{\partial Z}{\partial T}\right)_{P'}\frac{dP'}{P'}", r"s^{R}=-R\int_{0}^{P}\left[(Z-1)+T\left(\frac{\partial Z}{\partial T}\right)_{P'}\right]\frac{dP'}{P'}",
     r"h^{ig}=h_0+\int_{T_0}^{T}C_P^{ig}\,dT,\quad s^{ig}=s_0+\int_{T_0}^{T}\frac{C_P^{ig}}{T}dT-R\ln\frac{P}{P_0}"],
 3: [r"\left(\frac{\partial s}{\partial v}\right)_T=\left(\frac{\partial P}{\partial T}\right)_v", r"u^{R}=\int_{v}^{\infty}\left[P-T\left(\frac{\partial P}{\partial T}\right)_v\right]dv",
     r"s^{R}(T,v)=\int_{v}^{\infty}\left[\frac{R}{v'}-\left(\frac{\partial P}{\partial T}\right)_v\right]dv'", r"h^{R}=u^{R}+RT(Z-1),\quad s^{R}(T,P)=s^{R}(T,v)+R\ln Z"]}
DERIVED = [r"u=h-Pv,\quad A=u-Ts,\quad G=h-Ts,\quad \rho=1/v"]
