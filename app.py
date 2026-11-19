"""THERMO – Streamlit frontend. No thermodynamic calculations here; everything goes through thermo.engine."""
import streamlit as st, pandas as pd, altair as alt
from thermo import engine, units, eos, routes, quantity_parser
from thermo.fluids import FLUIDS, list_fluids, get_fluid
from thermo.errors import ThermoError

st.set_page_config(page_title="THERMO", page_icon="🜂", layout="wide")
st.markdown("""<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
html,body,[class*="css"]{font-family:'Inter',sans-serif}
.block-container{max-width:1150px;padding-top:1.2rem}
.hero{background:linear-gradient(120deg,#F0ECFF 0%,#F7F5FF 60%,#fff 100%);border:1px solid #E8E4F2;border-radius:18px;padding:22px 28px;margin-bottom:14px}
.hero h1{margin:0;color:#6554C0;letter-spacing:.18em;font-size:2rem} .hero p{margin:4px 0 0;color:#77758A}
.steps{display:flex;gap:8px;margin:6px 0 16px} .step{flex:1;text-align:center;padding:8px;border-radius:12px;border:1px solid #E8E4F2;color:#77758A;font-size:.85rem;background:#fff}
.step.on{background:#6554C0;color:#fff;border-color:#6554C0;font-weight:600} .step.done{background:#F0ECFF;color:#6554C0}
div[data-testid="stVerticalBlockBorderWrapper"]{border-radius:16px!important;border-color:#E8E4F2!important;background:#fff}
h2,h3,h4{color:#25243A} .sec{color:#6554C0;font-weight:600;letter-spacing:.06em;font-size:.8rem;text-transform:uppercase;margin-bottom:4px}
.card{background:#F7F5FF;border:1px solid #E8E4F2;border-radius:14px;padding:14px 16px;margin-bottom:10px}
.cl{color:#77758A;font-size:.78rem;text-transform:uppercase;letter-spacing:.05em} .cv{color:#25243A;font-size:1.55rem;font-weight:700} .cu{color:#8B7CE8;font-size:.85rem}
.pill{display:inline-block;padding:2px 10px;border-radius:99px;background:#F0ECFF;color:#6554C0;font-size:.78rem;margin-right:6px}
.ok{background:#E7F6EC;color:#1E7A3C}.warn{background:#FFF4DB;color:#8A5A00}
div.stButton>button,div.stDownloadButton>button{background:#8B7CE8;color:#fff;border:0;border-radius:12px;padding:.55rem 1.5rem;font-weight:600}
div.stButton>button:hover,div.stDownloadButton>button:hover{background:#6554C0;color:#fff}
button[data-baseweb="tab"]{font-weight:600} .stTabs [aria-selected="true"]{color:#6554C0}
div[data-testid="stExpander"]{border-radius:12px;border-color:#E8E4F2;background:#FBFAFF}
</style>""", unsafe_allow_html=True)

S = st.session_state; S.setdefault("page", 1)
KEYS = ["T", "P", "v", "rho", "u", "h", "s", "A", "G"]
name = lambda k: f"{units.PROPS[k][0]} ({k})"
fmt = lambda x: f"{x:.6g}"
card = lambda l, v, u="": f"<div class='card'><div class='cl'>{l}</div><div class='cv'>{v}</div><div class='cu'>{u}</div></div>"
def go(n): S.page = n; st.rerun()

st.markdown("<div class='hero'><h1>THERMO</h1><p><b>Computational Engine for the Calculus of Thermodynamics</b> — an intelligent thermodynamic calculation and analysis environment.</p></div>", unsafe_allow_html=True)
labs = ["1 · Define problem", "2 · Validation", "3 · Route & numerics", "4 · Results"]
st.markdown("<div class='steps'>" + "".join(f"<div class='step {'on' if i+1 == S.page else 'done' if i+1 < S.page else ''}'>{l}</div>" for i, l in enumerate(labs)) + "</div>", unsafe_allow_html=True)

PRESETS = {
 "1 · Methane 300 K, 1 MPa (success)": dict(fluid="Methane", model="PR", props=["T", "P"], val_T=300.0, val_P=1e6, targets=["v", "h", "s"], phase="Any Phase"),
 "2 · Only T given (validation fails)": dict(fluid="Methane", model="PR", props=["T"], val_T=300.0, targets=["v", "h", "s"], phase="Any Phase"),
 "3 · Propane 300 K, 0.98 MPa (multiple roots)": dict(fluid="Propane", model="PR", props=["T", "P"], val_T=300.0, val_P=9.8e5, targets=["v", "rho", "h", "s"], phase="Any Phase"),
 "4 · Supercritical CO₂ 350 K, 15 MPa": dict(fluid="Carbon dioxide", model="SRK", props=["T", "P"], val_T=350.0, val_P=1.5e7, targets=["rho", "h", "s", "G"], phase="Supercritical"),
 "5 · Water vapour 500 K, 1 MPa (vdW)": dict(fluid="Water", model="vdW", props=["T", "P"], val_T=500.0, val_P=1e6, targets=["v", "h", "s"], phase="Vapour"),
}
def load(n):
    for k, v in PRESETS[n].items(): S[k] = v
    S.system = "SI"; S["u_T_SI"] = "K"; S["u_P_SI"] = "Pa"
def fill_from_text():
    """Callback: sentence -> fluid / input properties / values. Never touches targets, EOS or phase."""
    r = quantity_parser.parse_text(S.get("nl_text", ""))
    S.nl_report = r
    if r["fluid"]: S.fluid = r["fluid"]
    if r["inputs"]:
        S.system = "SI"; S.props = [k for k, _, _ in r["inputs"]]
        for k, val, unit in r["inputs"]: S[f"val_{k}"] = val; S[f"u_{k}_SI"] = unit
for k, v in dict(fluid="Methane", model="PR", props=["T", "P"], targets=["v", "h", "s"], phase="Any Phase", system="SI").items(): S.setdefault(k, v)

if S.page == 1:
    with st.container(border=True):
        st.markdown("<div class='sec'>Quick start</div>", unsafe_allow_html=True)
        c = st.columns([3, 1]); pick = c[0].selectbox("Demo scenarios", list(PRESETS), label_visibility="collapsed")
        c[1].button("Load scenario", on_click=load, args=(pick,))
        c = st.columns([3, 1]); c[0].text_input("Describe the state in words", key="nl_text", placeholder="e.g. Methane at 300 K and 1 MPa", label_visibility="collapsed")
        c[1].button("Fill from text", on_click=fill_from_text)
        rep = S.get("nl_report")
        if rep:
            if rep["inputs"] or rep["fluid"]: st.caption("Filled: " + ", ".join(([rep["fluid"]] if rep["fluid"] else []) + [f"{k} = {v:.6g} {u}" for k, v, u in rep["inputs"]]) + ". Check them below; targets, EOS and phase are not set from text.")
            for w in rep["warnings"]: st.warning(w)
            for x in rep["conflicts"]: st.warning(f"{x['prop']} stated more than once with different values ({', '.join(x['candidates'])}); not filled - enter it manually.")
            if rep["unassigned"]: st.caption("Ignored: " + "; ".join(f"{x['text']} ({x['reason']})" for x in rep["unassigned"]))
            if not (rep["inputs"] or rep["fluid"]): st.info("Nothing usable found in that text.")
    with st.container(border=True):
        st.markdown("<div class='sec'>A · Fluid</div>", unsafe_allow_html=True)
        available_fluids = list_fluids()
        selected_fluid = st.selectbox(
            "Search / select fluid",
            options=available_fluids,
            key="fluid",
        )
        f = get_fluid(selected_fluid)
        st.markdown(f"<span class='pill'>Pure fluid</span><span class='pill'>Tc = {f['Tc']} K</span><span class='pill'>Pc = {f['Pc']/1e6:.4g} MPa</span><span class='pill'>ω = {f['w']}</span><span class='pill'>M = {f['M']} g/mol</span>", unsafe_allow_html=True)
    with st.container(border=True):
        st.markdown("<div class='sec'>B · Unit system</div>", unsafe_allow_html=True)
        st.radio("System", list(units.UNITS), key="system", horizontal=True, label_visibility="collapsed")
        st.caption("Unit choices below follow this system. The backend converts everything to SI base units (K, Pa, m³/kg, J/kg, J/(kg·K)) with Pint before any calculation.")
    with st.container(border=True):
        st.markdown("<div class='sec'>C · Input properties</div>", unsafe_allow_html=True)
        st.multiselect("Add / remove properties", KEYS, key="props", format_func=name)
        st.caption("Provide any available independent properties; the engine evaluates whether they are sufficient.")
        rows = []
        for k in S.props:
            a, b = st.columns(2); kind = units.PROPS[k][1]
            val = a.number_input(name(k), value=None, key=f"val_{k}", format="%.6g", placeholder="enter value")
            rows.append((k, val, b.selectbox("Unit", units.UNITS[S.system][kind], key=f"u_{k}_{S.system}")))
    with st.container(border=True):
        st.markdown("<div class='sec'>D · Targets & phase</div>", unsafe_allow_html=True)
        st.multiselect("Target properties", KEYS, key="targets", format_func=name)
        st.radio("Phase preference", ["Liquid", "Vapour", "Supercritical", "Any Phase"], key="phase", horizontal=True)
    with st.container(border=True):
        st.markdown("<div class='sec'>E · Reference state</div>", unsafe_allow_html=True)
        mode = st.radio("Convention", ["Engine-defined (ideal gas, 298.15 K, 1 bar: h = s = 0)", "Specified reference"], key="refmode")
        S.ref = engine.ENGINE_REF
        if mode.startswith("Spec"):
            r = st.columns(4); S.ref = (r[0].number_input("T_ref [K]", value=298.15), r[1].number_input("P_ref [Pa]", value=1e5), r[2].number_input("h_ref [J/kg]", value=0.0), r[3].number_input("s_ref [J/(kg·K)]", value=0.0))
        st.caption("Absolute h, s, u, A, G depend on the reference convention; differences are the physically meaningful quantities. The convention is reported with results.")
    with st.container(border=True):
        st.markdown("<div class='sec'>F · Calculation model (EOS) — your choice, never auto-switched</div>", unsafe_allow_html=True)
        st.radio("EOS", list(eos.MODELS), key="model", horizontal=True, format_func=lambda m: eos.MODELS[m]["name"], label_visibility="collapsed")
        st.caption(eos.MODELS[S.model]["info"]); st.latex(eos.EOS_LATEX[S.model][0])
    if st.button("VALIDATE INPUTS"):
        S.req = dict(fluid=S.fluid, model=S.model, system=S.system, inputs=rows, targets=S.targets, phase=S.phase)
        S.val = engine.validate(S.req); go(2)

elif S.page == 2:
    with st.container(border=True):
        st.markdown("<div class='sec'>Validation report</div>", unsafe_allow_html=True)
        for i in S.val["issues"]:
            (st.error if i["severity"] == "error" else st.warning)(f"**{i['category']}** — `{i['field']}`: {i['message']}")
        if S.val["ok"]:
            st.success("All validation checks passed." if not S.val["issues"] else "Validation passed with warnings.")
            st.markdown("Converted inputs (SI):"); st.json(S.val["si"], expanded=False)
    if S.val["ok"]:
        if st.button("CONTINUE TO CALCULATION CONFIGURATION"): go(3)
    elif st.button("BACK TO INPUTS"): go(1)

elif S.page == 3:
    nm = dict(engine.DEFAULT_NM)
    try: stt = engine.prepare(S.req, S.val["si"], nm)
    except ThermoError as e:
        st.error(f"**{e.category}**: {e.message}")
        if st.button("Back to inputs"): go(1)
        st.stop()
    p = stt["p"]
    with st.container(border=True):
        st.markdown("<div class='sec'>Validated summary & EOS parameters</div>", unsafe_allow_html=True)
        c = st.columns(4)
        for col, (l, v, u) in zip(c, [("State", stt["label"], f"Z = {stt['Z']:.5f}"), ("T", fmt(stt["T"]), "K"), ("P", fmt(stt["P"]), "Pa"), ("v", fmt(stt["v"]), "m³/kg")]): col.markdown(card(l, v, u), unsafe_allow_html=True)
        st.dataframe(pd.DataFrame([dict(model=eos.MODELS[S.model]["name"], Tc=p.Tc, Pc=p.Pc, omega=p.omega, a=p.a, b=p.b, kappa=p.k)]), hide_index=True)
    need = bool(set(S.req["targets"]) & engine.ENERGY); route = None
    feas = routes.feasibility(p, stt["T"], stt["P"], stt["Z"], p.cp_func is not None, need)
    if need:
        with st.container(border=True):
            st.markdown("<div class='sec'>Computational routes</div>", unsafe_allow_html=True)
            dflt = stt.get("plan", {}).get("recommended_route") or next((r for r in (1, 3, 2) if feas[r][0]), 1); cols = st.columns(3)
            for r, col in zip((1, 2, 3), cols):
                ok = feas[r][0]; col.markdown(f"<div class='card'><span class='pill {'ok' if ok else 'warn'}'>{'feasible' if ok else 'unsupported'}</span>{'<span class=pill>default</span>' if r == dflt else ''}<br><b>{r} · {routes.ROUTES[r][0]}</b><br><span style='color:#77758A;font-size:.85rem'>{feas[r][1]}</span></div>", unsafe_allow_html=True)
            st.caption(f"Default = route {dflt}: analytic/closed form where available (lowest cost, no quadrature error). Override allowed; an unsupported route errors out and is never silently substituted.")
            route = st.selectbox("Route", [1, 2, 3], index=dflt - 1, format_func=lambda r: f"{r} · {routes.ROUTES[r][0]}")
    with st.container(border=True):
        st.markdown("<div class='sec'>Numerical methods</div>", unsafe_allow_html=True)
        if "P" in S.val["si"]: nm["root"] = st.selectbox("Root finding", ["Brent's method", "Bisection", "Newton-Raphson", "Secant method"])
        if route == 2:
            nm["diff"] = st.selectbox("Differentiation", ["Central finite difference", "Forward finite difference", "Five-point central difference"]); nm["fd_rel"] = st.number_input("FD relative step", value=1e-4, format="%.1e")
        if route in (2, 3):
            nm["integ"] = st.selectbox("Integration", ["Adaptive quadrature", "Simpson's rule"]); nm["simpson_n"] = int(st.number_input("Simpson intervals", 20, 20000, 200, 2))
        c = st.columns(3)
        nm["tol"] = c[0].number_input("Step / integration tolerance", value=1e-10, min_value=1e-14, format="%.1e")
        nm["residual_tol"] = c[1].number_input("Residual tolerance", value=1e-8, min_value=1e-14, format="%.1e")
        nm["maxit"] = int(c[2].number_input("Max iterations", 1, 1000, 100))
        S.cross = st.checkbox("Cross-check against other feasible routes (internal consistency only)")
    c = st.columns([1, 1, 4])
    if c[0].button("CALCULATE NOW"):
        try:
            with st.spinner("Calculating…"):
                stt = engine.prepare(S.req, S.val["si"], nm); S.res = engine.calculate(S.req, stt, route, nm, S.ref, S.cross); S.stt, S.err = stt, None
        except ThermoError as e: S.res, S.err = None, e
        go(4)
    if c[1].button("Back"): go(1)

else:
    if S.err:
        st.error(f"**{S.err.category}**: {S.err.message}")
    else:
        r, stt = S.res, S.stt; V = r["values"]
        st.success(f"Calculation completed · {eos.MODELS[S.model]['name']} (own implementation) · {r['info'].get('route_name', 'EOS root solution only')} · state: {r['label']}")
        t1, t2, t3, t4 = st.tabs(["Results", "Phase & roots", "Reference & limits", "Calculation details"])
        with t1:
            sysout = st.radio("Output units", list(units.UNITS), index=list(units.UNITS).index(S.system), horizontal=True)
            tab = []; cols = st.columns(3)
            for i, k in enumerate(S.req["targets"]):
                kind = units.PROPS[k][1]; un = units.UNITS[sysout][kind][0]; val = units.convert(V[k], un, kind, False)
                tab.append(dict(property=name(k), value=val, unit=un)); cols[i % 3].markdown(card(units.PROPS[k][0], fmt(val), un), unsafe_allow_html=True)
            st.dataframe(pd.DataFrame(tab), hide_index=True); st.download_button("Download CSV results", pd.DataFrame(tab).to_csv(index=False), "results.csv")
        with t2:
            if stt["cands"]: st.dataframe(pd.DataFrame([{k: c.get(k) for k in ("Z", "status", "reason", "label", "iters", "evals", "residual")} for c in stt["cands"]]), hide_index=True)
            else: st.info("T and v were given, so no root solving or phase elimination was required.")
        with t3:
            st.write(dict(T0=r["ref"][0], P0=r["ref"][1], h0=r["ref"][2], s0=r["ref"][3]))
            st.info("MODEL LIMITATIONS: results come from the selected EOS and implemented numerical methods; they are not experimental measurements. CoolProp is an independent reference backend and does not replace the selected EOS.")
            if st.checkbox("Show CoolProp reference properties (separate from our EOS)"):
                st.json(engine.coolprop_reference(S.fluid, stt["T"], stt["P"]), expanded=False)
            if st.checkbox("Show CoolProp fluid constants / EOS metadata"):
                meta = engine.fluid_database(S.fluid)
                st.json(meta, expanded=False)
        with t4:
            for t, c in r["trace"]:
                with st.expander(t):
                    if isinstance(c, dict) and "latex" in c:
                        for e in c["latex"]: st.latex(e)
                    elif isinstance(c, list) and c and isinstance(c[0], dict): st.dataframe(pd.DataFrame(c))
                    else: st.write(c)
    a, b, c = st.columns(3)
    if a.button("Edit inputs"): go(1)
    if b.button("Back to configuration"): go(3)
    if c.button("Reset all"): S.clear(); st.rerun()
