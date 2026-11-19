import numpy as np, pytest
from thermo import numerics as nx, units, eos, engine, routes
from thermo.errors import ThermoError
f = lambda x: x * x - 2

def test_roots():
    for r in (nx.bisection(f, 0, 2, 1e-12), nx.brent(f, 0, 2, 1e-12), nx.newton(f, lambda x: 2 * x, 1.0)):
        assert abs(r["root"] - 2**0.5) < 1e-10 and r["converged"]
def test_derivatives():
    assert abs(nx.fd_central(np.sin, 1.0, 1e-5) - np.cos(1.0)) < 1e-9
    assert abs(nx.fd_forward(np.sin, 1.0, 1e-7) - np.cos(1.0)) < 1e-6
def test_integration():
    for m in ("Adaptive quadrature", "Simpson's rule"):
        assert abs(nx.integrate_1d(np.sin, 0, np.pi, m, 1e-10, 200)[0] - 2) < 1e-8
def test_units():
    assert abs(units.convert(25, "degC", "temperature") - 298.15) < 1e-9
    assert abs(units.convert(1, "bar", "pressure") - 1e5) < 1e-6
    assert abs(units.convert(1, "g/cm^3", "density") - 1000) < 1e-9
    with pytest.raises(units.UnitError): units.convert(1, "kg", "pressure")
@pytest.mark.parametrize("m,zc", [("PR", 0.3074), ("SRK", 1 / 3), ("RK", 1 / 3), ("vdW", 0.375)])
def test_critical_Zc(m, zc):   # known analytical benchmark for cubic EOS
    p = eos.make_params("Methane", m); c, A, B = eos.z_poly(p, p.Tc, p.Pc)
    assert abs(np.poly1d(c)(zc)) < 1e-4
def test_methane_vs_coolprop():   # independent reference; EOS model deviation tolerated (1%)
    pytest.importorskip("CoolProp")
    p = eos.make_params("Methane", "PR"); ref = engine.coolprop_reference("Methane", 300, 1e6)
    assert abs(eos.zmax(p, 300, 1e6) - ref["Z"]) / ref["Z"] < 0.01
def _req(**k):
    d = dict(fluid="Carbon dioxide", model="PR", system="SI", inputs=[("T", 400, "K"), ("P", 5, "MPa")], targets=["h", "s", "v"], phase="Any Phase")
    d.update(k); return d
def test_validation_categories():
    assert engine.validate(_req())["ok"]
    assert engine.validate(_req(inputs=[("T", 400, "K")]))["issues"][0]["category"] == "Insufficient independent properties"
    assert engine.validate(_req(inputs=[("T", -5, "K"), ("P", 1, "bar")]))["issues"][0]["category"] == "Physically inadmissible state"
    from thermo.fluids import FLUIDS
    saved = FLUIDS["Water"]["cp"]; FLUIDS["Water"]["cp"] = None
    try: assert any(i["category"] == "Missing parameters" for i in engine.validate(_req(fluid="Water"))["issues"])
    finally: FLUIDS["Water"]["cp"] = saved
    assert any(i["category"] == "Unsupported fluid-model combination" for i in engine.validate(_req(model="BWR"))["issues"])
def test_routes_consistent():   # internal consistency only, NOT independent verification
    req = _req(); nm = dict(engine.DEFAULT_NM); st = engine.prepare(req, engine.validate(req)["si"], nm)
    res = {r: engine.calculate(req, st, r, nm)["values"] for r in (1, 2, 3)}
    for r in (2, 3):
        assert abs(res[r]["h"] - res[1]["h"]) < 1e-3 * abs(res[1]["h"]) + 5.0
        assert abs(res[r]["s"] - res[1]["s"]) < 1e-3 * abs(res[1]["s"]) + 0.5
def test_route2_unsupported_below_Tc_no_substitution():
    req = _req(fluid="Methane", inputs=[("T", 150, "K"), ("P", 1, "bar")]); nm = dict(engine.DEFAULT_NM)
    st = engine.prepare(req, engine.validate(req)["si"], nm)
    assert not routes.feasibility(st["p"], 150.0, 1e5, st["Z"], True, True)[2][0]
    with pytest.raises(ThermoError): engine.calculate(req, st, 2, nm)

def test_knowledge_base_dependencies():
    from thermo.knowledge_base import dependencies, THERMODYNAMIC_RELATIONS
    d = dependencies(["h", "u", "A"])
    assert "Cp ideal-gas data" in d["h"]["needs"]
    assert d["u"]["relation"] == "u = h - P v"
    assert "maxwell_1" in THERMODYNAMIC_RELATIONS

def test_uncertainty_and_tolerance():
    from thermo.uncertainty import validate_tolerances, propagate_independent_uncertainties
    assert validate_tolerances(1e-8, 50)["max_iterations"] == 50
    with pytest.raises(ValueError): validate_tolerances(-1, 10)
    out = propagate_independent_uncertainties({"T": 2.0, "P": 3.0}, {"T": 0.1, "P": 0.2})
    assert abs(out["standard_uncertainty"] - (0.2**2 + 0.6**2)**0.5) < 1e-12

def test_five_point_derivative():
    assert abs(nx.fd_five_point(np.sin, 1.0, 1e-2) - np.cos(1.0)) < 1e-8

def test_secant_root():
    r = nx.secant(f, 1.0, 2.0, 1e-12, 100)
    assert r["converged"] and abs(r["root"] - 2**0.5) < 1e-10

def test_fluid_catalogue_fallback():
    from thermo.fluids import list_fluids
    assert "Methane" in list_fluids()

def test_structured_request_parser():
    from thermo.request_parser import parse_request
    req = _req()
    assert parse_request(req)["fluid"] == "Carbon dioxide"
    with pytest.raises(ThermoError): parse_request({"fluid": "Water"})


def test_coolprop_catalogue_fluid_can_build_eos_params(monkeypatch):
    """A fluid outside the bundled 22 must be usable when CoolProp supplies its constants."""
    from thermo import eos
    dynamic_record = {
        "Tc": 304.1282, "Pc": 7.3773e6, "w": 0.22394, "M": 44.01,
        "cp": None, "cool": "CarbonDioxide",
        "source": "CoolProp fluid database",
    }
    monkeypatch.setattr(eos, "get_fluid", lambda name: dynamic_record)
    p = eos.make_params("CoolProp test fluid", "PR")
    assert p.fluid == "CoolProp test fluid"
    assert p.data_source == "CoolProp fluid database"
    assert p.cp is None
    assert callable(p.cp_func)


def test_coolprop_cp_failure_is_reported_as_thermo_error(monkeypatch):
    from thermo import eos, coolprop_data
    dynamic_record = {
        "Tc": 304.1282, "Pc": 7.3773e6, "w": 0.22394, "M": 44.01,
        "cp": None, "cool": "CoolPropOnlyFluid",
        "source": "CoolProp fluid database",
    }
    monkeypatch.setattr(eos, "get_fluid", lambda name: dynamic_record)
    monkeypatch.setattr(
        coolprop_data, "ideal_cp_mass",
        lambda fluid, T: (_ for _ in ()).throw(ValueError("CP0MASS unavailable")),
    )
    p = eos.make_params("CoolProp test fluid", "PR")
    with pytest.raises(ThermoError, match="CoolProp could not evaluate ideal-gas Cp"):
        p.cp_func(300.0)
