import pytest
from thermo import engine
from thermo.quantity_parser import parse_text, extract_quantities, detect_fluids


def test_basic_state():
    r = parse_text("Methane at 300 K and 1 MPa.")
    assert r["fluid"] == "Methane" and r["inputs"] == [("T", 300.0, "K"), ("P", 1.0, "MPa")]


def test_glued_and_equals_style():
    assert parse_text("T=350K, P=2e6 Pa")["inputs"] == [("T", 350.0, "K"), ("P", 2e6, "Pa")]


def test_bare_C_is_celsius_not_coulomb():
    assert parse_text("25 C and 10 bar")["inputs"][0] == ("T", 25.0, "degC")
    assert parse_text("77F, 1 atm")["inputs"][0] == ("T", 77.0, "degF")


def test_conflicting_temperatures_are_not_guessed():
    r = parse_text("We started at 298.15 K, 1 bar, then heated to 500 K.")
    assert [c["prop"] for c in r["conflicts"]] == ["T"] and [i[0] for i in r["inputs"]] == ["P"]


def test_repeated_same_value_is_not_a_conflict():
    assert parse_text("300 K ... again 300 K and 1 bar")["conflicts"] == []


def test_extensive_quantities_not_treated_as_specific():
    r = parse_text("mass of 2 kg and a volume of 1.5 liter")
    assert r["inputs"] == [] and len(r["unassigned"]) == 2


def test_molar_volume_needs_fluid_then_converts():
    assert parse_text("0.025 m^3/mol at 350K")["inputs"] == [("T", 350.0, "K")]
    r = parse_text("CO2 with molar volume 0.025 m^3/mol at 350 K")
    v = dict((k, (val, u)) for k, val, u in r["inputs"])["v"]
    assert v[1] == "m^3/kg" and abs(v[0] - 0.025 / 0.04401) < 1e-3


def test_energy_needs_label():
    assert parse_text("300 kJ/kg")["inputs"] == []
    assert parse_text("enthalpy 300 kJ/kg")["inputs"][0][0] == "h"
    assert parse_text("1.5 kJ/(kg K)")["inputs"][0][0] == "s"


def test_stopwords_not_swallowed_into_units():
    assert [q["unit"] for q in extract_quantities("at 350K in 2 bar")] == ["K", "bar"]


def test_fluid_detection():
    assert detect_fluids("carbon dioxide at 350 K") == ["Carbon dioxide"]
    assert detect_fluids("CO2 and methane") == ["Carbon dioxide", "Methane"]
    assert detect_fluids("a 5 mm gap") == []          # 'MM' is a fluid; must not match 'mm'


def test_output_feeds_engine_validate():
    r = parse_text("Methane at 300 K and 1 MPa")
    req = dict(fluid=r["fluid"], model="PR", system="SI", inputs=r["inputs"], targets=["v"], phase="Any Phase")
    assert engine.validate(req)["ok"]
