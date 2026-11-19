"""Free-text quantity extraction -> THERMO request inputs.

Pipeline:  text --extract_quantities--> numbers + Pint units
                --parse_text---------> (property, value, unit) tuples in the schema
                                       that ``request_parser`` / ``engine.validate`` expect.

This is a deliberately conservative, rule-based prototype, not NLP:

* A number only becomes an engine input if its *dimensions* identify exactly one
  property (K -> T, Pa -> P, m^3/kg -> v, kg/m^3 -> rho, J/(kg K) -> s).
* Anything ambiguous (two different temperatures, J/kg with no label, 2 kg, 1.5 L)
  is reported back in ``unassigned`` / ``conflicts`` -- never guessed.
* Everything returned still goes through ``engine.validate``; this module does not
  replace validation.
"""
import re
from functools import lru_cache

import pint

from .units import ureg, UNITS, BASE, PROPS, convert, UnitError
from .fluids import FLUIDS, list_fluids, get_fluid

# plain number: 400   -3.2   2e6   0.025
NUMBER_RE = re.compile(r"^[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?$")

# number with its unit glued on, no space: 350K   2MPa   25°C
NUMBER_GLUED_UNIT_RE = re.compile(
    r"^([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)([a-zA-Z°℃℉/^*²³]+[a-zA-Z°/^*²³0-9()-]*)$")

# Ordinary English words that are ALSO valid Pint unit symbols ("at" = atmosphere,
# "a" = are, "in" = inch ...). Without this list "... at 350K" swallows "at" into
# the unit. A blacklist is a stop-gap; a real parser would use grammar.
STOPWORDS = {
    "a", "at", "in", "is", "as", "of", "to", "was", "were", "and", "the",
    "with", "from", "by", "on", "for", "it", "that", "this", "we", "started",
    "then", "had", "took", "about", "are", "be", "or", "an", "per",
}

# In Pint "C" is coulomb and "F" is farad. In a thermodynamics sentence "25 C" is
# a temperature; silently reading it as coulombs would be a wrong-but-plausible
# failure, so these tokens are remapped explicitly.
UNIT_ALIASES = {"C": "degC", "°C": "degC", "℃": "degC", "F": "degF", "°F": "degF",
                "℉": "degF", "psia": "psi", "bara": "bar"}

# words that may sit between a label and its number: "enthalpy is 300 kJ/kg"
_FILLER = {"is", "was", "are", "of", "at", "=", ":", "equals", "="}
# label -> property, used only where dimensions alone cannot decide (J/kg is h, u, A or G)
LABELS = {
    "h": "h", "enthalpy": "h", "u": "u", "internal": "u", "energy": "u",
    "a": "A", "helmholtz": "A", "g": "G", "gibbs": "G", "s": "s", "entropy": "s",
    "t": "T", "temperature": "T", "p": "P", "pressure": "P",
    "v": "v", "rho": "rho", "density": "rho",
}

# common formulae -> catalogue names (case-sensitive; ambiguous ones like CO/Ar omitted)
FORMULAE = {"CO2": "Carbon dioxide", "CH4": "Methane", "H2O": "Water", "C2H6": "Ethane",
            "C3H8": "Propane", "NH3": "Ammonia", "H2S": "Hydrogen sulfide",
            "SO2": "Sulfur dioxide", "N2": "Nitrogen", "O2": "Oxygen"}


# ----------------------------------------------------------------------------
# unit helpers
# ----------------------------------------------------------------------------
@lru_cache(maxsize=2048)
def _resolve_unit(text):
    """Return a Pint-parsable unit string for ``text`` or None."""
    text = UNIT_ALIASES.get(text, text)
    if not text or NUMBER_RE.match(text):
        return None
    try:
        parsed = ureg.parse_units(text)
    except Exception:       # pint raises several unrelated error types here
        return None
    return None if str(parsed) == "dimensionless" else text


def _is_valid_unit(text):
    """True if Pint can parse `text` as a (non-dimensionless) unit."""
    return _resolve_unit(text) is not None


def _balanced(text):
    """Drop a trailing ')' that closes a '(' from outside the candidate, e.g. '(300 K)'."""
    while text.endswith(")") and text.count(")") > text.count("("):
        text = text[:-1]
    return text


def _strip_leading_paren(word):
    return re.sub(r"^\(+(?=[-+\d.])", "", word)


# ----------------------------------------------------------------------------
# extraction
# ----------------------------------------------------------------------------
def extract_quantities(sentence, max_unit_words=2):
    """Find every number in `sentence` and the unit (if any) attached to it.

    Returns a list of dicts: text, number, unit (string or None), quantity
    (a Pint Quantity in the project's registry), has_unit, label (the word
    before the number, lower-cased, or None).
    """
    # Remove sentence punctuation WITHOUT touching decimal points, and split
    # "T=300K" into "T", "300K". Only , . ; : not followed by a digit are removed.
    cleaned = re.sub(r"[,.;:](?!\d)", " ", sentence).replace("=", " = ")
    words = cleaned.split()

    found, i = [], 0
    while i < len(words):
        word = _strip_leading_paren(words[i])
        number_str = word if NUMBER_RE.match(word) else None
        glued_unit = None
        if number_str is None:
            m = NUMBER_GLUED_UNIT_RE.match(_balanced(word))
            if m:
                number_str, glued_unit = m.group(1), m.group(2)
        if number_str is None:
            i += 1
            continue

        value = float(number_str)
        unit_text, advance = None, 1
        if glued_unit:
            unit_text = _resolve_unit(glued_unit)
        else:
            # try the next 2 words, then 1, as the unit; stop at stopwords/numbers
            # because those mean we've walked past the unit into the next clause
            for n in range(max_unit_words, 0, -1):
                chunk = words[i + 1: i + 1 + n]
                if len(chunk) < n:
                    continue
                if any(w.lower() in STOPWORDS or NUMBER_RE.match(_strip_leading_paren(w)) for w in chunk):
                    continue
                candidate = _balanced(" ".join(chunk))
                if _is_valid_unit(candidate):
                    unit_text, advance = candidate, 1 + n
                    break

        label = None
        j = i - 1
        while j >= 0 and words[j].lower() in _FILLER:
            j -= 1
        if j >= 0:
            label = words[j].lower().strip("()")
            if label.endswith("+") or label.endswith("-"):
                label = None

        if unit_text:
            unit_text = UNIT_ALIASES.get(unit_text, unit_text)
            qty = ureg.Quantity(value, unit_text)
            found.append({"text": f"{number_str} {unit_text}", "number": value, "unit": unit_text,
                          "quantity": qty, "has_unit": True, "label": label})
        else:
            # bare number: still reported (flagged), never silently dropped
            found.append({"text": number_str, "number": value, "unit": None,
                          "quantity": ureg.Quantity(value, "dimensionless"),
                          "has_unit": False, "label": label})
        i += advance
    return found


# ----------------------------------------------------------------------------
# fluid detection
# ----------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _fluid_patterns():
    pats = []
    for name in list_fluids():
        spaced = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", name)              # CarbonDioxide -> Carbon Dioxide
        parts = [re.escape(p) for p in re.split(r"[\s_\-]+", spaced) if p]
        body = r"[\s_-]?".join(parts)
        flags = 0 if len(name) <= 4 else re.IGNORECASE                  # 'MM', 'D4', 'R11': exact case only
        pats.append((name, re.compile(r"(?<![A-Za-z0-9])" + body + r"(?![A-Za-z0-9])", flags)))
    return pats


def detect_fluids(text):
    """Distinct catalogue fluids mentioned in `text`. Longest match wins at a given spot;
    among duplicates ('Carbon dioxide' vs 'CarbonDioxide') the bundled-data entry wins."""
    hits = []                                   # (start, end, name)
    for name, pat in _fluid_patterns():
        for m in pat.finditer(text):
            hits.append((m.start(), m.end(), name))
    for formula, name in FORMULAE.items():
        if name in list_fluids():
            for m in re.finditer(r"(?<![A-Za-z0-9])" + formula + r"(?![A-Za-z0-9])", text):
                hits.append((m.start(), m.end(), name))
    hits.sort(key=lambda h: (-(h[1] - h[0]), h[2] not in FLUIDS))
    taken, chosen = [], []
    for s, e, name in hits:
        if any(s < te and ts < e for ts, te in taken):
            continue
        taken.append((s, e)); chosen.append((s, name))
    ordered = [n for _, n in sorted(chosen)]            # order of appearance in the text
    return list(dict.fromkeys(ordered))


# ----------------------------------------------------------------------------
# quantity -> engine property
# ----------------------------------------------------------------------------
def _dim(expr):
    return ureg.get_dimensionality(expr)


_DIM_TO_PROP = {_dim(BASE["temperature"]): "T", _dim(BASE["pressure"]): "P",
                _dim(BASE["volume"]): "v", _dim(BASE["density"]): "rho",
                _dim(BASE["entropy"]): "s"}
_ENERGY_DIM = _dim(BASE["energy"])
_MOLAR_VOLUME_DIM, _MOLAR_DENSITY_DIM = _dim("m^3/mol"), _dim("mol/m^3")


def _classify(item, molar_mass_g_mol):
    """-> (prop, value, unit, None) on success, or (None, None, None, reason)."""
    q = item["quantity"]
    if not item["has_unit"]:
        return None, None, None, "bare number, no unit recognised nearby"
    dim = q.dimensionality
    prop = _DIM_TO_PROP.get(dim)
    if prop is None and dim == _ENERGY_DIM:
        prop = LABELS.get(item["label"] or "")
        if prop not in ("h", "u", "A", "G"):
            return None, None, None, "J/kg could be h, u, A or G - put a label before it (e.g. 'enthalpy 300 kJ/kg')"
    if prop is not None:
        kind = PROPS[prop][1]
        for u in UNITS["SI"][kind]:                       # keep the user's own unit if the UI offers it
            if q.units == ureg.Unit(u):
                return prop, float(q.magnitude), u, None
        try:
            return prop, convert(q.magnitude, str(q.units), kind), BASE[kind], None
        except UnitError as e:
            return None, None, None, str(e)
    if dim in (_MOLAR_VOLUME_DIM, _MOLAR_DENSITY_DIM):
        if molar_mass_g_mol is None:
            return None, None, None, "molar quantity needs a fluid (molar mass) to convert to specific v / rho"
        M = molar_mass_g_mol / 1000.0                    # kg/mol
        if dim == _MOLAR_VOLUME_DIM:
            return "v", float(q.to("m^3/mol").magnitude) / M, BASE["volume"], None
        return "rho", float(q.to("mol/m^3").magnitude) * M, BASE["density"], None
    return None, None, None, f"dimensions {dim} do not match a THERMO input property (extensive quantities like mass or volume are not converted)"


# ----------------------------------------------------------------------------
# public entry point
# ----------------------------------------------------------------------------
def parse_text(text, fluid=None):
    """Turn a sentence into pieces of a THERMO request.

    Returns dict:
      fluid        catalogue name, or None (given `fluid` wins; ambiguity -> None + warning)
      inputs       [(prop, value, unit)] - only unambiguous, dimension-identified values
      unassigned   [{"text", "reason"}]  numbers that could not be used
      conflicts    [{"prop", "candidates"}]  same property stated with different values (all dropped)
      warnings     [str]
    """
    warnings = []
    fl = fluid
    if fl is None:
        names = detect_fluids(text)
        if len(names) == 1:
            fl = names[0]
        elif len(names) > 1:
            warnings.append("Several fluids mentioned (" + ", ".join(names) + "); choose one manually.")
    M = get_fluid(fl)["M"] if fl else None

    per_prop, unassigned = {}, []
    for item in extract_quantities(text):
        prop, val, unit, why = _classify(item, M)
        if prop is None:
            unassigned.append({"text": item["text"], "reason": why})
        else:
            per_prop.setdefault(prop, []).append((val, unit, item["text"]))

    inputs, conflicts = [], []
    for prop, cands in per_prop.items():
        # same state written twice (300 K ... 300 K) is not a conflict
        base = [convert(v, u, PROPS[prop][1]) for v, u, _ in cands]
        if all(abs(b - base[0]) <= 1e-9 * max(1.0, abs(base[0])) for b in base):
            inputs.append((prop, cands[0][0], cands[0][1]))
        else:
            conflicts.append({"prop": prop, "candidates": [c[2] for c in cands]})
    return dict(fluid=fl, inputs=inputs, unassigned=unassigned, conflicts=conflicts, warnings=warnings)


# ----------------------------------------------------------------------------
# demo
# ----------------------------------------------------------------------------
def pretty_print(sentence):
    """Demo helper: shows what was found, with SI base units so it's clear the
    Quantity is real and not just a string."""
    print(f'\n"{sentence}"')
    results = extract_quantities(sentence)
    if not results:
        print("  (no numbers found)")
        return
    for r in results:
        q = r["quantity"]
        if r["has_unit"]:
            print(f"  found: {r['text']:<15s} ->  {q:~P}   (= {q.to_base_units():~P})")
        else:
            print(f"  found: {r['text']:<15s} ->  bare number, no unit recognised nearby")
    p = parse_text(sentence)
    print(f"  request pieces: fluid={p['fluid']} inputs={p['inputs']}")
    for k in ("conflicts", "unassigned", "warnings"):
        if p[k]:
            print(f"  {k}: {p[k]}")


if __name__ == "__main__":
    for s in [
        "The temperature is 400 K and the pressure is 2 MPa.",
        "Assume a molar volume of 0.025 m^3/mol at 350K.",
        "Methane at 25 C and 10 bar.",
        "Cp0 was measured as 29.1 J/mol/K.",
        "We started at 298.15 K, 1 bar, then heated to 500 K.",
        "The sample had a mass of 2 kg and a volume of 1.5 liter.",
    ]:
        pretty_print(s)
