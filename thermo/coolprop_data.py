"""CoolProp-backed pure-fluid catalogue and reference-property access.

CoolProp is an independent property backend/reference source. Its fluid JSON
metadata is exposed as-is where available; this project does not claim that
CoolProp's internal Helmholtz equations are the same as the cubic EOS models
implemented in thermo.eos.
"""
from __future__ import annotations
from functools import lru_cache
import json
import re

try:
    import CoolProp.CoolProp as CP
    COOLPROP_AVAILABLE = True
except Exception:  # Optional at import time; a clear error is raised on use.
    CP = None
    COOLPROP_AVAILABLE = False


def require_coolprop():
    if not COOLPROP_AVAILABLE:
        raise RuntimeError("CoolProp is not installed. Install it with: pip install CoolProp")


@lru_cache(maxsize=1)
def fluid_names():
    """Return all pure-fluid names exposed by the installed CoolProp build.

    CoolProp builds/versions may return the catalogue separated by commas,
    ampersands, or semicolons. Split on all supported separators so the
    complete catalogue is never treated as one fluid name.
    """
    require_coolprop()
    raw = CP.get_global_param_string("fluids_list")
    if not isinstance(raw, str):
        raise TypeError(
            "CoolProp returned an unexpected fluids_list value: "
            f"{type(raw).__name__}"
        )

    names = {
        item.strip()
        for item in re.split(r"[&,;]", raw)
        if item.strip()
    }
    if not names:
        raise RuntimeError("CoolProp returned an empty fluid catalogue.")
    return sorted(names)


@lru_cache(maxsize=256)
def fluid_data(name: str):
    """Retrieve critical constants and available model metadata for one fluid."""
    require_coolprop()

    if not isinstance(name, str):
        raise TypeError(
            f"Expected one fluid name as a string, got {type(name).__name__}."
        )
    name = name.strip()
    if not name:
        raise ValueError("A fluid name must be selected before loading properties.")
    if any(separator in name for separator in (",", "&", ";")):
        raise ValueError(
            "Expected one fluid name, but received what looks like a "
            "whole fluid catalogue. Check fluid_names() separator handling."
        )

    aliases = [name, name.replace(" ", "")]
    # Resolve the human-readable labels used by this project's original UI.
    try:
        from .fluids import FLUIDS
        if name in FLUIDS:
            aliases.insert(0, FLUIDS[name].get("cool", name))
    except Exception:
        pass
    aliases = tuple(dict.fromkeys(aliases))
    last = None
    for fluid in aliases:
        try:
            get = lambda key: float(CP.PropsSI(key, fluid))
            data = {
                "coolprop_name": fluid,
                "Tc": get("Tcrit"), "Pc": get("pcrit"),
                "rhoc": get("rhomass_critical"),
                "M_kg_per_mol": get("molar_mass"),
                "omega": get("acentric"),
                "Tmin": get("Tmin"), "Tmax": get("Tmax"),
                "ptriple": get("ptriple"),
            }
            # The fluid JSON often contains EOS coefficients, reducing values,
            # alpha functions, ancillaries and other constants. Keep it raw so
            # we do not misrepresent its schema or silently discard fields.
            try:
                raw_json = CP.get_fluid_param_string(fluid, "JSON")
                data["json"] = json.loads(raw_json)
            except Exception as exc:
                data["json_error"] = str(exc)
            return data
        except Exception as exc:
            last = exc
    raise ValueError(f"CoolProp could not load fluid '{name}': {last}")


def props(fluid: str, T: float, P: float):
    """Independent CoolProp state properties (SI); phase errors are propagated."""
    require_coolprop()
    d = fluid_data(fluid)
    f = d["coolprop_name"]
    out = {"source": "CoolProp", "fluid": f, "T_K": float(T), "P_Pa": float(P)}
    for label, key in (("rho_kg_m3", "Dmass"), ("v_m3_kg", "Dmass"),
                       ("h_J_kg", "Hmass"), ("s_J_kgK", "Smass"),
                       ("cp_J_kgK", "Cpmass"), ("cv_J_kgK", "Cvmass"),
                       ("Z", "Z")):
        try:
            value = float(CP.PropsSI(key, "T", T, "P", P, f))
            out[label] = 1.0 / value if label == "v_m3_kg" else value
        except Exception as exc:
            out[label + "_error"] = str(exc)
    return out


def ideal_cp_mass(fluid: str, T: float):
    """Ideal-gas mass Cp where supported by CoolProp's CP0MASS output."""
    require_coolprop()
    f = fluid_data(fluid)["coolprop_name"]
    return float(CP.PropsSI("CP0MASS", "T", float(T), "P", 1e5, f))


def export_database(path: str, include_raw_json: bool = True):
    """Export the complete installed CoolProp pure-fluid catalogue to JSON.

    The output is generated from the exact CoolProp version installed in the
    user's environment, so the file also records the CoolProp version string.
    Individual fluid failures are captured rather than aborting the export.
    """
    require_coolprop()
    import pathlib
    version = CP.get_global_param_string("version")
    records, failures = {}, {}
    for name in fluid_names():
        try:
            d = dict(fluid_data(name))
            if not include_raw_json:
                d.pop("json", None)
            records[name] = d
        except Exception as exc:
            failures[name] = str(exc)
    payload = {"source": "CoolProp", "coolprop_version": version,
               "fluid_count_available": len(fluid_names()),
               "fluid_count_exported": len(records), "fluids": records,
               "failures": failures}
    target = pathlib.Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return {"path": str(target), "exported": len(records), "failed": len(failures),
            "coolprop_version": version}
