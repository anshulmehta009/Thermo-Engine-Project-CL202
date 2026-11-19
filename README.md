# THERMO
## Computational Engine for the Calculus of Thermodynamics

THERMO is a Streamlit-based thermodynamic calculation application. It accepts a pure-fluid state specification, validates the inputs, converts units to an internal SI basis, evaluates a user-selected equation of state (EOS), determines an admissible compressibility-factor root, and calculates requested thermodynamic properties.

The project separates the user interface from the calculation backend. The Streamlit interface is in `app.py`; thermodynamic models and algorithms are implemented in the `thermo/` package.

> **Model boundary:** THERMO's selected EOS calculations are performed by this project's own implementations of van der Waals, Redlich–Kwong, Soave–Redlich–Kwong, Peng–Robinson (1976), and Peng–Robinson (1978). CoolProp is an optional, independent source of fluid data and reference properties. It does not silently replace the selected EOS.

---

# 1. Project folder structure

```text
thermo_engine/
│
├── app.py
│   └── Streamlit user interface and four-step workflow
│
├── requirements.txt
│   └── Python packages required to run the application
│
├── README.md
│   └── Main project documentation and execution instructions
│
├── thermo/
│   ├── __init__.py
│   ├── engine.py
│   ├── eos.py
│   ├── fluids.py
│   ├── coolprop_data.py
│   ├── routes.py
│   ├── planner.py
│   ├── numerics.py
│   ├── units.py
│   ├── knowledge_base.py
│   ├── symbolic.py
│   ├── uncertainty.py
│   ├── request_parser.py
│   ├── quantity_parser.py
│   └── errors.py
│
├── scripts/
│   └── export_coolprop_database.py
│
├── data/
│   └── README.md
│
├── docs/
│   └── PROJECT_EXECUTION_ARCHITECTURE.md
│
├── tests/
│   ├── test_all.py
│   └── test_quantity_parser.py
│
└── .streamlit/
    └── config.toml
```

## 1.1 Main execution sequence

```text
User enters fluid, EOS, state inputs and target properties
                         |
                         v
                 app.py (Streamlit UI)
                         |
                         v
              request_parser.py
          checks request structure and fields
                         |
                         v
                  units.py
          converts entered units into SI
                         |
                         v
                 engine.validate()
       checks physical and thermodynamic inputs
                         |
                         v
                  fluids.py
       loads fluid constants and Cp information
                         |
                         v
                   eos.py
       constructs EOS parameters and Z equation
                         |
                         v
           root solving and root filtering
                         |
                         v
             planner.py / routes.py
       checks route feasibility and recommends
                         |
                         v
                 engine.calculate()
       evaluates requested thermodynamic properties
                         |
                         v
              Results shown in app.py
```

The interface does not contain the EOS calculation equations. It collects the request and calls the backend functions in `thermo/engine.py`.

---

# 2. Description of every code file

## 2.1 `app.py` — Streamlit application and user interface

This is the main file used to launch THERMO. It provides the interactive interface and controls the sequence of operations.

The application is organised into four workflow stages:

1. **Define problem:** select a fluid, unit system, input properties, target properties, phase preference, reference-state convention and EOS.
2. **Validation:** display validation issues and show the converted SI inputs when validation succeeds.
3. **Route and numerical configuration:** display the calculated state and EOS parameters, show feasible computational routes, and allow the user to select numerical settings.
4. **Results:** display calculated properties, phase/root information, reference-state information and calculation details.

The interface includes five Quick Start presets:

| Preset | Fluid | EOS | State | Purpose |
|---|---|---|---|---|
| 1 | Methane | PR | 300 K, 1 MPa | Successful calculation example |
| 2 | Methane | PR | Temperature only | Demonstrates insufficient independent inputs |
| 3 | Propane | PR | 300 K, 0.98 MPa | Demonstrates multiple mathematical roots |
| 4 | Carbon dioxide | SRK | 350 K, 15 MPa | Supercritical-state example |
| 5 | Water | vdW | 500 K, 1 MPa | Vapour-state example |

The user can request properties such as specific volume, density, internal energy, enthalpy, entropy, Helmholtz energy and Gibbs energy.

The Results view provides a table of requested values and units, and a CSV download for the current calculation. The other result sections show retained/root-solver information, reference-state details, optional CoolProp reference data and the calculation trace.

The current interface does not display the residuals/convergence dashboard or the P–v diagram.

## 2.2 `thermo/fluids.py` — bundled pure-fluid database

This file contains the project's bundled pure-fluid records. Each record stores:

- Critical temperature, `Tc` (K)
- Critical pressure, `Pc` (Pa)
- Acentric factor, `w` or omega
- Molecular weight, `M` (g/mol)
- Ideal-gas heat-capacity polynomial coefficients
- Corresponding CoolProp fluid identifier
- Data-source label

The 22 bundled fluid records are:

1. Methane
2. Ethane
3. Propane
4. n-Butane
5. Isobutane
6. n-Pentane
7. n-Hexane
8. Ethylene
9. Propylene
10. Acetylene
11. Benzene
12. Carbon dioxide
13. Carbon monoxide
14. Nitrogen
15. Oxygen
16. Argon
17. Water
18. Ammonia
19. Hydrogen sulfide
20. Sulfur dioxide
21. Methanol
22. Ethanol

### Heat-capacity representation

For bundled records, the heat-capacity coefficients are used in the implementation through the form:

\[
C_p^{ig}(T)=R_s\left(A+BT+CT^2+\frac{D}{T^2}\right)
\]

where `R_s` is the specific gas constant and `A`, `B`, `C`, and `D` are the stored coefficients.

The bundled correlation range is documented in the code as approximately 250–1500 K. The application issues a warning when an energy-property calculation uses a bundled polynomial outside that range. This range should be treated as a documented applicability guide, not as a guarantee of experimental accuracy for every fluid.

### Main functions

- `list_fluids()` returns bundled names and, when CoolProp is available, also adds names exposed by the installed CoolProp catalogue.
- `get_fluid(name)` returns the bundled record for a bundled fluid. For a non-bundled name, it attempts to retrieve metadata through `coolprop_data.py`.

**CoolProp-fluid EOS support:** `eos.make_params()` resolves bundled fluids first and then uses `get_fluid()` to load non-bundled fluids from the installed CoolProp catalogue. The same project-owned vdW, RK, SRK, PR and PR78 equations are used for these fluids; CoolProp supplies fluid constants and (when needed) ideal-gas Cp data, not the selected EOS result. Fluids for which CoolProp cannot provide the required constants or ideal-gas Cp are reported as unsupported for the relevant calculation.

## 2.3 `thermo/coolprop_data.py` — CoolProp integration and independent reference data

This module connects THERMO to the external CoolProp property library.

CoolProp is used as an independent data/reference source. It is separate from the project's own cubic EOS implementation.

### Fluid catalogue

`fluid_names()` reads the installed CoolProp catalogue using:

```python
CP.get_global_param_string("fluids_list")
```

It splits the returned catalogue into individual names and caches the result. The separator handling supports ampersand, comma and semicolon delimiters.

### Fluid constants and metadata

`fluid_data(name)` requests the following information from CoolProp:

- Critical temperature (`Tcrit`)
- Critical pressure (`pcrit`)
- Critical mass density (`rhomass_critical`)
- Molar mass (`molar_mass`)
- Acentric factor (`acentric`)
- Minimum and maximum temperature limits
- Triple-point pressure

It also attempts to retrieve the fluid's raw JSON metadata through `get_fluid_param_string(fluid, "JSON")`. The JSON may contain model coefficients, reducing parameters, ancillary correlations and other CoolProp-specific metadata. It is preserved as returned rather than being treated as a directly interchangeable cubic-EOS parameter set.

### Reference state properties

`props(fluid, T, P)` requests independent CoolProp state properties at a specified temperature and pressure, including:

- Density
- Specific volume
- Enthalpy
- Entropy
- Cp
- Cv
- Compressibility factor

If an individual property cannot be evaluated, the function records an error for that property.

`ideal_cp_mass(fluid, T)` requests ideal-gas mass heat capacity through CoolProp's `CP0MASS` output. This is used by the fluid parameter construction path for a CoolProp-backed record when that path is supported.

### Database export

`export_database(path, include_raw_json=True)` loops through the installed CoolProp pure-fluid catalogue and writes a JSON snapshot containing the CoolProp version, available fluid count, successfully exported records and per-fluid failures.

CoolProp must be installed for these operations. The CoolProp version and catalogue available on the executing computer determine the exported content.

## 2.4 `thermo/eos.py` — equation-of-state models and state/root calculations

This is the core thermodynamic model module. It implements the selected EOS equations rather than delegating the selected model calculation to CoolProp.

The general cubic form used is:

\[
P=\frac{R_sT}{v-b}
-\frac{a\alpha(T)}{v^2+ubv+wb^2}
\]

The compressibility factor is:

\[
Z=\frac{Pv}{R_sT}
\]

The dimensionless parameters are:

\[
A=\frac{a\alpha P}{R_s^2T^2},
\qquad
B=\frac{bP}{R_sT}
\]

The corresponding cubic equation in `Z` is solved to identify possible state roots.

### Implemented EOS models

| Model | Model name | Main characteristics in this implementation |
|---|---|---|
| vdW | van der Waals | Original cubic EOS; constant alpha |
| RK | Redlich–Kwong | Temperature-dependent alpha, \(T_r^{-1/2}\) |
| SRK | Soave–Redlich–Kwong | Soave alpha function and acentric-factor correction |
| PR | Peng–Robinson (1976) | Soave alpha function; commonly used cubic EOS |
| PR78 | Peng–Robinson (1978) | PR model with an extended kappa correlation for high acentric factor |

The models use their critical-property-based constants and alpha-function definitions as coded in `MODELS` and `EOS_LATEX`.

### Root handling

The implementation:

1. Forms the cubic polynomial in `Z`.
2. Uses NumPy polynomial roots to identify mathematical roots.
3. Rejects roots with a significant imaginary component.
4. Refines real candidates using the selected root solver.
5. Checks convergence and the cubic residual tolerance.
6. Rejects roots for which `Z <= B`, corresponding to a volume at or below the EOS co-volume.
7. Rejects mechanically unstable roots where \((\partial P/\partial v)_T \geq 0\).
8. Labels retained roots by phase logic.
9. Applies the requested phase preference.
10. Among admissible roots, selects the root with the lowest logarithmic fugacity coefficient, corresponding to the minimum Gibbs-energy criterion at fixed temperature and pressure.

The phase/root table in the interface makes the root classification and selection trace inspectable.

## 2.5 `thermo/routes.py` — enthalpy and entropy calculation routes

This module contains three different methods for calculating residual enthalpy and entropy contributions. Each route combines ideal-gas contributions with residual contributions from the selected EOS.

### Route 1 — Direct EOS / property route

This route evaluates the closed-form residual enthalpy and entropy expressions for the selected cubic EOS.

It uses the EOS alpha function and its temperature derivative, together with the selected compressibility root.

**Characteristics:**
- Analytical/closed-form residual-property expressions.
- No numerical integration is required for the residual terms.
- The route requires a real non-zero value of the discriminant expression \(u^2-4w\) in this implementation.
- The van der Waals parameter combination has \(u^2-4w=0\), so this closed-form implementation does not mark Route 1 feasible for vdW.

### Route 2 — Cp integration + pressure correction

This route calculates ideal-gas enthalpy and entropy by numerical integration of Cp and Cp/T. It then calculates residual corrections through pressure-path integrals involving the temperature derivative of the compressibility factor.

The pressure-path implementation follows an isothermal path from a low-pressure approximation toward the specified pressure.

**Feasibility conditions:**
- Ideal-gas Cp data must be available.
- The state temperature must be at or above the fluid critical temperature.
- The selected root must match the vapour-side root tracked by the pressure path.

For a subcritical temperature, the pressure path may cross a saturation/two-phase region. Since a two-phase path is not implemented, Route 2 is marked infeasible in that case. It is also rejected when the selected root is not the vapour-side root.

### Route 3 — Maxwell relations + EOS

This route uses Maxwell relations and derivatives of the selected EOS to calculate residual internal-energy and entropy integrals over specific volume. The residual enthalpy and entropy are then assembled from these integrals and compressibility-factor terms.

The EOS derivatives are generated using the symbolic expressions in `eos.py`/`symbolic.py`; the integrals are evaluated numerically.

**Characteristics:**
- Based on thermodynamic Maxwell relations.
- Uses EOS derivatives with respect to temperature and volume.
- Uses adaptive quadrature or Simpson's rule, depending on the selected integration method.
- The route is marked feasible when Cp data are available; the current feasibility function does not impose the Route 2 vapour-path restriction on Route 3.

### How the route is selected

Route selection is feasibility-first. The planner checks the state and target properties, then marks each route as feasible or infeasible with a reason.

The planner uses this preference order among feasible routes:

1. Route 1 — direct analytical/closed-form calculation
2. Route 3 — Maxwell relations and EOS volume integrals
3. Route 2 — Cp integration and pressure correction

This order is encoded as cost ranks in `thermo/planner.py`:

- Route 1: cost rank 1
- Route 3: cost rank 2
- Route 2: cost rank 3

The recommended route is the feasible route with the lowest cost rank. This means the route is not selected simply because it has the lowest route number: it must first satisfy its feasibility conditions.

The interface displays the feasibility reason for each route. The user may choose another route, but if that route is infeasible, the calculation raises an error. **The engine does not silently switch to a different route.**

Routes are only needed when an energy-type target is requested: `h`, `s`, `u`, `A` or `G`. If the requested targets are only state properties such as `v`, `rho` or `Z`, the energy-property routes are not required.

## 2.6 `thermo/planner.py` — feasibility-based route planner

This module coordinates route selection.

It:
- Checks whether the requested targets require energy-property calculations.
- Requests route feasibility information from `routes.py`.
- Includes property dependency information from `knowledge_base.py`.
- Assigns a cost rank to each route.
- Recommends the lowest-cost feasible route.
- Returns a reason for every route's feasibility status.

The selection basis is explicitly recorded as: feasibility first, then analytical/closed-form methods before numerical-integral methods.

## 2.7 `thermo/numerics.py` — numerical algorithms

This module provides numerical algorithms used by the EOS and property routes.

The available methods include:

**Root finding**
- Bisection
- Brent's method
- Newton–Raphson
- Secant method

**Numerical differentiation**
- Central finite difference
- Forward finite difference
- Five-point central finite difference

**Numerical integration**
- Adaptive quadrature
- Simpson's rule

The user can configure the root method, differentiation method, integration method, step/integration tolerance, residual tolerance, maximum iterations, finite-difference relative step and Simpson interval count through the interface.

## 2.8 `thermo/units.py` — unit conversion

This module uses Pint to convert user-entered quantities into internal SI base units before calculation.

The internal base units are:

| Quantity | Internal unit |
|---|---|
| Temperature | K |
| Pressure | Pa |
| Specific volume | m³/kg |
| Density | kg/m³ |
| Energy per unit mass | J/kg |
| Entropy per unit mass | J/(kg·K) |

The interface offers SI, CGS, US Customary and Engineering (metric) unit systems. The conversion function supports conversion to the internal base unit and conversion from the base unit to the selected display unit.

Unit incompatibility is raised as a `UnitError` and reported during validation.

## 2.9 `thermo/knowledge_base.py` — thermodynamic relationships and dependencies

This module stores the thermodynamic relationships used to describe property dependencies and route requirements.

It includes:
- First-law differential relation
- Enthalpy relation
- Helmholtz and Gibbs energy definitions
- Maxwell relations
- Cp–Cv relation
- Differential enthalpy and entropy relations
- Target-property dependency information

For example:

\[
u=h-Pv,\qquad A=u-Ts,\qquad G=h-Ts
\]

The planner uses the dependency registry to explain which data or intermediate properties are required for requested targets.

## 2.10 `thermo/symbolic.py` — symbolic mathematics

This module uses SymPy for symbolic operations.

It provides:
- Differentiation
- Integration
- Simplification
- Symbolic temperature and volume derivatives of an EOS expression

The EOS module uses symbolic derivatives to construct expressions needed by the Maxwell-relation route.

## 2.11 `thermo/uncertainty.py` — tolerance checks and uncertainty helper

This module validates numerical settings such as positive finite tolerances and a valid maximum iteration count.

It also provides a first-order uncertainty propagation helper for independent input uncertainties. The method uses root-sum-square propagation:

\[
\sigma_y=\sqrt{\sum_i\left(\frac{\partial y}{\partial x_i}\sigma_{x_i}\right)^2}
\]

The helper requires derivative and standard-uncertainty values to be supplied. The application does not invent measurement uncertainties or automatically attach uncertainty estimates to every calculated property.

## 2.12 `thermo/request_parser.py` — request structure validation

This module checks that a calculation request contains the required fields and that its input entries have the expected structure.

It normalises fluid and model names, converts target collections into lists and checks that each input entry contains a property, value and unit.

It raises a structured `ThermoError` when the request format is invalid.

## 2.13 `thermo/errors.py` — structured error handling

This module defines `ThermoError`, which carries an error category, a message and an optional field name.

The categories include:
- Invalid input
- Insufficient independent properties
- Unsupported fluid-model combination
- Missing parameters
- Mathematical failure
- Physically inadmissible state
- Numerical non-convergence

This allows the application to display meaningful validation and calculation messages rather than exposing only an unclassified exception.

## 2.14 `thermo/__init__.py`

This file marks `thermo` as a Python package so its modules can be imported using package paths such as `from thermo import engine`.

---

# 3. Input requirements and validation

The current state-calculation workflow is designed for a single-phase pure-fluid state.

The engine requires:
- One selected fluid.
- One implemented EOS model.
- At least one requested target property.
- Temperature plus at least one of pressure, specific volume or density.

The state can therefore be specified using:
- `T` and `P`
- `T` and `v`
- `T` and `rho`

When density is supplied, the engine calculates specific volume using \(v=1/\rho\).

The current implementation does not provide general inverse state solving from enthalpy, entropy or other energy properties. It also checks consistency if both specific volume and density are supplied.

The validator checks finite numerical inputs, compatible units, positive absolute temperature/pressure/volume/density, EOS co-volume restrictions, pressure-volume consistency and selected phase conditions.

For a supercritical phase preference, both \(T>T_c\) and \(P>P_c\) are required by the current validator.

---

# 4. Thermodynamic properties

The interface can request the following properties:

| Symbol | Property | Meaning |
|---|---|---|
| T | Temperature | Thermodynamic temperature |
| P | Pressure | Thermodynamic pressure |
| v | Specific volume | Volume per unit mass |
| rho | Density | Mass per unit volume |
| u | Internal energy | Specific internal energy |
| h | Enthalpy | Specific enthalpy |
| s | Entropy | Specific entropy |
| A | Helmholtz energy | Specific Helmholtz free energy |
| G | Gibbs energy | Specific Gibbs free energy |

The engine calculates:

\[
\rho=\frac{1}{v}
\]

\[
u=h-Pv
\]

\[
A=u-Ts
\]

\[
G=h-Ts
\]

The energy and entropy values depend on the selected reference-state convention.

The default reference state is:

- \(T_0=298.15\ \mathrm{K}\)
- \(P_0=100000\ \mathrm{Pa}\) (1 bar)
- \(h_0=0\ \mathrm{J/kg}\)
- \(s_0=0\ \mathrm{J/(kg\,K)}\)

The user can alternatively specify a reference state through the interface.

---

# 5. Output and calculation trace

After a successful calculation, the application displays the requested target properties in the selected output unit system.

The Results section includes:
- Property name
- Calculated value
- Display unit
- Downloadable CSV for the current calculation

The Phase & roots section shows the mathematical root candidates and their status, reason, phase label and available solver information. If the state was specified using temperature and specific volume, pressure is obtained from the EOS and no compressibility-root solving is required.

The Reference & limits section reports the reference-state values and model limitations. It also provides optional independent CoolProp reference properties and fluid metadata when CoolProp is installed and the selected fluid is supported by that data path.

The Calculation details section exposes the calculation trace, including converted SI inputs, EOS parameters and equations, fluid-data provenance, root filtering, property dependencies, route analysis, selected route and numerical settings.

The CSV download in the interface is for the current calculation. It is not an automatic batch export of every possible fluid, EOS and state combination.

---

# 6. Installation and execution instructions

## 6.1 Requirements

Use Python 3.10 or later where possible. The project dependencies are listed in `requirements.txt`:

- Streamlit
- Altair
- NumPy
- SciPy
- SymPy
- CoolProp
- Pint
- pandas
- pytest

## 6.2 Extract the project

Extract the ZIP file. Open a terminal or Command Prompt in the directory that contains `app.py` and the `thermo/` folder.

## 6.3 Create and activate a virtual environment (recommended)

Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
```

macOS/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

## 6.4 Install dependencies

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

CoolProp is optional when using only the 22 bundled fluids. Install CoolProp to enable its extended fluid catalogue, CoolProp-backed fluid constants and ideal-gas Cp, independent reference properties and database export. Pint is required for unit conversion.

## 6.5 Run the application

From the project root:

```bash
streamlit run app.py
```

Streamlit will print a local URL, normally:

```text
http://localhost:8501
```

Open that URL in a browser.

## 6.6 Basic operation

1. Choose a Quick Start scenario or select a fluid manually.
2. Select the unit system.
3. Add the independent input properties and enter their values and units.
4. Select the target properties.
5. Choose the phase preference.
6. Select the reference-state convention.
7. Select the EOS model.
8. Click **VALIDATE INPUTS**.
9. Review the validation report.
10. Continue to calculation configuration if validation succeeds.
11. Review route feasibility and the recommended route.
12. Select the root solver, differentiation method, integration method and numerical tolerances as required.
13. Click **CALCULATE NOW**.
14. Review the result sections and download the current result CSV if needed.

---

# 7. Running tests

The automated tests are in `tests/test_all.py`.

From the project root, run:

```bash
pytest -q
```

The test suite covers numerical root finding, derivatives, integration, unit conversion, critical compressibility-factor benchmarks, validation categories, route consistency, route infeasibility handling, knowledge-base dependencies, tolerance checks, uncertainty propagation, fluid-catalogue handling and structured request parsing.

The CoolProp comparison test uses `pytest.importorskip("CoolProp")`. If CoolProp is not installed, that particular test is skipped rather than treated as a successful independent reference comparison.

A passing internal route-consistency test checks agreement between implemented calculation paths within its stated tolerance. It is not an independent experimental validation of the EOS.

---

# 8. Exporting the CoolProp fluid database

The script `scripts/export_coolprop_database.py` exports the installed CoolProp catalogue and its available constants/metadata.

From the project root, run:

```bash
python scripts/export_coolprop_database.py
```

The script writes:

```text
data/coolprop_fluid_database.json
```

The JSON file records the CoolProp version, number of available fluids, successfully exported fluid records and any per-fluid failures. The output depends on the installed CoolProp build and is generated locally; it is not a fixed universal database.

---

# 9. Numerical method and route selection summary

| Task | Available implementation |
|---|---|
| Cubic polynomial roots | NumPy polynomial roots |
| Root refinement | Bisection, Brent, Newton–Raphson, Secant |
| Derivatives | Central, forward, five-point central finite differences |
| Integration | Adaptive quadrature, Simpson's rule |
| Symbolic EOS derivatives | SymPy |
| Unit conversion | Pint |
| Route feasibility | `routes.feasibility()` |
| Recommended route | `planner.analyse()` |
| Independent property reference | CoolProp, when installed |

The route planner does not use a fixed route regardless of state. It first evaluates route-specific feasibility. It then recommends the feasible route with the highest priority according to the configured cost ranking: Route 1, then Route 3, then Route 2. The user can override the recommendation, but an infeasible route produces an explicit error and is not silently replaced.

---

# 10. Important assumptions and limitations

1. The current state workflow is for single-phase pure-fluid calculations.
2. The five EOS models are the project's own cubic-EOS implementations. They are not CoolProp EOS replacements.
3. The bundled fluid records contain approximate/project-supplied constants and Cp polynomial coefficients. Their provenance and applicability should be considered when interpreting results.
4. The bundled Cp polynomial range is documented as approximately 250–1500 K; the engine warns when an energy-property calculation with bundled Cp data is outside that range.
5. The CoolProp catalogue depends on the installed CoolProp version.
6. CoolProp raw JSON is preserved as external metadata; it is not automatically transformed into the project's cubic EOS coefficients.
7. Non-bundled fluids from the installed CoolProp catalogue can be used with the project's cubic EOS parameter builder when CoolProp supplies the required constants. Their energy-property routes also require CoolProp ideal-gas Cp support.
8. Route 2 is restricted by its implemented vapour-side pressure path and is not used when the path would cross a subcritical saturation region.
9. Route 1's current closed-form residual expression is not marked feasible for the van der Waals parameter combination.
10. The current engine does not implement general inverse state solving from energy properties.
11. A successful route cross-check is an internal consistency check, not an independent validation.
12. CoolProp reference values are separate from the selected EOS result and do not overwrite it.
13. The uncertainty module provides a helper; it does not invent measurement uncertainty inputs.
14. Numerical convergence and physical-model accuracy are different matters. A converged numerical solution does not prove that the selected EOS accurately represents experimental behaviour.

---

# 11. File-by-file quick reference

| File | Main responsibility |
|---|---|
| `app.py` | Streamlit user interface and workflow |
| `thermo/engine.py` | Validation, state preparation, planning and calculation orchestration |
| `thermo/eos.py` | Cubic EOS equations, parameters, roots and phase/root selection |
| `thermo/fluids.py` | Bundled 22-fluid data, CoolProp-extended catalogue and fluid lookup |
| `thermo/coolprop_data.py` | CoolProp catalogue, properties, metadata and export |
| `thermo/routes.py` | Three enthalpy/entropy computational routes |
| `thermo/planner.py` | Route feasibility ranking and recommendation |
| `thermo/numerics.py` | Root solvers, finite differences and integration algorithms |
| `thermo/units.py` | Pint-based unit conversion |
| `thermo/knowledge_base.py` | Thermodynamic relations and property dependencies |
| `thermo/symbolic.py` | SymPy symbolic operations |
| `thermo/uncertainty.py` | Tolerance validation and uncertainty propagation helper |
| `thermo/request_parser.py` | Calculation request structure checks |
| `thermo/errors.py` | Structured error categories and exceptions |
| `scripts/export_coolprop_database.py` | CoolProp database export command |
| `tests/test_all.py` | Automated unit and integration-oriented tests |
| `docs/PROJECT_EXECUTION_ARCHITECTURE.md` | Additional architecture and execution notes |
| `requirements.txt` | Python dependency list |
| `.streamlit/config.toml` | Streamlit application configuration |

---

## Project summary

THERMO combines a Streamlit interface, a bundled pure-fluid dataset, five cubic equations of state, compressibility-root analysis, phase/root filtering, unit conversion, thermodynamic property relationships, three computational routes for enthalpy and entropy, numerical algorithms and optional CoolProp reference data.

Its main design principle is to keep the selected EOS calculation explicit, evaluate route feasibility before calculation, explain why a route is recommended, and report unsupported conditions rather than silently substituting another model or method.


# Free-text input (quantity parser)

`thermo/quantity_parser.py` turns a sentence such as `Methane at 25 C and 10 bar` into the `(property, value, unit)` tuples and fluid name that `engine.validate` expects. In the UI it is the *Fill from text* box under *Quick start*.

Rules it follows:

- A number becomes an input only if its **dimensions** identify one property (K -> T, Pa -> P, m^3/kg -> v, kg/m^3 -> rho, J/(kg K) -> s). Bare `C` / `F` are read as degC / degF (Pint would otherwise read coulomb / farad).
- Molar volume or molar density is converted to v / rho only when a fluid is identified (uses its molar mass).
- J/kg needs a label (`enthalpy 300 kJ/kg`); masses, total volumes and anything else are reported as *ignored*, never guessed.
- The same property stated with different values (e.g. two temperatures) is **not** filled; the user is told and enters it manually.
- Targets, EOS and phase are never inferred from text, and all parsed values still go through `engine.validate`.

Limits: rule-based, not NLP. It cannot tell *which* of two temperatures is the state point, ignores comma-grouped numbers (`1,000`), and uses a stopword list to stop words like "at"/"in" being read as units.
