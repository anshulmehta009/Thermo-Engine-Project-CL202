"""Export every pure fluid and its available constants/EOS metadata from CoolProp.

Run from the project root:
    python scripts/export_coolprop_database.py
"""
from thermo.coolprop_data import export_database

if __name__ == "__main__":
    print(export_database("data/coolprop_fluid_database.json", include_raw_json=True))
