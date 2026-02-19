#!/usr/bin/env python3
"""
Build index files from element and reaction data.

Generates:
  - indexes/periodic-table.json  (lightweight: number, symbol, name, mass, category)
  - indexes/by-category.json     (elements grouped by classification)
  - indexes/by-property.json     (min/max/avg for each numeric property)
"""

import json
import statistics
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ELEMENTS_DIR = PROJECT_ROOT / "elements"
INDEXES_DIR = PROJECT_ROOT / "indexes"


def load_all_elements():
    """Load all element JSON files, return list sorted by atomic number."""
    elements = []
    for path in sorted(ELEMENTS_DIR.glob("*.json")):
        with open(path, "r", encoding="utf-8") as f:
            elements.append(json.load(f))
    elements.sort(key=lambda e: e["atomic_number"])
    return elements


def build_periodic_table_index(elements):
    """Build lightweight periodic table index."""
    table = []
    for e in elements:
        table.append({
            "atomic_number": e["atomic_number"],
            "symbol": e["symbol"],
            "name": e["name"],
            "atomic_mass_u": e["atomic_mass_u"],
            "group": e["classification"]["group"],
            "period": e["classification"]["period"],
            "block": e["classification"]["block"],
            "category": e["classification"]["category"],
            "phase_at_stp": e["physical_properties"]["phase_at_stp"],
            "electronegativity_pauling": e["atomic_structure"].get("electronegativity_pauling"),
            "radioactive": e["nuclear_properties"]["radioactive"],
        })
    return table


def build_by_category_index(elements):
    """Group elements by their classification category."""
    categories = {}
    for e in elements:
        cat = e["classification"]["category"]
        if cat not in categories:
            categories[cat] = []
        categories[cat].append({
            "atomic_number": e["atomic_number"],
            "symbol": e["symbol"],
            "name": e["name"],
        })
    return categories


def build_by_property_index(elements):
    """Compute min/max/mean for each numeric property."""
    property_paths = {
        "atomic_mass_u": lambda e: e["atomic_mass_u"],
        "electronegativity_pauling": lambda e: e["atomic_structure"].get("electronegativity_pauling"),
        "first_ionization_energy_kj_mol": lambda e: e["atomic_structure"]["ionization_energies_kj_mol"][0] if e["atomic_structure"]["ionization_energies_kj_mol"] else None,
        "electron_affinity_kj_mol": lambda e: e["atomic_structure"].get("electron_affinity_kj_mol"),
        "atomic_radius_pm": lambda e: e["atomic_structure"].get("atomic_radius_pm"),
        "covalent_radius_pm": lambda e: e["atomic_structure"].get("covalent_radius_pm"),
        "melting_point_k": lambda e: e["physical_properties"].get("melting_point_k"),
        "boiling_point_k": lambda e: e["physical_properties"].get("boiling_point_k"),
        "density_kg_m3": lambda e: e["physical_properties"].get("density_kg_m3"),
        "thermal_conductivity_w_m_k": lambda e: e["physical_properties"].get("thermal_conductivity_w_m_k"),
        "heat_of_fusion_kj_mol": lambda e: e["physical_properties"].get("heat_of_fusion_kj_mol"),
        "heat_of_vaporization_kj_mol": lambda e: e["physical_properties"].get("heat_of_vaporization_kj_mol"),
        "molar_heat_capacity_j_mol_k": lambda e: e["physical_properties"].get("molar_heat_capacity_j_mol_k"),
    }

    props = {}
    for prop_name, getter in property_paths.items():
        values = []
        elements_with_val = []
        for e in elements:
            val = getter(e)
            if val is not None:
                values.append(val)
                elements_with_val.append({
                    "value": val,
                    "symbol": e["symbol"],
                    "name": e["name"],
                })

        if not values:
            continue

        elements_with_val.sort(key=lambda x: x["value"])

        props[prop_name] = {
            "min": {
                "value": elements_with_val[0]["value"],
                "element": elements_with_val[0]["symbol"],
            },
            "max": {
                "value": elements_with_val[-1]["value"],
                "element": elements_with_val[-1]["symbol"],
            },
            "mean": round(statistics.mean(values), 4),
            "median": round(statistics.median(values), 4),
            "count": len(values),
        }

    return props


def main():
    INDEXES_DIR.mkdir(parents=True, exist_ok=True)
    elements = load_all_elements()
    print(f"Loaded {len(elements)} elements")

    # Periodic table index
    pt = build_periodic_table_index(elements)
    pt_path = INDEXES_DIR / "periodic-table.json"
    with open(pt_path, "w", encoding="utf-8") as f:
        json.dump(pt, f, indent=2, ensure_ascii=False)
    print(f"  periodic-table.json ({len(pt)} entries)")

    # By category
    cats = build_by_category_index(elements)
    cat_path = INDEXES_DIR / "by-category.json"
    with open(cat_path, "w", encoding="utf-8") as f:
        json.dump(cats, f, indent=2, ensure_ascii=False)
    cat_summary = ", ".join(f"{k}: {len(v)}" for k, v in sorted(cats.items()))
    print(f"  by-category.json ({cat_summary})")

    # By property
    props = build_by_property_index(elements)
    prop_path = INDEXES_DIR / "by-property.json"
    with open(prop_path, "w", encoding="utf-8") as f:
        json.dump(props, f, indent=2, ensure_ascii=False)
    print(f"  by-property.json ({len(props)} properties)")

    print("\nDone!")


if __name__ == "__main__":
    main()
