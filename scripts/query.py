#!/usr/bin/env python3
"""
CLI query tool for the Periodically periodic table database.

Usage:
  python query.py element <name_or_symbol>
  python query.py search [--property-above VALUE] [--property-below VALUE] [--category CAT]
  python query.py reactions [--element SYMBOL] [--category CAT]
  python query.py compare <symbol1> <symbol2> [symbol3...]
  python query.py stats
"""

import argparse
import json
import sys
from pathlib import Path

try:
    from tabulate import tabulate
except ImportError:
    tabulate = None

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ELEMENTS_DIR = PROJECT_ROOT / "elements"
REACTIONS_DIR = PROJECT_ROOT / "reactions"
INDEXES_DIR = PROJECT_ROOT / "indexes"


def load_element(identifier):
    """Load element by name, symbol, or atomic number."""
    identifier = str(identifier).strip()

    # Try atomic number
    try:
        z = int(identifier)
        for path in ELEMENTS_DIR.glob(f"{z:03d}-*.json"):
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
    except ValueError:
        pass

    # Try symbol or name (case-insensitive)
    for path in sorted(ELEMENTS_DIR.glob("*.json")):
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if (data["symbol"].lower() == identifier.lower() or
                data["name"].lower() == identifier.lower()):
            return data

    return None


def load_all_elements():
    """Load all elements."""
    elements = []
    for path in sorted(ELEMENTS_DIR.glob("*.json")):
        with open(path, "r", encoding="utf-8") as f:
            elements.append(json.load(f))
    return elements


def load_reactions():
    """Load all reaction category files."""
    reactions = []
    for path in sorted(REACTIONS_DIR.glob("*.json")):
        if path.name == "index.json":
            continue
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, list):
            reactions.extend(data)
        elif isinstance(data, dict) and "reactions" in data:
            reactions.extend(data["reactions"])
    return reactions


def format_value(val, unit=""):
    """Format a value for display."""
    if val is None:
        return "—"
    if isinstance(val, float):
        if abs(val) >= 1000:
            return f"{val:,.1f}{unit}"
        elif abs(val) < 0.01 and val != 0:
            return f"{val:.2e}{unit}"
        else:
            return f"{val}{unit}"
    return str(val) + unit


def cmd_element(args):
    """Display full element data."""
    elem = load_element(args.identifier)
    if not elem:
        print(f"Element not found: {args.identifier}")
        return 1

    z = elem["atomic_number"]
    print(f"\n{'=' * 60}")
    print(f"  {elem['name']} ({elem['symbol']}) — Element {z}")
    print(f"{'=' * 60}")

    # Classification
    c = elem["classification"]
    print(f"\n  Category:   {c['category']}")
    print(f"  Group:      {c['group'] or 'N/A'}, Period: {c['period']}, Block: {c['block']}")
    print(f"  Occurrence: {c.get('natural_occurrence', 'N/A')}")

    # Atomic structure
    a = elem["atomic_structure"]
    print(f"\n  Electron config:  {a['electron_configuration']}")
    print(f"  Shells:           {a['electron_shells']}")
    print(f"  Valence e-:       {a.get('valence_electrons', 'N/A')}")
    print(f"  Oxidation states: {a.get('oxidation_states', [])}")
    print(f"  Electronegativity:{format_value(a.get('electronegativity_pauling'), ' (Pauling)')}")
    if a.get("ionization_energies_kj_mol"):
        print(f"  1st ionization:   {a['ionization_energies_kj_mol'][0]} kJ/mol")
    print(f"  Electron affinity:{format_value(a.get('electron_affinity_kj_mol'), ' kJ/mol')}")
    print(f"  Atomic radius:    {format_value(a.get('atomic_radius_pm'), ' pm')}")

    # Physical properties
    p = elem["physical_properties"]
    print(f"\n  Phase at STP:     {p['phase_at_stp']}")
    print(f"  Melting point:    {format_value(p.get('melting_point_k'), ' K')}")
    print(f"  Boiling point:    {format_value(p.get('boiling_point_k'), ' K')}")
    print(f"  Density:          {format_value(p.get('density_kg_m3'), ' kg/m3')}")
    print(f"  Molar heat cap:   {format_value(p.get('molar_heat_capacity_j_mol_k'), ' J/(mol·K)')}")
    print(f"  Crystal structure:{format_value(p.get('crystal_structure'))}")
    print(f"  Magnetic ordering:{format_value(p.get('magnetic_ordering'))}")
    print(f"  Thermal cond:     {format_value(p.get('thermal_conductivity_w_m_k'), ' W/(m·K)')}")

    # Nuclear
    n = elem["nuclear_properties"]
    print(f"\n  Radioactive:      {'Yes' if n['radioactive'] else 'No'}")
    if n["radioactive"]:
        print(f"  Half-life:        {n.get('half_life', 'N/A')}")
        print(f"  Decay mode:       {n.get('decay_mode', 'N/A')}")
    stable = n.get("stable_isotopes", [])
    print(f"  Stable isotopes:  {stable if stable else 'None'}")

    # Discovery
    d = elem["discovery"]
    print(f"\n  Discovered:       {d.get('year') or 'Ancient'}")
    print(f"  By:               {', '.join(d.get('discoverers', []))}")
    print(f"  Name origin:      {d.get('name_origin', 'N/A')}")

    # Applications
    apps = elem.get("applications", [])
    if apps:
        print(f"\n  Applications ({len(apps)}):")
        for app in apps[:8]:
            print(f"    - {app}")
        if len(apps) > 8:
            print(f"    ... and {len(apps) - 8} more")

    # Reactions
    rxns = elem.get("reactions", [])
    if rxns:
        print(f"\n  Reactions ({len(rxns)}):")
        for rxn in rxns[:5]:
            print(f"    [{rxn['category']}] {rxn['name']}: {rxn['equation']}")
        if len(rxns) > 5:
            print(f"    ... and {len(rxns) - 5} more")

    print()
    return 0


def cmd_search(args):
    """Search elements by property filters."""
    elements = load_all_elements()
    results = elements

    # Property filters
    prop_map = {
        "melting_point": ("physical_properties", "melting_point_k"),
        "boiling_point": ("physical_properties", "boiling_point_k"),
        "density": ("physical_properties", "density_kg_m3"),
        "electronegativity": ("atomic_structure", "electronegativity_pauling"),
        "atomic_mass": (None, "atomic_mass_u"),
        "atomic_radius": ("atomic_structure", "atomic_radius_pm"),
        "thermal_conductivity": ("physical_properties", "thermal_conductivity_w_m_k"),
    }

    if args.above:
        prop_name, threshold = args.above.split("=") if "=" in args.above else (args.above.rsplit("-", 1)[0].replace("-", "_"), args.above.rsplit("-", 1)[1] if "-" in args.above else None)
        # Parse --above prop=value format
        parts = args.above.split("=")
        if len(parts) == 2:
            prop_key = parts[0].replace("-", "_")
            threshold = float(parts[1])
        else:
            print("Usage: --above property=value (e.g., --above melting_point=3000)")
            return 1

        if prop_key in prop_map:
            section, field = prop_map[prop_key]
            filtered = []
            for e in results:
                if section:
                    val = e.get(section, {}).get(field)
                else:
                    val = e.get(field)
                if val is not None and val > threshold:
                    filtered.append(e)
            results = filtered
        else:
            print(f"Unknown property: {prop_key}")
            print(f"Available: {', '.join(prop_map.keys())}")
            return 1

    if args.below:
        parts = args.below.split("=")
        if len(parts) == 2:
            prop_key = parts[0].replace("-", "_")
            threshold = float(parts[1])
        else:
            print("Usage: --below property=value")
            return 1

        if prop_key in prop_map:
            section, field = prop_map[prop_key]
            filtered = []
            for e in results:
                if section:
                    val = e.get(section, {}).get(field)
                else:
                    val = e.get(field)
                if val is not None and val < threshold:
                    filtered.append(e)
            results = filtered

    if args.category:
        results = [e for e in results if e["classification"]["category"] == args.category]

    if args.phase:
        results = [e for e in results if e["physical_properties"]["phase_at_stp"] == args.phase]

    if args.block:
        results = [e for e in results if e["classification"]["block"] == args.block]

    if args.radioactive is not None:
        is_radio = args.radioactive.lower() in ("true", "yes", "1")
        results = [e for e in results if e["nuclear_properties"]["radioactive"] == is_radio]

    # Display results
    if not results:
        print("No elements match the criteria.")
        return 0

    headers = ["Z", "Symbol", "Name", "Category", "Mass (u)", "MP (K)", "BP (K)", "Density"]
    rows = []
    for e in results:
        rows.append([
            e["atomic_number"],
            e["symbol"],
            e["name"],
            e["classification"]["category"],
            format_value(e["atomic_mass_u"]),
            format_value(e["physical_properties"].get("melting_point_k")),
            format_value(e["physical_properties"].get("boiling_point_k")),
            format_value(e["physical_properties"].get("density_kg_m3")),
        ])

    if tabulate:
        print(tabulate(rows, headers=headers, tablefmt="simple"))
    else:
        # Fallback: simple formatted output
        print("  ".join(f"{h:>10}" for h in headers))
        print("-" * (12 * len(headers)))
        for row in rows:
            print("  ".join(f"{str(v):>10}" for v in row))

    print(f"\n{len(results)} element(s) found")
    return 0


def cmd_reactions(args):
    """Search and display reactions."""
    reactions = load_reactions()

    if args.element:
        reactions = [r for r in reactions
                     if args.element.upper() in [e.upper() for e in r.get("elements_involved", [])]]

    if args.category:
        reactions = [r for r in reactions if r.get("category") == args.category]

    if args.type:
        reactions = [r for r in reactions if r.get("type") == args.type]

    if not reactions:
        print("No reactions match the criteria.")
        return 0

    for rxn in reactions:
        thermo = rxn.get("thermodynamics", {})
        dh = thermo.get("delta_h_kj")
        exo = "exothermic" if thermo.get("exothermic") else "endothermic" if thermo.get("exothermic") is False else ""

        print(f"\n[{rxn['id']}] {rxn['name']}")
        print(f"  {rxn['equation']}")
        print(f"  Type: {rxn['type']} | Category: {rxn['category']}")
        print(f"  Elements: {', '.join(rxn.get('elements_involved', []))}")
        if dh is not None:
            print(f"  dH = {dh} kJ {exo}")
        cond = rxn.get("conditions", {})
        cond_parts = []
        if cond.get("temperature_k"):
            cond_parts.append(f"{cond['temperature_k']} K")
        if cond.get("pressure_atm") and cond["pressure_atm"] != 1.0:
            cond_parts.append(f"{cond['pressure_atm']} atm")
        if cond.get("catalyst"):
            cond_parts.append(f"catalyst: {cond['catalyst']}")
        if cond.get("other"):
            cond_parts.append(cond["other"])
        if cond_parts:
            print(f"  Conditions: {'; '.join(cond_parts)}")
        print(f"  {rxn.get('description', '')}")

    print(f"\n{len(reactions)} reaction(s) found")
    return 0


def cmd_compare(args):
    """Side-by-side comparison of elements."""
    elements = []
    for ident in args.elements:
        elem = load_element(ident)
        if not elem:
            print(f"Element not found: {ident}")
            return 1
        elements.append(elem)

    properties = [
        ("Atomic number", lambda e: e["atomic_number"]),
        ("Atomic mass (u)", lambda e: e["atomic_mass_u"]),
        ("Category", lambda e: e["classification"]["category"]),
        ("Group", lambda e: e["classification"]["group"]),
        ("Period", lambda e: e["classification"]["period"]),
        ("Block", lambda e: e["classification"]["block"]),
        ("Phase at STP", lambda e: e["physical_properties"]["phase_at_stp"]),
        ("Melting point (K)", lambda e: e["physical_properties"].get("melting_point_k")),
        ("Boiling point (K)", lambda e: e["physical_properties"].get("boiling_point_k")),
        ("Density (kg/m3)", lambda e: e["physical_properties"].get("density_kg_m3")),
        ("Electronegativity", lambda e: e["atomic_structure"].get("electronegativity_pauling")),
        ("1st ionization (kJ/mol)", lambda e: e["atomic_structure"]["ionization_energies_kj_mol"][0] if e["atomic_structure"]["ionization_energies_kj_mol"] else None),
        ("Electron affinity (kJ/mol)", lambda e: e["atomic_structure"].get("electron_affinity_kj_mol")),
        ("Atomic radius (pm)", lambda e: e["atomic_structure"].get("atomic_radius_pm")),
        ("Electron config", lambda e: e["atomic_structure"]["electron_configuration"]),
        ("Oxidation states", lambda e: str(e["atomic_structure"].get("common_oxidation_states", []))),
        ("Crystal structure", lambda e: e["physical_properties"].get("crystal_structure")),
        ("Magnetic ordering", lambda e: e["physical_properties"].get("magnetic_ordering")),
        ("Thermal conductivity", lambda e: e["physical_properties"].get("thermal_conductivity_w_m_k")),
        ("Radioactive", lambda e: "Yes" if e["nuclear_properties"]["radioactive"] else "No"),
    ]

    headers = ["Property"] + [f"{e['symbol']} ({e['name']})" for e in elements]
    rows = []
    for prop_name, getter in properties:
        row = [prop_name]
        for elem in elements:
            val = getter(elem)
            row.append(format_value(val))
        rows.append(row)

    if tabulate:
        print(tabulate(rows, headers=headers, tablefmt="grid"))
    else:
        col_widths = [max(len(str(row[i])) for row in [headers] + rows) for i in range(len(headers))]
        header_str = " | ".join(str(headers[i]).ljust(col_widths[i]) for i in range(len(headers)))
        print(header_str)
        print("-" * len(header_str))
        for row in rows:
            print(" | ".join(str(row[i]).ljust(col_widths[i]) for i in range(len(row))))

    return 0


def cmd_stats(args):
    """Show database statistics."""
    elements = load_all_elements()

    print(f"\n{'=' * 40}")
    print(f"  Periodically — Database Statistics")
    print(f"{'=' * 40}")
    print(f"\n  Elements: {len(elements)}")

    # Category breakdown
    cats = {}
    for e in elements:
        cat = e["classification"]["category"]
        cats[cat] = cats.get(cat, 0) + 1
    print(f"\n  Categories:")
    for cat, count in sorted(cats.items()):
        print(f"    {cat}: {count}")

    # Phase breakdown
    phases = {}
    for e in elements:
        phase = e["physical_properties"]["phase_at_stp"]
        phases[phase] = phases.get(phase, 0) + 1
    print(f"\n  Phases at STP:")
    for phase, count in sorted(phases.items()):
        print(f"    {phase}: {count}")

    # Radioactive count
    radioactive = sum(1 for e in elements if e["nuclear_properties"]["radioactive"])
    print(f"\n  Radioactive: {radioactive}")
    print(f"  Stable: {len(elements) - radioactive}")

    # Reaction count
    total_rxns = sum(len(e.get("reactions", [])) for e in elements)
    print(f"\n  Element-linked reactions: {total_rxns}")

    # Reaction files
    if REACTIONS_DIR.exists():
        rxn_files = list(REACTIONS_DIR.glob("*.json"))
        rxn_count = 0
        for path in rxn_files:
            if path.name == "index.json":
                continue
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            if isinstance(data, list):
                rxn_count += len(data)
            elif isinstance(data, dict) and "reactions" in data:
                rxn_count += len(data["reactions"])
        print(f"  Reaction files: {len(rxn_files)}")
        print(f"  Total reactions: {rxn_count}")

    # Coverage stats
    null_fields = ["electronegativity_pauling", "atomic_radius_pm", "melting_point_k",
                   "boiling_point_k", "density_kg_m3", "thermal_conductivity_w_m_k"]
    print(f"\n  Data coverage:")
    for field in null_fields:
        if field in ("electronegativity_pauling", "atomic_radius_pm"):
            has_val = sum(1 for e in elements if e["atomic_structure"].get(field) is not None)
        else:
            has_val = sum(1 for e in elements if e["physical_properties"].get(field) is not None)
        pct = has_val / len(elements) * 100
        print(f"    {field}: {has_val}/{len(elements)} ({pct:.0f}%)")

    print()
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Periodically — Periodic Table Query Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python query.py element hydrogen
  python query.py element Fe
  python query.py element 79
  python query.py search --above melting_point=3000
  python query.py search --category "noble gas"
  python query.py search --phase gas
  python query.py reactions --element Fe --category industrial
  python query.py compare H He Li
  python query.py stats""",
    )

    subparsers = parser.add_subparsers(dest="command", help="Command")

    # element
    p_elem = subparsers.add_parser("element", help="Display element details")
    p_elem.add_argument("identifier", help="Element name, symbol, or atomic number")

    # search
    p_search = subparsers.add_parser("search", help="Search elements by property")
    p_search.add_argument("--above", help="Filter property above value (property=value)")
    p_search.add_argument("--below", help="Filter property below value (property=value)")
    p_search.add_argument("--category", help="Filter by category")
    p_search.add_argument("--phase", help="Filter by phase at STP")
    p_search.add_argument("--block", help="Filter by block (s, p, d, f)")
    p_search.add_argument("--radioactive", help="Filter by radioactivity (true/false)")

    # reactions
    p_rxn = subparsers.add_parser("reactions", help="Search reactions")
    p_rxn.add_argument("--element", help="Filter by element symbol")
    p_rxn.add_argument("--category", help="Filter by category")
    p_rxn.add_argument("--type", help="Filter by reaction type")

    # compare
    p_cmp = subparsers.add_parser("compare", help="Compare elements side by side")
    p_cmp.add_argument("elements", nargs="+", help="Element names, symbols, or numbers")

    # stats
    subparsers.add_parser("stats", help="Show database statistics")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 0

    commands = {
        "element": cmd_element,
        "search": cmd_search,
        "reactions": cmd_reactions,
        "compare": cmd_compare,
        "stats": cmd_stats,
    }

    return commands[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
