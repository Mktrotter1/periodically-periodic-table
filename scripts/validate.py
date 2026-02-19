#!/usr/bin/env python3
"""
Validate all element and reaction JSON files against their schemas.

Checks:
  - Schema compliance for all element files
  - Schema compliance for all reaction category files
  - Cross-references (reaction element refs match actual elements)
  - Missing required fields, null counts, coverage stats
"""

import json
import sys
from pathlib import Path

try:
    from jsonschema import validate, ValidationError, Draft202012Validator
except ImportError:
    print("ERROR: jsonschema not installed. Run: pip install jsonschema")
    sys.exit(1)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ELEMENTS_DIR = PROJECT_ROOT / "elements"
REACTIONS_DIR = PROJECT_ROOT / "reactions"
SCHEMA_DIR = PROJECT_ROOT / "schema"


def load_schema(name):
    path = SCHEMA_DIR / name
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def validate_elements(element_schema):
    """Validate all element files against schema."""
    errors = []
    stats = {"total": 0, "valid": 0, "null_counts": {}}
    elements = {}

    for path in sorted(ELEMENTS_DIR.glob("*.json")):
        stats["total"] += 1
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        elements[data["symbol"]] = data

        try:
            validate(instance=data, schema=element_schema)
            stats["valid"] += 1
        except ValidationError as e:
            errors.append(f"  {path.name}: {e.message} (path: {'.'.join(str(p) for p in e.absolute_path)})")

        # Count null fields for coverage reporting
        count_nulls(data, "", stats["null_counts"])

    return errors, stats, elements


def count_nulls(obj, prefix, counts):
    """Recursively count null values in nested dicts."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            full_key = f"{prefix}.{k}" if prefix else k
            if v is None:
                counts[full_key] = counts.get(full_key, 0) + 1
            elif isinstance(v, dict):
                count_nulls(v, full_key, counts)
    elif isinstance(obj, list):
        for item in obj:
            if isinstance(item, dict):
                count_nulls(item, prefix, counts)


def validate_reactions(reaction_schema):
    """Validate all reaction category files."""
    errors = []
    stats = {"total": 0, "valid": 0, "categories": {}}
    all_reaction_ids = set()
    all_element_refs = set()

    for path in sorted(REACTIONS_DIR.glob("*.json")):
        if path.name == "index.json":
            continue

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        category = path.stem
        reactions = data if isinstance(data, list) else data.get("reactions", [])
        stats["categories"][category] = len(reactions)

        for rxn in reactions:
            stats["total"] += 1
            try:
                validate(instance=rxn, schema=reaction_schema)
                stats["valid"] += 1
            except ValidationError as e:
                errors.append(f"  {path.name}/{rxn.get('id', '?')}: {e.message}")

            rxn_id = rxn.get("id")
            if rxn_id:
                if rxn_id in all_reaction_ids:
                    errors.append(f"  Duplicate reaction ID: {rxn_id}")
                all_reaction_ids.add(rxn_id)

            for elem in rxn.get("elements_involved", []):
                all_element_refs.add(elem)

    return errors, stats, all_element_refs


def check_cross_references(element_refs, known_symbols):
    """Check that all element references in reactions point to real elements."""
    errors = []
    for ref in sorted(element_refs):
        if ref not in known_symbols:
            errors.append(f"  Reaction references unknown element: {ref}")
    return errors


def main():
    print("=" * 60)
    print("Periodically — Data Validation Report")
    print("=" * 60)

    all_errors = []

    # Validate elements
    print("\n--- Element Validation ---")
    element_schema = load_schema("element.schema.json")
    elem_errors, elem_stats, elements = validate_elements(element_schema)
    all_errors.extend(elem_errors)

    print(f"Files: {elem_stats['total']}")
    print(f"Valid: {elem_stats['valid']}/{elem_stats['total']}")
    if elem_errors:
        print(f"Errors ({len(elem_errors)}):")
        for e in elem_errors:
            print(e)

    # Null coverage report
    null_counts = elem_stats["null_counts"]
    if null_counts:
        top_nulls = sorted(null_counts.items(), key=lambda x: -x[1])[:15]
        print(f"\nTop null fields (of {len(null_counts)} with nulls):")
        for field, count in top_nulls:
            pct = count / elem_stats["total"] * 100
            print(f"  {field}: {count}/{elem_stats['total']} ({pct:.0f}% missing)")

    # Validate reactions
    print("\n--- Reaction Validation ---")
    reaction_schema_path = SCHEMA_DIR / "reaction.schema.json"
    if reaction_schema_path.exists():
        reaction_schema = load_schema("reaction.schema.json")

        if REACTIONS_DIR.exists() and any(REACTIONS_DIR.glob("*.json")):
            rxn_errors, rxn_stats, element_refs = validate_reactions(reaction_schema)
            all_errors.extend(rxn_errors)

            print(f"Reactions: {rxn_stats['total']}")
            print(f"Valid: {rxn_stats['valid']}/{rxn_stats['total']}")
            for cat, count in sorted(rxn_stats["categories"].items()):
                print(f"  {cat}: {count} reactions")

            if rxn_errors:
                print(f"Errors ({len(rxn_errors)}):")
                for e in rxn_errors:
                    print(e)

            # Cross-reference check
            known_symbols = set(elements.keys())
            xref_errors = check_cross_references(element_refs, known_symbols)
            if xref_errors:
                all_errors.extend(xref_errors)
                print(f"\nCross-reference errors ({len(xref_errors)}):")
                for e in xref_errors:
                    print(e)
            else:
                print(f"\nCross-references: All {len(element_refs)} element refs valid")
        else:
            print("No reaction files found (run generate_reactions.py first)")
    else:
        print("No reaction schema found")

    # Summary
    print("\n" + "=" * 60)
    if all_errors:
        print(f"VALIDATION FAILED: {len(all_errors)} error(s)")
        return 1
    else:
        print("VALIDATION PASSED")
        return 0


if __name__ == "__main__":
    sys.exit(main())
