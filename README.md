# Periodically — Machine-Readable Periodic Table Knowledge Base

Comprehensive JSON database of all 118 chemical elements with atomic structure, physical/chemical/thermodynamic/nuclear properties, isotopes, and 154 categorized chemical reactions.

## Purpose

1. **LLM context injection** — load one element file, get complete knowledge for chemistry/materials questions
2. **Programmatic queries** — search elements by property, find reactions, compare materials via Python CLI

## Quick Start

```bash
python -m venv .venv
source .venv/bin/activate  # .venv\Scripts\activate on Windows
pip install -r requirements.txt

# Query an element
python scripts/query.py element iron
python scripts/query.py element 79

# Search by property
python scripts/query.py search --above melting_point=3000
python scripts/query.py search --category "noble gas"
python scripts/query.py search --phase gas --block p

# Find reactions
python scripts/query.py reactions --element Fe --category industrial
python scripts/query.py reactions --type combustion

# Compare elements
python scripts/query.py compare H He Li

# Database stats
python scripts/query.py stats
```

## Data Coverage

| Category | Count |
|----------|-------|
| Elements | 118 |
| Reactions | 154 (industrial: 52, laboratory: 36, notable: 30, biological: 18, environmental: 18) |
| Elements with reactions | 45 |
| Properties per element | ~50 fields |
| Isotope data | Major isotopes for all elements |

## Regenerating Data

The Python scripts in `scripts/` are the single source of truth:

```bash
python scripts/generate_elements.py   # Regenerate all 118 element files
python scripts/generate_reactions.py   # Regenerate reactions + inject into elements
python scripts/build_indexes.py        # Rebuild index files
python scripts/validate.py             # Validate everything against schemas
```

## Structure

```
elements/           118 JSON files (001-hydrogen.json ... 118-oganesson.json)
reactions/          Category files + master index
indexes/            Lightweight lookup tables
schema/             JSON Schema for validation
scripts/            Generation, validation, and query tools
```

## Unit Conventions

- Temperature: Kelvin (K)
- Mass: atomic mass units (u), density in kg/m3
- Energy: kJ/mol
- Distance: picometers (pm) for radii
- Pressure: kPa/MPa (critical points), atm (reactions)
- Null for unknown/unmeasured properties
