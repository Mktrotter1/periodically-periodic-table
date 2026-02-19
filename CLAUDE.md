# Periodically — Machine-Readable Periodic Table Knowledge Base

## What This Is
Comprehensive JSON database of all 118 chemical elements with atomic structure, physical/chemical/thermodynamic/nuclear properties, isotopes, and categorized chemical reactions.

## Dual Purpose
1. **LLM context injection** — load one element file, get complete knowledge for chemistry questions
2. **Programmatic queries** — search elements by property, find reactions, compare materials via Python CLI

## Architecture
- `elements/NNN-name.json` — One file per element (118 total), single source of truth after generation
- `reactions/*.json` — Reactions grouped by category (industrial, laboratory, biological, environmental, notable)
- `indexes/*.json` — Lightweight lookup tables for fast filtering
- `schema/*.json` — JSON Schema for validation
- `scripts/` — Generation, validation, and query tools

## Data Pipeline
```
generate_elements.py  →  elements/*.json  (118 files)
generate_reactions.py →  reactions/*.json  (patches element files too)
build_indexes.py      →  indexes/*.json   (derived from element data)
validate.py           →  checks everything against schemas
query.py              →  CLI search/compare tool
```

## Key Commands
```bash
# Activate venv first
source .venv/bin/activate  # or .venv\Scripts\activate on Windows

# Regenerate everything
python scripts/generate_elements.py
python scripts/generate_reactions.py
python scripts/build_indexes.py

# Validate
python scripts/validate.py

# Query
python scripts/query.py element hydrogen
python scripts/query.py search --melting-point-above 3000
python scripts/query.py reactions --element Fe --category industrial
python scripts/query.py compare H He Li
```

## Unit Conventions
- Temperature: Kelvin (K)
- Mass: atomic mass units (u) for atoms, kg/m3 for density
- Energy: kJ/mol (thermodynamics)
- Distance: picometers (pm) for atomic radii
- Pressure: kPa or MPa (critical points), atm (reactions)
- Abundance: percent or ppm as standard for each context
- Null values for unknown/unmeasured properties

## Tech Stack
- Python 3.13
- jsonschema (validation)
- tabulate (CLI output formatting)

## File Naming
- Elements: `NNN-name.json` (e.g., `001-hydrogen.json`, `026-iron.json`)
- Zero-padded 3-digit atomic number + lowercase element name
