#!/usr/bin/env python3
"""
Generate comprehensive chemical reaction JSON files for all categories.

This is the single source of truth for reaction data. Edit reactions here, then
regenerate files with: python scripts/generate_reactions.py

Output:
  reactions/industrial.json     - Industrial/manufacturing reactions
  reactions/laboratory.json     - Laboratory/analytical reactions
  reactions/biological.json     - Biochemical reactions
  reactions/environmental.json  - Environmental/geochemical reactions
  reactions/notable.json        - Notable/historic/extreme reactions
  reactions/index.json          - Master index mapping IDs to categories/elements
  elements/NNN-name.json        - Each element file gets its reactions[] populated

Data sources: NIST Chemistry WebBook, CRC Handbook of Chemistry and Physics,
Atkins' Physical Chemistry, Housecroft & Sharpe Inorganic Chemistry,
Shriver & Atkins' Inorganic Chemistry, Voet & Voet Biochemistry.
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ELEMENTS_DIR = PROJECT_ROOT / "elements"
REACTIONS_DIR = PROJECT_ROOT / "reactions"

# =============================================================================
# ELEMENT LOOKUP TABLES
# =============================================================================

ELEMENTS = {
    1: ("Hydrogen", "H"), 2: ("Helium", "He"), 3: ("Lithium", "Li"),
    4: ("Beryllium", "Be"), 5: ("Boron", "B"), 6: ("Carbon", "C"),
    7: ("Nitrogen", "N"), 8: ("Oxygen", "O"), 9: ("Fluorine", "F"),
    10: ("Neon", "Ne"), 11: ("Sodium", "Na"), 12: ("Magnesium", "Mg"),
    13: ("Aluminium", "Al"), 14: ("Silicon", "Si"), 15: ("Phosphorus", "P"),
    16: ("Sulfur", "S"), 17: ("Chlorine", "Cl"), 18: ("Argon", "Ar"),
    19: ("Potassium", "K"), 20: ("Calcium", "Ca"), 21: ("Scandium", "Sc"),
    22: ("Titanium", "Ti"), 23: ("Vanadium", "V"), 24: ("Chromium", "Cr"),
    25: ("Manganese", "Mn"), 26: ("Iron", "Fe"), 27: ("Cobalt", "Co"),
    28: ("Nickel", "Ni"), 29: ("Copper", "Cu"), 30: ("Zinc", "Zn"),
    31: ("Gallium", "Ga"), 32: ("Germanium", "Ge"), 33: ("Arsenic", "As"),
    34: ("Selenium", "Se"), 35: ("Bromine", "Br"), 36: ("Krypton", "Kr"),
    37: ("Rubidium", "Rb"), 38: ("Strontium", "Sr"), 39: ("Yttrium", "Y"),
    40: ("Zirconium", "Zr"), 41: ("Niobium", "Nb"), 42: ("Molybdenum", "Mo"),
    43: ("Technetium", "Tc"), 44: ("Ruthenium", "Ru"), 45: ("Rhodium", "Rh"),
    46: ("Palladium", "Pd"), 47: ("Silver", "Ag"), 48: ("Cadmium", "Cd"),
    49: ("Indium", "In"), 50: ("Tin", "Sn"), 51: ("Antimony", "Sb"),
    52: ("Tellurium", "Te"), 53: ("Iodine", "I"), 54: ("Xenon", "Xe"),
    55: ("Caesium", "Cs"), 56: ("Barium", "Ba"), 57: ("Lanthanum", "La"),
    58: ("Cerium", "Ce"), 59: ("Praseodymium", "Pr"), 60: ("Neodymium", "Nd"),
    61: ("Promethium", "Pm"), 62: ("Samarium", "Sm"), 63: ("Europium", "Eu"),
    64: ("Gadolinium", "Gd"), 65: ("Terbium", "Tb"), 66: ("Dysprosium", "Dy"),
    67: ("Holmium", "Ho"), 68: ("Erbium", "Er"), 69: ("Thulium", "Tm"),
    70: ("Ytterbium", "Yb"), 71: ("Lutetium", "Lu"), 72: ("Hafnium", "Hf"),
    73: ("Tantalum", "Ta"), 74: ("Tungsten", "W"), 75: ("Rhenium", "Re"),
    76: ("Osmium", "Os"), 77: ("Iridium", "Ir"), 78: ("Platinum", "Pt"),
    79: ("Gold", "Au"), 80: ("Mercury", "Hg"), 81: ("Thallium", "Tl"),
    82: ("Lead", "Pb"), 83: ("Bismuth", "Bi"), 84: ("Polonium", "Po"),
    85: ("Astatine", "At"), 86: ("Radon", "Rn"), 87: ("Francium", "Fr"),
    88: ("Radium", "Ra"), 89: ("Actinium", "Ac"), 90: ("Thorium", "Th"),
    91: ("Protactinium", "Pa"), 92: ("Uranium", "U"), 93: ("Neptunium", "Np"),
    94: ("Plutonium", "Pu"), 95: ("Americium", "Am"), 96: ("Curium", "Cm"),
    97: ("Berkelium", "Bk"), 98: ("Californium", "Cf"), 99: ("Einsteinium", "Es"),
    100: ("Fermium", "Fm"), 101: ("Mendelevium", "Md"), 102: ("Nobelium", "No"),
    103: ("Lawrencium", "Lr"), 104: ("Rutherfordium", "Rf"),
    105: ("Dubnium", "Db"), 106: ("Seaborgium", "Sg"), 107: ("Bohrium", "Bh"),
    108: ("Hassium", "Hs"), 109: ("Meitnerium", "Mt"),
    110: ("Darmstadtium", "Ds"), 111: ("Roentgenium", "Rg"),
    112: ("Copernicium", "Cn"), 113: ("Nihonium", "Nh"),
    114: ("Flerovium", "Fl"), 115: ("Moscovium", "Mc"),
    116: ("Livermorium", "Lv"), 117: ("Tennessine", "Ts"),
    118: ("Oganesson", "Og"),
}

# Reverse lookup: symbol -> (atomic_number, name)
SYMBOL_TO_ELEMENT = {sym: (z, name) for z, (name, sym) in ELEMENTS.items()}


def _r(rxn_id, name, equation, equation_latex, rxn_type, category,
       elements_involved, reactants, products, thermodynamics,
       conditions, reversible, description):
    """Build a reaction dict with full schema."""
    return {
        "id": rxn_id,
        "name": name,
        "equation": equation,
        "equation_latex": equation_latex,
        "type": rxn_type,
        "category": category,
        "elements_involved": elements_involved,
        "reactants": reactants,
        "products": products,
        "thermodynamics": thermodynamics,
        "conditions": conditions,
        "reversible": reversible,
        "description": description,
    }


def _thermo(delta_h=None, delta_g=None, delta_s=None, exothermic=None):
    """Build thermodynamics dict."""
    return {
        "delta_h_kj": delta_h,
        "delta_g_kj": delta_g,
        "delta_s_j_k": delta_s,
        "exothermic": exothermic,
    }


def _cond(temp_k=None, pressure_atm=None, catalyst=None, other=None):
    """Build conditions dict."""
    return {
        "temperature_k": temp_k,
        "pressure_atm": pressure_atm,
        "catalyst": catalyst,
        "other": other,
    }


def _reactant(formula, moles, state):
    """Build a reactant/product dict."""
    return {"formula": formula, "moles": moles, "state": state}


# =============================================================================
# INDUSTRIAL REACTIONS
# =============================================================================

INDUSTRIAL_REACTIONS = [
    # --- Hydrogen & Ammonia ---
    _r("H-industrial-001", "Hydrogen combustion",
       "2 H2(g) + O2(g) -> 2 H2O(l)",
       "2\\,\\text{H}_2(g) + \\text{O}_2(g) \\rightarrow 2\\,\\text{H}_2\\text{O}(l)",
       "combustion", "industrial", ["H", "O"],
       [_reactant("H2", 2, "g"), _reactant("O2", 1, "g")],
       [_reactant("H2O", 2, "l")],
       _thermo(-572, -474.4, -326.7, True),
       _cond(other="ignition or spark required"),
       False,
       "Highly exothermic combustion producing water."),

    _r("N-industrial-001", "Haber process (ammonia synthesis)",
       "N2(g) + 3 H2(g) <=> 2 NH3(g)",
       "\\text{N}_2(g) + 3\\,\\text{H}_2(g) \\rightleftharpoons 2\\,\\text{NH}_3(g)",
       "synthesis", "industrial", ["N", "H"],
       [_reactant("N2", 1, "g"), _reactant("H2", 3, "g")],
       [_reactant("NH3", 2, "g")],
       _thermo(-92.2, -33.3, -198.7, True),
       _cond(temp_k=723, pressure_atm=200, catalyst="iron with K2O/Al2O3 promoters"),
       True,
       "Primary industrial route to ammonia. ~150 million tonnes/year globally."),

    _r("H-industrial-002", "Steam methane reforming",
       "CH4(g) + H2O(g) -> CO(g) + 3 H2(g)",
       "\\text{CH}_4(g) + \\text{H}_2\\text{O}(g) \\rightarrow \\text{CO}(g) + 3\\,\\text{H}_2(g)",
       "reforming", "industrial", ["C", "H", "O"],
       [_reactant("CH4", 1, "g"), _reactant("H2O", 1, "g")],
       [_reactant("CO", 1, "g"), _reactant("H2", 3, "g")],
       _thermo(206, 142, 215, False),
       _cond(temp_k=1100, pressure_atm=20, catalyst="nickel on alumina"),
       True,
       "Primary industrial hydrogen production method. Produces syngas."),

    _r("H-industrial-003", "Water-gas shift reaction",
       "CO(g) + H2O(g) <=> CO2(g) + H2(g)",
       "\\text{CO}(g) + \\text{H}_2\\text{O}(g) \\rightleftharpoons \\text{CO}_2(g) + \\text{H}_2(g)",
       "shift", "industrial", ["C", "H", "O"],
       [_reactant("CO", 1, "g"), _reactant("H2O", 1, "g")],
       [_reactant("CO2", 1, "g"), _reactant("H2", 1, "g")],
       _thermo(-41.2, -28.6, -42.1, True),
       _cond(temp_k=623, catalyst="iron oxide (HTS) or Cu-ZnO-Al2O3 (LTS)"),
       True,
       "Converts CO and steam to CO2 and H2. Used after steam reforming."),

    _r("H-industrial-004", "Water electrolysis",
       "2 H2O(l) -> 2 H2(g) + O2(g)",
       "2\\,\\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{H}_2(g) + \\text{O}_2(g)",
       "electrolysis", "industrial", ["H", "O"],
       [_reactant("H2O", 2, "l")],
       [_reactant("H2", 2, "g"), _reactant("O2", 1, "g")],
       _thermo(572, 474.4, 326.7, False),
       _cond(other="electrolysis; PEM or alkaline electrolyzer"),
       False,
       "Green hydrogen production via electrolysis. Growing rapidly for decarbonization."),

    # --- Sulfuric Acid (Contact Process) ---
    _r("S-industrial-001", "Contact process step 1 (sulfur combustion)",
       "S(s) + O2(g) -> SO2(g)",
       "\\text{S}(s) + \\text{O}_2(g) \\rightarrow \\text{SO}_2(g)",
       "combustion", "industrial", ["S", "O"],
       [_reactant("S", 1, "s"), _reactant("O2", 1, "g")],
       [_reactant("SO2", 1, "g")],
       _thermo(-296.8, -300.1, 11.1, True),
       _cond(),
       False,
       "First step of the Contact process for sulfuric acid production."),

    _r("S-industrial-002", "Contact process step 2 (SO2 oxidation)",
       "2 SO2(g) + O2(g) <=> 2 SO3(g)",
       "2\\,\\text{SO}_2(g) + \\text{O}_2(g) \\rightleftharpoons 2\\,\\text{SO}_3(g)",
       "oxidation", "industrial", ["S", "O"],
       [_reactant("SO2", 2, "g"), _reactant("O2", 1, "g")],
       [_reactant("SO3", 2, "g")],
       _thermo(-198, -140, -187, True),
       _cond(temp_k=723, catalyst="vanadium(V) oxide (V2O5)"),
       True,
       "Key catalytic step of the Contact process. ~99.5% conversion."),

    _r("S-industrial-003", "Contact process step 3 (oleum formation)",
       "SO3(g) + H2SO4(l) -> H2S2O7(l)",
       "\\text{SO}_3(g) + \\text{H}_2\\text{SO}_4(l) \\rightarrow \\text{H}_2\\text{S}_2\\text{O}_7(l)",
       "absorption", "industrial", ["S", "O", "H"],
       [_reactant("SO3", 1, "g"), _reactant("H2SO4", 1, "l")],
       [_reactant("H2S2O7", 1, "l")],
       _thermo(-113, None, None, True),
       _cond(),
       False,
       "SO3 absorbed into concentrated H2SO4 to form oleum (fuming sulfuric acid)."),

    _r("S-industrial-004", "Oleum dilution to sulfuric acid",
       "H2S2O7(l) + H2O(l) -> 2 H2SO4(l)",
       "\\text{H}_2\\text{S}_2\\text{O}_7(l) + \\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{H}_2\\text{SO}_4(l)",
       "hydrolysis", "industrial", ["S", "O", "H"],
       [_reactant("H2S2O7", 1, "l"), _reactant("H2O", 1, "l")],
       [_reactant("H2SO4", 2, "l")],
       _thermo(-80, None, None, True),
       _cond(),
       False,
       "Oleum diluted with water to produce concentrated sulfuric acid."),

    # --- Nitric Acid (Ostwald Process) ---
    _r("N-industrial-002", "Ostwald process step 1 (ammonia oxidation)",
       "4 NH3(g) + 5 O2(g) -> 4 NO(g) + 6 H2O(g)",
       "4\\,\\text{NH}_3(g) + 5\\,\\text{O}_2(g) \\rightarrow 4\\,\\text{NO}(g) + 6\\,\\text{H}_2\\text{O}(g)",
       "oxidation", "industrial", ["N", "H", "O"],
       [_reactant("NH3", 4, "g"), _reactant("O2", 5, "g")],
       [_reactant("NO", 4, "g"), _reactant("H2O", 6, "g")],
       _thermo(-905.2, None, None, True),
       _cond(temp_k=1100, catalyst="platinum-rhodium gauze (90% Pt / 10% Rh)"),
       False,
       "First step of Ostwald process for nitric acid. ~95% yield with Pt-Rh catalyst."),

    _r("N-industrial-003", "Ostwald process step 2 (NO oxidation)",
       "2 NO(g) + O2(g) -> 2 NO2(g)",
       "2\\,\\text{NO}(g) + \\text{O}_2(g) \\rightarrow 2\\,\\text{NO}_2(g)",
       "oxidation", "industrial", ["N", "O"],
       [_reactant("NO", 2, "g"), _reactant("O2", 1, "g")],
       [_reactant("NO2", 2, "g")],
       _thermo(-114, -70.4, -146.5, True),
       _cond(),
       False,
       "Spontaneous oxidation of nitric oxide to nitrogen dioxide."),

    _r("N-industrial-004", "Ostwald process step 3 (nitric acid formation)",
       "3 NO2(g) + H2O(l) -> 2 HNO3(aq) + NO(g)",
       "3\\,\\text{NO}_2(g) + \\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{HNO}_3(aq) + \\text{NO}(g)",
       "disproportionation", "industrial", ["N", "O", "H"],
       [_reactant("NO2", 3, "g"), _reactant("H2O", 1, "l")],
       [_reactant("HNO3", 2, "aq"), _reactant("NO", 1, "g")],
       _thermo(-138, None, None, True),
       _cond(),
       False,
       "Final step producing nitric acid. The NO byproduct is recycled to step 2."),

    # --- Solvay Process (sodium carbonate) ---
    _r("Na-industrial-001", "Solvay process (overall)",
       "2 NaCl(aq) + CaCO3(s) -> Na2CO3(s) + CaCl2(aq)",
       "2\\,\\text{NaCl}(aq) + \\text{CaCO}_3(s) \\rightarrow \\text{Na}_2\\text{CO}_3(s) + \\text{CaCl}_2(aq)",
       "double displacement", "industrial", ["Na", "Cl", "Ca", "C", "O"],
       [_reactant("NaCl", 2, "aq"), _reactant("CaCO3", 1, "s")],
       [_reactant("Na2CO3", 1, "s"), _reactant("CaCl2", 1, "aq")],
       _thermo(None, None, None, None),
       _cond(other="multi-step process using NH3 as intermediate"),
       False,
       "Overall Solvay process for soda ash (Na2CO3). Key step involves NaHCO3 precipitation."),

    _r("Na-industrial-002", "Solvay process (key precipitation step)",
       "NaCl(aq) + NH3(aq) + CO2(g) + H2O(l) -> NaHCO3(s) + NH4Cl(aq)",
       "\\text{NaCl}(aq) + \\text{NH}_3(aq) + \\text{CO}_2(g) + \\text{H}_2\\text{O}(l) \\rightarrow \\text{NaHCO}_3(s) + \\text{NH}_4\\text{Cl}(aq)",
       "precipitation", "industrial", ["Na", "Cl", "N", "H", "C", "O"],
       [_reactant("NaCl", 1, "aq"), _reactant("NH3", 1, "aq"),
        _reactant("CO2", 1, "g"), _reactant("H2O", 1, "l")],
       [_reactant("NaHCO3", 1, "s"), _reactant("NH4Cl", 1, "aq")],
       _thermo(None, None, None, True),
       _cond(temp_k=288, other="CO2 bubbled through ammoniacal brine"),
       False,
       "NaHCO3 precipitates due to low solubility. Filtered and calcined to Na2CO3."),

    # --- Aluminium Production (Hall-Heroult) ---
    _r("Al-industrial-001", "Hall-Heroult process",
       "2 Al2O3(l) + 3 C(s) -> 4 Al(l) + 3 CO2(g)",
       "2\\,\\text{Al}_2\\text{O}_3(l) + 3\\,\\text{C}(s) \\rightarrow 4\\,\\text{Al}(l) + 3\\,\\text{CO}_2(g)",
       "electrolysis", "industrial", ["Al", "O", "C"],
       [_reactant("Al2O3", 2, "l"), _reactant("C", 3, "s")],
       [_reactant("Al", 4, "l"), _reactant("CO2", 3, "g")],
       _thermo(2170, None, None, False),
       _cond(temp_k=1233, other="electrolysis in molten cryolite (Na3AlF6) bath; carbon anodes consumed"),
       False,
       "Primary aluminium smelting. Consumes ~15 kWh per kg Al. Dissolved in cryolite at ~960C."),

    _r("Al-industrial-002", "Bayer process (alumina extraction)",
       "Al2O3 * 3H2O(s) + 2 NaOH(aq) -> 2 NaAlO2(aq) + 4 H2O(l)",
       "\\text{Al}_2\\text{O}_3 \\cdot 3\\text{H}_2\\text{O}(s) + 2\\,\\text{NaOH}(aq) \\rightarrow 2\\,\\text{NaAlO}_2(aq) + 4\\,\\text{H}_2\\text{O}(l)",
       "dissolution", "industrial", ["Al", "O", "Na", "H"],
       [_reactant("Al2O3*3H2O", 1, "s"), _reactant("NaOH", 2, "aq")],
       [_reactant("NaAlO2", 2, "aq"), _reactant("H2O", 4, "l")],
       _thermo(None, None, None, None),
       _cond(temp_k=523, pressure_atm=30, other="digestion of bauxite in hot caustic soda"),
       False,
       "Bayer process extracts alumina from bauxite ore. Precursor to Hall-Heroult."),

    # --- Iron & Steel ---
    _r("Fe-industrial-001", "Blast furnace iron smelting (overall)",
       "Fe2O3(s) + 3 CO(g) -> 2 Fe(l) + 3 CO2(g)",
       "\\text{Fe}_2\\text{O}_3(s) + 3\\,\\text{CO}(g) \\rightarrow 2\\,\\text{Fe}(l) + 3\\,\\text{CO}_2(g)",
       "reduction", "industrial", ["Fe", "O", "C"],
       [_reactant("Fe2O3", 1, "s"), _reactant("CO", 3, "g")],
       [_reactant("Fe", 2, "l"), _reactant("CO2", 3, "g")],
       _thermo(-24.8, None, None, True),
       _cond(temp_k=1773, other="blast furnace with coke as reductant"),
       False,
       "Primary iron ore reduction. CO generated in situ from coke + hot air."),

    _r("Fe-industrial-002", "Coke combustion in blast furnace",
       "2 C(s) + O2(g) -> 2 CO(g)",
       "2\\,\\text{C}(s) + \\text{O}_2(g) \\rightarrow 2\\,\\text{CO}(g)",
       "combustion", "industrial", ["C", "O"],
       [_reactant("C", 2, "s"), _reactant("O2", 1, "g")],
       [_reactant("CO", 2, "g")],
       _thermo(-221, -274.4, 179.4, True),
       _cond(temp_k=1773, other="blast furnace tuyere zone"),
       False,
       "Limited combustion of coke produces CO reductant gas in blast furnace."),

    _r("Fe-industrial-003", "Basic oxygen steelmaking",
       "2 C(s) + O2(g) -> 2 CO(g)",
       "2\\,\\text{C}(s) + \\text{O}_2(g) \\rightarrow 2\\,\\text{CO}(g)",
       "oxidation", "industrial", ["C", "O", "Fe"],
       [_reactant("C", 2, "s"), _reactant("O2", 1, "g")],
       [_reactant("CO", 2, "g")],
       _thermo(-221, None, None, True),
       _cond(temp_k=1923, other="pure O2 blown into molten pig iron; removes excess C as CO"),
       False,
       "Basic oxygen furnace (BOF) converts pig iron to steel by blowing O2. ~70% of world steel."),

    _r("Fe-industrial-004", "Slag formation in steelmaking",
       "CaO(s) + SiO2(s) -> CaSiO3(l)",
       "\\text{CaO}(s) + \\text{SiO}_2(s) \\rightarrow \\text{CaSiO}_3(l)",
       "fusion", "industrial", ["Ca", "Si", "O"],
       [_reactant("CaO", 1, "s"), _reactant("SiO2", 1, "s")],
       [_reactant("CaSiO3", 1, "l")],
       _thermo(-89, None, None, True),
       _cond(temp_k=1873, other="limestone flux in blast furnace"),
       False,
       "CaO flux combines with silica impurities to form slag."),

    # --- Chloralkali ---
    _r("Cl-industrial-001", "Chloralkali electrolysis",
       "2 NaCl(aq) + 2 H2O(l) -> Cl2(g) + 2 NaOH(aq) + H2(g)",
       "2\\,\\text{NaCl}(aq) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow \\text{Cl}_2(g) + 2\\,\\text{NaOH}(aq) + \\text{H}_2(g)",
       "electrolysis", "industrial", ["Na", "Cl", "H", "O"],
       [_reactant("NaCl", 2, "aq"), _reactant("H2O", 2, "l")],
       [_reactant("Cl2", 1, "g"), _reactant("NaOH", 2, "aq"), _reactant("H2", 1, "g")],
       _thermo(422, None, None, False),
       _cond(other="membrane cell electrolysis at ~3.5 V"),
       False,
       "Produces chlorine, sodium hydroxide, and hydrogen simultaneously."),

    # --- Titanium (Kroll Process) ---
    _r("Ti-industrial-001", "Kroll process step 1 (chlorination)",
       "TiO2(s) + 2 Cl2(g) + 2 C(s) -> TiCl4(l) + 2 CO(g)",
       "\\text{TiO}_2(s) + 2\\,\\text{Cl}_2(g) + 2\\,\\text{C}(s) \\rightarrow \\text{TiCl}_4(l) + 2\\,\\text{CO}(g)",
       "chlorination", "industrial", ["Ti", "Cl", "O", "C"],
       [_reactant("TiO2", 1, "s"), _reactant("Cl2", 2, "g"), _reactant("C", 2, "s")],
       [_reactant("TiCl4", 1, "l"), _reactant("CO", 2, "g")],
       _thermo(-258, None, None, True),
       _cond(temp_k=1173),
       False,
       "Rutile ore converted to titanium tetrachloride. First step of Kroll process."),

    _r("Ti-industrial-002", "Kroll process step 2 (reduction)",
       "TiCl4(l) + 2 Mg(l) -> Ti(s) + 2 MgCl2(l)",
       "\\text{TiCl}_4(l) + 2\\,\\text{Mg}(l) \\rightarrow \\text{Ti}(s) + 2\\,\\text{MgCl}_2(l)",
       "reduction", "industrial", ["Ti", "Cl", "Mg"],
       [_reactant("TiCl4", 1, "l"), _reactant("Mg", 2, "l")],
       [_reactant("Ti", 1, "s"), _reactant("MgCl2", 2, "l")],
       _thermo(-512, None, None, True),
       _cond(temp_k=1073, other="argon atmosphere to prevent oxidation"),
       False,
       "Magnesium reduces TiCl4 to titanium sponge. Batch process, very energy-intensive."),

    # --- Copper ---
    _r("Cu-industrial-001", "Copper matte smelting",
       "2 Cu2S(s) + 3 O2(g) -> 2 Cu2O(s) + 2 SO2(g)",
       "2\\,\\text{Cu}_2\\text{S}(s) + 3\\,\\text{O}_2(g) \\rightarrow 2\\,\\text{Cu}_2\\text{O}(s) + 2\\,\\text{SO}_2(g)",
       "roasting", "industrial", ["Cu", "S", "O"],
       [_reactant("Cu2S", 2, "s"), _reactant("O2", 3, "g")],
       [_reactant("Cu2O", 2, "s"), _reactant("SO2", 2, "g")],
       _thermo(-773, None, None, True),
       _cond(temp_k=1473),
       False,
       "Partial roasting of copper matte in converting furnace."),

    _r("Cu-industrial-002", "Copper blister formation",
       "Cu2S(s) + 2 Cu2O(s) -> 6 Cu(l) + SO2(g)",
       "\\text{Cu}_2\\text{S}(s) + 2\\,\\text{Cu}_2\\text{O}(s) \\rightarrow 6\\,\\text{Cu}(l) + \\text{SO}_2(g)",
       "reduction", "industrial", ["Cu", "S", "O"],
       [_reactant("Cu2S", 1, "s"), _reactant("Cu2O", 2, "s")],
       [_reactant("Cu", 6, "l"), _reactant("SO2", 1, "g")],
       _thermo(None, None, None, True),
       _cond(temp_k=1473),
       False,
       "Self-reduction of Cu2S by Cu2O to produce blister copper (~99% Cu)."),

    _r("Cu-industrial-003", "Copper electrorefining",
       "Cu(s, impure) -> Cu(s, pure)",
       "\\text{Cu}(s,\\,\\text{impure}) \\rightarrow \\text{Cu}(s,\\,\\text{pure})",
       "electrolysis", "industrial", ["Cu"],
       [_reactant("Cu", 1, "s")],
       [_reactant("Cu", 1, "s")],
       _thermo(None, None, None, None),
       _cond(other="electrolysis in CuSO4/H2SO4; impure Cu anode, pure Cu cathode; ~0.3V"),
       False,
       "Electrorefining produces 99.99% pure copper. Precious metals recovered from anode slime."),

    # --- Zinc ---
    _r("Zn-industrial-001", "Zinc roasting",
       "2 ZnS(s) + 3 O2(g) -> 2 ZnO(s) + 2 SO2(g)",
       "2\\,\\text{ZnS}(s) + 3\\,\\text{O}_2(g) \\rightarrow 2\\,\\text{ZnO}(s) + 2\\,\\text{SO}_2(g)",
       "roasting", "industrial", ["Zn", "S", "O"],
       [_reactant("ZnS", 2, "s"), _reactant("O2", 3, "g")],
       [_reactant("ZnO", 2, "s"), _reactant("SO2", 2, "g")],
       _thermo(-879, None, None, True),
       _cond(temp_k=1173),
       False,
       "Roasting sphalerite (ZnS) to zinc oxide for subsequent reduction."),

    _r("Zn-industrial-002", "Zinc smelting (carbothermic reduction)",
       "ZnO(s) + C(s) -> Zn(g) + CO(g)",
       "\\text{ZnO}(s) + \\text{C}(s) \\rightarrow \\text{Zn}(g) + \\text{CO}(g)",
       "reduction", "industrial", ["Zn", "O", "C"],
       [_reactant("ZnO", 1, "s"), _reactant("C", 1, "s")],
       [_reactant("Zn", 1, "g"), _reactant("CO", 1, "g")],
       _thermo(237, None, None, False),
       _cond(temp_k=1373, other="Imperial Smelting Process or retort"),
       False,
       "Carbothermic reduction of ZnO. Zinc vapor condensed and collected."),

    # --- Sulfur Recovery ---
    _r("S-industrial-005", "Claus process (main reaction)",
       "2 H2S(g) + SO2(g) -> 3 S(s) + 2 H2O(g)",
       "2\\,\\text{H}_2\\text{S}(g) + \\text{SO}_2(g) \\rightarrow 3\\,\\text{S}(s) + 2\\,\\text{H}_2\\text{O}(g)",
       "redox", "industrial", ["S", "H", "O"],
       [_reactant("H2S", 2, "g"), _reactant("SO2", 1, "g")],
       [_reactant("S", 3, "s"), _reactant("H2O", 2, "g")],
       _thermo(-146, None, None, True),
       _cond(temp_k=573, catalyst="alumina or titanium dioxide"),
       False,
       "Recovers elemental sulfur from H2S in natural gas/refinery streams. >95% recovery."),

    _r("S-industrial-006", "Frasch process sulfur melting",
       "S(s) -> S(l)",
       "\\text{S}(s) \\rightarrow \\text{S}(l)",
       "phase change", "industrial", ["S"],
       [_reactant("S", 1, "s")],
       [_reactant("S", 1, "l")],
       _thermo(1.7, None, None, False),
       _cond(temp_k=388, other="superheated water (165C) injected underground to melt sulfur deposits"),
       False,
       "Frasch process mines sulfur by melting underground deposits with superheated water."),

    # --- Fischer-Tropsch & Methanol ---
    _r("C-industrial-001", "Fischer-Tropsch synthesis (general)",
       "n CO(g) + (2n+1) H2(g) -> CnH(2n+2)(l) + n H2O(g)",
       "n\\,\\text{CO}(g) + (2n+1)\\,\\text{H}_2(g) \\rightarrow \\text{C}_n\\text{H}_{2n+2}(l) + n\\,\\text{H}_2\\text{O}(g)",
       "synthesis", "industrial", ["C", "H", "O"],
       [_reactant("CO", 1, "g"), _reactant("H2", 2, "g")],
       [_reactant("CH2", 1, "l"), _reactant("H2O", 1, "g")],
       _thermo(-165, None, None, True),
       _cond(temp_k=523, pressure_atm=25, catalyst="iron or cobalt"),
       False,
       "Converts syngas to liquid hydrocarbons. Used for gas-to-liquids (GTL) and coal-to-liquids."),

    _r("C-industrial-002", "Methanol synthesis",
       "CO(g) + 2 H2(g) <=> CH3OH(g)",
       "\\text{CO}(g) + 2\\,\\text{H}_2(g) \\rightleftharpoons \\text{CH}_3\\text{OH}(g)",
       "synthesis", "industrial", ["C", "H", "O"],
       [_reactant("CO", 1, "g"), _reactant("H2", 2, "g")],
       [_reactant("CH3OH", 1, "g")],
       _thermo(-90.5, -25.3, -219, True),
       _cond(temp_k=523, pressure_atm=80, catalyst="Cu/ZnO/Al2O3"),
       True,
       "Industrial methanol synthesis from syngas. ~100 million tonnes/year globally."),

    # --- Ethylene & Polymers ---
    _r("C-industrial-003", "Ethylene production (steam cracking)",
       "C2H6(g) -> C2H4(g) + H2(g)",
       "\\text{C}_2\\text{H}_6(g) \\rightarrow \\text{C}_2\\text{H}_4(g) + \\text{H}_2(g)",
       "cracking", "industrial", ["C", "H"],
       [_reactant("C2H6", 1, "g")],
       [_reactant("C2H4", 1, "g"), _reactant("H2", 1, "g")],
       _thermo(136.4, 101, 119, False),
       _cond(temp_k=1123, other="steam cracking; residence time <1 second"),
       False,
       "Thermal dehydrogenation of ethane to ethylene. Largest-volume petrochemical."),

    _r("C-industrial-004", "Polyethylene polymerization",
       "n C2H4(g) -> (C2H4)n(s)",
       "n\\,\\text{C}_2\\text{H}_4(g) \\rightarrow (\\text{C}_2\\text{H}_4)_n(s)",
       "polymerization", "industrial", ["C", "H"],
       [_reactant("C2H4", 1, "g")],
       [_reactant("(C2H4)n", 1, "s")],
       _thermo(-93, None, None, True),
       _cond(catalyst="Ziegler-Natta (TiCl4/AlEt3) or metallocene", other="various pressures by process type"),
       False,
       "Polymerization of ethylene to polyethylene. Most-produced plastic globally."),

    # --- Cement & Glass ---
    _r("Ca-industrial-001", "Cement clinker formation (limestone calcination)",
       "CaCO3(s) -> CaO(s) + CO2(g)",
       "\\text{CaCO}_3(s) \\rightarrow \\text{CaO}(s) + \\text{CO}_2(g)",
       "decomposition", "industrial", ["Ca", "C", "O"],
       [_reactant("CaCO3", 1, "s")],
       [_reactant("CaO", 1, "s"), _reactant("CO2", 1, "g")],
       _thermo(178, 130.4, 160.6, False),
       _cond(temp_k=1173, other="rotary kiln"),
       False,
       "Calcination of limestone. Major source of industrial CO2 emissions."),

    _r("Ca-industrial-002", "Portland cement clinker (belite formation)",
       "2 CaO(s) + SiO2(s) -> Ca2SiO4(s)",
       "2\\,\\text{CaO}(s) + \\text{SiO}_2(s) \\rightarrow \\text{Ca}_2\\text{SiO}_4(s)",
       "synthesis", "industrial", ["Ca", "Si", "O"],
       [_reactant("CaO", 2, "s"), _reactant("SiO2", 1, "s")],
       [_reactant("Ca2SiO4", 1, "s")],
       _thermo(-127, None, None, True),
       _cond(temp_k=1723, other="rotary kiln; belite (C2S) phase"),
       False,
       "Calcium silicate (belite) formation in cement kiln. Key hydraulic phase."),

    _r("Si-industrial-001", "Soda-lime glass production",
       "Na2CO3(s) + SiO2(s) -> Na2SiO3(l) + CO2(g)",
       "\\text{Na}_2\\text{CO}_3(s) + \\text{SiO}_2(s) \\rightarrow \\text{Na}_2\\text{SiO}_3(l) + \\text{CO}_2(g)",
       "fusion", "industrial", ["Na", "Si", "C", "O"],
       [_reactant("Na2CO3", 1, "s"), _reactant("SiO2", 1, "s")],
       [_reactant("Na2SiO3", 1, "l"), _reactant("CO2", 1, "g")],
       _thermo(None, None, None, False),
       _cond(temp_k=1773, other="glass furnace"),
       False,
       "Sodium carbonate lowers melting point of silica in glass manufacture."),

    _r("Si-industrial-002", "Silicon production (carbothermic reduction)",
       "SiO2(s) + 2 C(s) -> Si(l) + 2 CO(g)",
       "\\text{SiO}_2(s) + 2\\,\\text{C}(s) \\rightarrow \\text{Si}(l) + 2\\,\\text{CO}(g)",
       "reduction", "industrial", ["Si", "O", "C"],
       [_reactant("SiO2", 1, "s"), _reactant("C", 2, "s")],
       [_reactant("Si", 1, "l"), _reactant("CO", 2, "g")],
       _thermo(689, None, None, False),
       _cond(temp_k=2173, other="electric arc furnace"),
       False,
       "Metallurgical-grade silicon production. Purified further for semiconductor use."),

    # --- Chlorine Chemistry ---
    _r("Cl-industrial-002", "PVC production (vinyl chloride)",
       "C2H4(g) + Cl2(g) -> C2H4Cl2(l)",
       "\\text{C}_2\\text{H}_4(g) + \\text{Cl}_2(g) \\rightarrow \\text{C}_2\\text{H}_4\\text{Cl}_2(l)",
       "addition", "industrial", ["C", "H", "Cl"],
       [_reactant("C2H4", 1, "g"), _reactant("Cl2", 1, "g")],
       [_reactant("C2H4Cl2", 1, "l")],
       _thermo(-218, None, None, True),
       _cond(catalyst="FeCl3"),
       False,
       "Direct chlorination of ethylene to 1,2-dichloroethane (EDC). Precursor to PVC."),

    # --- Phosphorus ---
    _r("P-industrial-001", "Phosphorus production (electric arc furnace)",
       "2 Ca3(PO4)2(s) + 6 SiO2(s) + 10 C(s) -> P4(g) + 6 CaSiO3(l) + 10 CO(g)",
       "2\\,\\text{Ca}_3(\\text{PO}_4)_2(s) + 6\\,\\text{SiO}_2(s) + 10\\,\\text{C}(s) \\rightarrow \\text{P}_4(g) + 6\\,\\text{CaSiO}_3(l) + 10\\,\\text{CO}(g)",
       "reduction", "industrial", ["P", "Ca", "Si", "O", "C"],
       [_reactant("Ca3(PO4)2", 2, "s"), _reactant("SiO2", 6, "s"), _reactant("C", 10, "s")],
       [_reactant("P4", 1, "g"), _reactant("CaSiO3", 6, "l"), _reactant("CO", 10, "g")],
       _thermo(None, None, None, False),
       _cond(temp_k=1773, other="electric arc furnace; P4 vapor condensed under water"),
       False,
       "White phosphorus production from phosphate rock via electric arc furnace."),

    _r("P-industrial-002", "Phosphoric acid production (wet process)",
       "Ca3(PO4)2(s) + 3 H2SO4(aq) -> 2 H3PO4(aq) + 3 CaSO4(s)",
       "\\text{Ca}_3(\\text{PO}_4)_2(s) + 3\\,\\text{H}_2\\text{SO}_4(aq) \\rightarrow 2\\,\\text{H}_3\\text{PO}_4(aq) + 3\\,\\text{CaSO}_4(s)",
       "acid digestion", "industrial", ["P", "Ca", "S", "O", "H"],
       [_reactant("Ca3(PO4)2", 1, "s"), _reactant("H2SO4", 3, "aq")],
       [_reactant("H3PO4", 2, "aq"), _reactant("CaSO4", 3, "s")],
       _thermo(None, None, None, True),
       _cond(temp_k=353),
       False,
       "Wet-process phosphoric acid for fertilizer production. Gypsum (CaSO4) byproduct."),

    # --- Sodium Hydroxide & Misc ---
    _r("Na-industrial-003", "Sodium hydroxide + CO2 (caustic scrubbing)",
       "2 NaOH(aq) + CO2(g) -> Na2CO3(aq) + H2O(l)",
       "2\\,\\text{NaOH}(aq) + \\text{CO}_2(g) \\rightarrow \\text{Na}_2\\text{CO}_3(aq) + \\text{H}_2\\text{O}(l)",
       "absorption", "industrial", ["Na", "O", "H", "C"],
       [_reactant("NaOH", 2, "aq"), _reactant("CO2", 1, "g")],
       [_reactant("Na2CO3", 1, "aq"), _reactant("H2O", 1, "l")],
       _thermo(-109.4, None, None, True),
       _cond(),
       False,
       "Caustic scrubbing of CO2 from flue gas or process streams."),

    # --- Petroleum Refining ---
    _r("C-industrial-005", "Catalytic cracking (general)",
       "C16H34(l) -> C8H18(l) + C8H16(g)",
       "\\text{C}_{16}\\text{H}_{34}(l) \\rightarrow \\text{C}_8\\text{H}_{18}(l) + \\text{C}_8\\text{H}_{16}(g)",
       "cracking", "industrial", ["C", "H"],
       [_reactant("C16H34", 1, "l")],
       [_reactant("C8H18", 1, "l"), _reactant("C8H16", 1, "g")],
       _thermo(None, None, None, False),
       _cond(temp_k=773, catalyst="zeolite (USY or ZSM-5)", other="FCC unit; 2-5 second contact time"),
       False,
       "Fluid catalytic cracking (FCC) of heavy petroleum fractions into gasoline."),

    _r("C-industrial-006", "Methane combustion",
       "CH4(g) + 2 O2(g) -> CO2(g) + 2 H2O(g)",
       "\\text{CH}_4(g) + 2\\,\\text{O}_2(g) \\rightarrow \\text{CO}_2(g) + 2\\,\\text{H}_2\\text{O}(g)",
       "combustion", "industrial", ["C", "H", "O"],
       [_reactant("CH4", 1, "g"), _reactant("O2", 2, "g")],
       [_reactant("CO2", 1, "g"), _reactant("H2O", 2, "g")],
       _thermo(-802.3, -800.8, -5.1, True),
       _cond(other="ignition source required"),
       False,
       "Natural gas combustion. Cleanest-burning fossil fuel."),

    # --- Magnesium ---
    _r("Mg-industrial-001", "Pidgeon process (magnesium production)",
       "2 MgO(s) + Si(s) -> 2 Mg(g) + SiO2(s)",
       "2\\,\\text{MgO}(s) + \\text{Si}(s) \\rightarrow 2\\,\\text{Mg}(g) + \\text{SiO}_2(s)",
       "reduction", "industrial", ["Mg", "O", "Si"],
       [_reactant("MgO", 2, "s"), _reactant("Si", 1, "s")],
       [_reactant("Mg", 2, "g"), _reactant("SiO2", 1, "s")],
       _thermo(None, None, None, False),
       _cond(temp_k=1473, pressure_atm=0.01, other="vacuum retort; Mg vapor condensed"),
       False,
       "Pidgeon process: silicothermic reduction of MgO. Dominant in China."),

    # --- Chromium ---
    _r("Cr-industrial-001", "Chromium production (aluminothermic)",
       "Cr2O3(s) + 2 Al(s) -> 2 Cr(s) + Al2O3(s)",
       "\\text{Cr}_2\\text{O}_3(s) + 2\\,\\text{Al}(s) \\rightarrow 2\\,\\text{Cr}(s) + \\text{Al}_2\\text{O}_3(s)",
       "aluminothermic", "industrial", ["Cr", "O", "Al"],
       [_reactant("Cr2O3", 1, "s"), _reactant("Al", 2, "s")],
       [_reactant("Cr", 2, "s"), _reactant("Al2O3", 1, "s")],
       _thermo(-536, None, None, True),
       _cond(other="thermite-type reaction; used for ferrochrome alloy"),
       False,
       "Aluminothermic reduction of chromium oxide for chromium metal production."),

    # --- Ammonia oxidation for urea ---
    _r("N-industrial-005", "Urea synthesis",
       "2 NH3(g) + CO2(g) <=> NH2CONH2(s) + H2O(l)",
       "2\\,\\text{NH}_3(g) + \\text{CO}_2(g) \\rightleftharpoons \\text{NH}_2\\text{CONH}_2(s) + \\text{H}_2\\text{O}(l)",
       "synthesis", "industrial", ["N", "H", "C", "O"],
       [_reactant("NH3", 2, "g"), _reactant("CO2", 1, "g")],
       [_reactant("NH2CONH2", 1, "s"), _reactant("H2O", 1, "l")],
       _thermo(-103.5, None, None, True),
       _cond(temp_k=458, pressure_atm=175, other="two-step via ammonium carbamate"),
       True,
       "Bosch-Meiser process. World's most-produced chemical by mass (~180 Mt/y). Fertilizer."),

    # --- Lithium ---
    _r("Li-industrial-001", "Lithium carbonate to lithium metal",
       "Li2CO3(s) -> Li2O(s) + CO2(g)",
       "\\text{Li}_2\\text{CO}_3(s) \\rightarrow \\text{Li}_2\\text{O}(s) + \\text{CO}_2(g)",
       "decomposition", "industrial", ["Li", "C", "O"],
       [_reactant("Li2CO3", 1, "s")],
       [_reactant("Li2O", 1, "s"), _reactant("CO2", 1, "g")],
       _thermo(234, None, None, False),
       _cond(temp_k=1583),
       False,
       "Thermal decomposition of lithium carbonate. Precursor to lithium metal electrolysis."),

    # --- Manganese ---
    _r("Mn-industrial-001", "Manganese dioxide reduction (ferromanganese)",
       "MnO2(s) + 2 C(s) -> Mn(s) + 2 CO(g)",
       "\\text{MnO}_2(s) + 2\\,\\text{C}(s) \\rightarrow \\text{Mn}(s) + 2\\,\\text{CO}(g)",
       "reduction", "industrial", ["Mn", "O", "C"],
       [_reactant("MnO2", 1, "s"), _reactant("C", 2, "s")],
       [_reactant("Mn", 1, "s"), _reactant("CO", 2, "g")],
       _thermo(None, None, None, False),
       _cond(temp_k=1573, other="electric arc furnace for ferromanganese"),
       False,
       "Carbothermic reduction of pyrolusite. 90% of Mn used in steel production."),

    # --- Nickel ---
    _r("Ni-industrial-001", "Mond process (nickel purification)",
       "Ni(s) + 4 CO(g) <=> Ni(CO)4(g)",
       "\\text{Ni}(s) + 4\\,\\text{CO}(g) \\rightleftharpoons \\text{Ni}(\\text{CO})_4(g)",
       "carbonyl formation", "industrial", ["Ni", "C", "O"],
       [_reactant("Ni", 1, "s"), _reactant("CO", 4, "g")],
       [_reactant("Ni(CO)4", 1, "g")],
       _thermo(-160, None, None, True),
       _cond(temp_k=323, other="Ni(CO)4 decomposes at 230C to deposit pure Ni; extremely toxic gas"),
       True,
       "Mond process: volatile Ni(CO)4 formed at 50C, decomposed at 230C for 99.99% pure Ni."),

    # --- Lead ---
    _r("Pb-industrial-001", "Lead smelting",
       "2 PbS(s) + 3 O2(g) -> 2 PbO(s) + 2 SO2(g)",
       "2\\,\\text{PbS}(s) + 3\\,\\text{O}_2(g) \\rightarrow 2\\,\\text{PbO}(s) + 2\\,\\text{SO}_2(g)",
       "roasting", "industrial", ["Pb", "S", "O"],
       [_reactant("PbS", 2, "s"), _reactant("O2", 3, "g")],
       [_reactant("PbO", 2, "s"), _reactant("SO2", 2, "g")],
       _thermo(-844, None, None, True),
       _cond(temp_k=1173),
       False,
       "Roasting of galena (PbS) to lead oxide, then reduced with carbon to lead metal."),

    # --- Tin ---
    _r("Sn-industrial-001", "Tin smelting",
       "SnO2(s) + 2 C(s) -> Sn(l) + 2 CO(g)",
       "\\text{SnO}_2(s) + 2\\,\\text{C}(s) \\rightarrow \\text{Sn}(l) + 2\\,\\text{CO}(g)",
       "reduction", "industrial", ["Sn", "O", "C"],
       [_reactant("SnO2", 1, "s"), _reactant("C", 2, "s")],
       [_reactant("Sn", 1, "l"), _reactant("CO", 2, "g")],
       _thermo(None, None, None, False),
       _cond(temp_k=1473, other="reverberatory furnace; one of the oldest smelting processes"),
       False,
       "Carbothermic reduction of cassiterite. Bronze Age technology. Tin used in solder and tinplate."),
]


# =============================================================================
# LABORATORY REACTIONS
# =============================================================================

LABORATORY_REACTIONS = [
    # --- Acid-Base ---
    _r("H-laboratory-001", "Strong acid-base neutralization (HCl + NaOH)",
       "HCl(aq) + NaOH(aq) -> NaCl(aq) + H2O(l)",
       "\\text{HCl}(aq) + \\text{NaOH}(aq) \\rightarrow \\text{NaCl}(aq) + \\text{H}_2\\text{O}(l)",
       "neutralization", "laboratory", ["H", "Cl", "Na", "O"],
       [_reactant("HCl", 1, "aq"), _reactant("NaOH", 1, "aq")],
       [_reactant("NaCl", 1, "aq"), _reactant("H2O", 1, "l")],
       _thermo(-57.1, -80.7, 79.1, True),
       _cond(),
       False,
       "Fundamental acid-base neutralization. Enthalpy is -57.1 kJ/mol for all strong acid-strong base pairs."),

    _r("H-laboratory-002", "Sulfuric acid + sodium hydroxide",
       "H2SO4(aq) + 2 NaOH(aq) -> Na2SO4(aq) + 2 H2O(l)",
       "\\text{H}_2\\text{SO}_4(aq) + 2\\,\\text{NaOH}(aq) \\rightarrow \\text{Na}_2\\text{SO}_4(aq) + 2\\,\\text{H}_2\\text{O}(l)",
       "neutralization", "laboratory", ["H", "S", "O", "Na"],
       [_reactant("H2SO4", 1, "aq"), _reactant("NaOH", 2, "aq")],
       [_reactant("Na2SO4", 1, "aq"), _reactant("H2O", 2, "l")],
       _thermo(-114.2, None, None, True),
       _cond(),
       False,
       "Diprotic acid neutralization. Two moles of NaOH required per mole of H2SO4."),

    _r("H-laboratory-003", "Acetic acid + NaOH (weak acid titration)",
       "CH3COOH(aq) + NaOH(aq) -> CH3COONa(aq) + H2O(l)",
       "\\text{CH}_3\\text{COOH}(aq) + \\text{NaOH}(aq) \\rightarrow \\text{CH}_3\\text{COONa}(aq) + \\text{H}_2\\text{O}(l)",
       "neutralization", "laboratory", ["C", "H", "O", "Na"],
       [_reactant("CH3COOH", 1, "aq"), _reactant("NaOH", 1, "aq")],
       [_reactant("CH3COONa", 1, "aq"), _reactant("H2O", 1, "l")],
       _thermo(-55.8, None, None, True),
       _cond(other="phenolphthalein indicator; endpoint pH ~8.7"),
       False,
       "Weak acid-strong base titration. Equivalence point is basic due to conjugate base hydrolysis."),

    # --- Redox Titrations ---
    _r("Mn-laboratory-001", "Permanganate titration of iron(II)",
       "MnO4-(aq) + 5 Fe2+(aq) + 8 H+(aq) -> Mn2+(aq) + 5 Fe3+(aq) + 4 H2O(l)",
       "\\text{MnO}_4^-(aq) + 5\\,\\text{Fe}^{2+}(aq) + 8\\,\\text{H}^+(aq) \\rightarrow \\text{Mn}^{2+}(aq) + 5\\,\\text{Fe}^{3+}(aq) + 4\\,\\text{H}_2\\text{O}(l)",
       "redox titration", "laboratory", ["Mn", "Fe", "O", "H"],
       [_reactant("KMnO4", 1, "aq"), _reactant("FeSO4", 5, "aq"), _reactant("H2SO4", 4, "aq")],
       [_reactant("MnSO4", 1, "aq"), _reactant("Fe2(SO4)3", 2.5, "aq"), _reactant("H2O", 4, "l")],
       _thermo(None, None, None, None),
       _cond(other="self-indicating; purple MnO4- decolorizes until endpoint"),
       False,
       "Classic redox titration. MnO4- is its own indicator (purple to colorless)."),

    _r("Cr-laboratory-001", "Dichromate titration of iron(II)",
       "Cr2O7 2-(aq) + 6 Fe2+(aq) + 14 H+(aq) -> 2 Cr3+(aq) + 6 Fe3+(aq) + 7 H2O(l)",
       "\\text{Cr}_2\\text{O}_7^{2-}(aq) + 6\\,\\text{Fe}^{2+}(aq) + 14\\,\\text{H}^+(aq) \\rightarrow 2\\,\\text{Cr}^{3+}(aq) + 6\\,\\text{Fe}^{3+}(aq) + 7\\,\\text{H}_2\\text{O}(l)",
       "redox titration", "laboratory", ["Cr", "Fe", "O", "H"],
       [_reactant("K2Cr2O7", 1, "aq"), _reactant("FeSO4", 6, "aq"), _reactant("H2SO4", 7, "aq")],
       [_reactant("Cr2(SO4)3", 1, "aq"), _reactant("Fe2(SO4)3", 3, "aq"), _reactant("H2O", 7, "l")],
       _thermo(None, None, None, None),
       _cond(other="diphenylamine or ferroin indicator"),
       False,
       "Dichromate is a primary standard (unlike permanganate). Orange to green color change."),

    _r("I-laboratory-001", "Iodometric titration (thiosulfate)",
       "I2(aq) + 2 Na2S2O3(aq) -> 2 NaI(aq) + Na2S4O6(aq)",
       "\\text{I}_2(aq) + 2\\,\\text{Na}_2\\text{S}_2\\text{O}_3(aq) \\rightarrow 2\\,\\text{NaI}(aq) + \\text{Na}_2\\text{S}_4\\text{O}_6(aq)",
       "redox titration", "laboratory", ["I", "Na", "S", "O"],
       [_reactant("I2", 1, "aq"), _reactant("Na2S2O3", 2, "aq")],
       [_reactant("NaI", 2, "aq"), _reactant("Na2S4O6", 1, "aq")],
       _thermo(None, None, None, None),
       _cond(other="starch indicator turns blue-black with I2"),
       False,
       "Iodometric back-titration. Starch indicator; blue to colorless at endpoint."),

    # --- Precipitation ---
    _r("Ag-laboratory-001", "Silver chloride precipitation",
       "AgNO3(aq) + NaCl(aq) -> AgCl(s) + NaNO3(aq)",
       "\\text{AgNO}_3(aq) + \\text{NaCl}(aq) \\rightarrow \\text{AgCl}(s) + \\text{NaNO}_3(aq)",
       "precipitation", "laboratory", ["Ag", "N", "O", "Na", "Cl"],
       [_reactant("AgNO3", 1, "aq"), _reactant("NaCl", 1, "aq")],
       [_reactant("AgCl", 1, "s"), _reactant("NaNO3", 1, "aq")],
       _thermo(-65.5, -55.7, -33.0, True),
       _cond(),
       False,
       "Classic precipitation. White curdy precipitate. Ksp = 1.77e-10. Used in chloride analysis."),

    _r("Ba-laboratory-001", "Barium sulfate precipitation (sulfate test)",
       "BaCl2(aq) + Na2SO4(aq) -> BaSO4(s) + 2 NaCl(aq)",
       "\\text{BaCl}_2(aq) + \\text{Na}_2\\text{SO}_4(aq) \\rightarrow \\text{BaSO}_4(s) + 2\\,\\text{NaCl}(aq)",
       "precipitation", "laboratory", ["Ba", "Cl", "Na", "S", "O"],
       [_reactant("BaCl2", 1, "aq"), _reactant("Na2SO4", 1, "aq")],
       [_reactant("BaSO4", 1, "s"), _reactant("NaCl", 2, "aq")],
       _thermo(None, None, None, True),
       _cond(other="add BaCl2 in dilute HCl to acidified sample"),
       False,
       "Confirmatory test for sulfate ions. White precipitate insoluble in dilute HCl. Ksp = 1.1e-10."),

    _r("Pb-laboratory-001", "Lead(II) iodide precipitation",
       "Pb(NO3)2(aq) + 2 KI(aq) -> PbI2(s) + 2 KNO3(aq)",
       "\\text{Pb}(\\text{NO}_3)_2(aq) + 2\\,\\text{KI}(aq) \\rightarrow \\text{PbI}_2(s) + 2\\,\\text{KNO}_3(aq)",
       "precipitation", "laboratory", ["Pb", "N", "O", "K", "I"],
       [_reactant("Pb(NO3)2", 1, "aq"), _reactant("KI", 2, "aq")],
       [_reactant("PbI2", 1, "s"), _reactant("KNO3", 2, "aq")],
       _thermo(None, None, None, True),
       _cond(other="bright yellow precipitate; 'golden rain' demo when hot solution cooled"),
       False,
       "Golden rain experiment. Bright yellow PbI2 crystallizes from hot solution on cooling."),

    _r("Ca-laboratory-001", "Calcium carbonate precipitation",
       "CaCl2(aq) + Na2CO3(aq) -> CaCO3(s) + 2 NaCl(aq)",
       "\\text{CaCl}_2(aq) + \\text{Na}_2\\text{CO}_3(aq) \\rightarrow \\text{CaCO}_3(s) + 2\\,\\text{NaCl}(aq)",
       "precipitation", "laboratory", ["Ca", "Cl", "Na", "C", "O"],
       [_reactant("CaCl2", 1, "aq"), _reactant("Na2CO3", 1, "aq")],
       [_reactant("CaCO3", 1, "s"), _reactant("NaCl", 2, "aq")],
       _thermo(None, None, None, True),
       _cond(),
       False,
       "White precipitate of calcium carbonate. Dissolves in dilute acid with effervescence."),

    # --- Gas Generation ---
    _r("H-laboratory-004", "Hydrogen generation (Zn + HCl)",
       "Zn(s) + 2 HCl(aq) -> ZnCl2(aq) + H2(g)",
       "\\text{Zn}(s) + 2\\,\\text{HCl}(aq) \\rightarrow \\text{ZnCl}_2(aq) + \\text{H}_2(g)",
       "single displacement", "laboratory", ["Zn", "H", "Cl"],
       [_reactant("Zn", 1, "s"), _reactant("HCl", 2, "aq")],
       [_reactant("ZnCl2", 1, "aq"), _reactant("H2", 1, "g")],
       _thermo(-153, None, None, True),
       _cond(),
       False,
       "Standard lab preparation of hydrogen gas. Collect by downward displacement of water."),

    _r("C-laboratory-001", "CO2 generation (marble + HCl)",
       "CaCO3(s) + 2 HCl(aq) -> CaCl2(aq) + H2O(l) + CO2(g)",
       "\\text{CaCO}_3(s) + 2\\,\\text{HCl}(aq) \\rightarrow \\text{CaCl}_2(aq) + \\text{H}_2\\text{O}(l) + \\text{CO}_2(g)",
       "acid-carbonate", "laboratory", ["Ca", "C", "O", "H", "Cl"],
       [_reactant("CaCO3", 1, "s"), _reactant("HCl", 2, "aq")],
       [_reactant("CaCl2", 1, "aq"), _reactant("H2O", 1, "l"), _reactant("CO2", 1, "g")],
       _thermo(-16, None, None, True),
       _cond(),
       False,
       "Reaction of marble chips (CaCO3) with hydrochloric acid. Common CO2 source in lab."),

    _r("O-laboratory-001", "Oxygen from hydrogen peroxide",
       "2 H2O2(aq) -> 2 H2O(l) + O2(g)",
       "2\\,\\text{H}_2\\text{O}_2(aq) \\rightarrow 2\\,\\text{H}_2\\text{O}(l) + \\text{O}_2(g)",
       "decomposition", "laboratory", ["H", "O"],
       [_reactant("H2O2", 2, "aq")],
       [_reactant("H2O", 2, "l"), _reactant("O2", 1, "g")],
       _thermo(-196, -233.6, 126, True),
       _cond(catalyst="MnO2 (catalase enzyme in biology)"),
       False,
       "Catalytic decomposition of hydrogen peroxide. Elephant toothpaste demo uses concentrated H2O2 + KI."),

    _r("Cl-laboratory-001", "Chlorine generation (MnO2 + HCl)",
       "MnO2(s) + 4 HCl(aq) -> MnCl2(aq) + 2 H2O(l) + Cl2(g)",
       "\\text{MnO}_2(s) + 4\\,\\text{HCl}(aq) \\rightarrow \\text{MnCl}_2(aq) + 2\\,\\text{H}_2\\text{O}(l) + \\text{Cl}_2(g)",
       "oxidation", "laboratory", ["Mn", "Cl", "H", "O"],
       [_reactant("MnO2", 1, "s"), _reactant("HCl", 4, "aq")],
       [_reactant("MnCl2", 1, "aq"), _reactant("H2O", 2, "l"), _reactant("Cl2", 1, "g")],
       _thermo(None, None, None, True),
       _cond(other="use concentrated HCl; fume hood required"),
       False,
       "Laboratory chlorine preparation. Scheele's original 1774 method."),

    _r("N-laboratory-001", "Ammonia from ammonium salt + base",
       "NH4Cl(s) + NaOH(s) -> NaCl(s) + H2O(l) + NH3(g)",
       "\\text{NH}_4\\text{Cl}(s) + \\text{NaOH}(s) \\rightarrow \\text{NaCl}(s) + \\text{H}_2\\text{O}(l) + \\text{NH}_3(g)",
       "base displacement", "laboratory", ["N", "H", "Cl", "Na", "O"],
       [_reactant("NH4Cl", 1, "s"), _reactant("NaOH", 1, "s")],
       [_reactant("NaCl", 1, "s"), _reactant("H2O", 1, "l"), _reactant("NH3", 1, "g")],
       _thermo(None, None, None, False),
       _cond(other="gentle heating"),
       False,
       "Lab preparation of ammonia. Detected by moist red litmus turning blue."),

    _r("S-laboratory-001", "Sulfur dioxide generation",
       "Na2SO3(s) + H2SO4(aq) -> Na2SO4(aq) + H2O(l) + SO2(g)",
       "\\text{Na}_2\\text{SO}_3(s) + \\text{H}_2\\text{SO}_4(aq) \\rightarrow \\text{Na}_2\\text{SO}_4(aq) + \\text{H}_2\\text{O}(l) + \\text{SO}_2(g)",
       "acid-salt", "laboratory", ["Na", "S", "O", "H"],
       [_reactant("Na2SO3", 1, "s"), _reactant("H2SO4", 1, "aq")],
       [_reactant("Na2SO4", 1, "aq"), _reactant("H2O", 1, "l"), _reactant("SO2", 1, "g")],
       _thermo(None, None, None, True),
       _cond(other="fume hood required; SO2 is toxic"),
       False,
       "Lab preparation of SO2 from sulfite salt and acid. Bleaches damp litmus."),

    # --- Organic Synthesis ---
    _r("C-laboratory-002", "Fischer esterification",
       "CH3COOH(l) + C2H5OH(l) <=> CH3COOC2H5(l) + H2O(l)",
       "\\text{CH}_3\\text{COOH}(l) + \\text{C}_2\\text{H}_5\\text{OH}(l) \\rightleftharpoons \\text{CH}_3\\text{COOC}_2\\text{H}_5(l) + \\text{H}_2\\text{O}(l)",
       "esterification", "laboratory", ["C", "H", "O"],
       [_reactant("CH3COOH", 1, "l"), _reactant("C2H5OH", 1, "l")],
       [_reactant("CH3COOC2H5", 1, "l"), _reactant("H2O", 1, "l")],
       _thermo(None, None, None, None),
       _cond(catalyst="concentrated H2SO4", other="reflux; equilibrium reaction"),
       True,
       "Acid-catalyzed esterification. Ethyl acetate has fruity smell. Common teaching reaction."),

    _r("C-laboratory-003", "Saponification (ester hydrolysis)",
       "CH3COOC2H5(l) + NaOH(aq) -> CH3COONa(aq) + C2H5OH(aq)",
       "\\text{CH}_3\\text{COOC}_2\\text{H}_5(l) + \\text{NaOH}(aq) \\rightarrow \\text{CH}_3\\text{COONa}(aq) + \\text{C}_2\\text{H}_5\\text{OH}(aq)",
       "hydrolysis", "laboratory", ["C", "H", "O", "Na"],
       [_reactant("CH3COOC2H5", 1, "l"), _reactant("NaOH", 1, "aq")],
       [_reactant("CH3COONa", 1, "aq"), _reactant("C2H5OH", 1, "aq")],
       _thermo(None, None, None, True),
       _cond(other="heating; irreversible (OH- drives to completion)"),
       False,
       "Base hydrolysis of an ester. Called saponification when applied to fats (soap making)."),

    _r("C-laboratory-004", "Grignard reaction (general)",
       "RMgBr(ether) + R'CHO(l) -> R'CH(OH)R(aq)",
       "\\text{RMgBr}(\\text{ether}) + \\text{R'CHO}(l) \\rightarrow \\text{R'CH(OH)R}(aq)",
       "addition", "laboratory", ["C", "H", "O", "Mg", "Br"],
       [_reactant("RMgBr", 1, "l"), _reactant("RCHO", 1, "l")],
       [_reactant("RCH(OH)R", 1, "l")],
       _thermo(None, None, None, True),
       _cond(other="anhydrous diethyl ether; strictly dry glassware; N2 atmosphere"),
       False,
       "Grignard reagent adds to carbonyl. Produces secondary alcohol from aldehyde. Nobel Prize 1912."),

    _r("Br-laboratory-001", "Bromine water test for unsaturation",
       "C2H4(g) + Br2(aq) -> C2H4Br2(l)",
       "\\text{C}_2\\text{H}_4(g) + \\text{Br}_2(aq) \\rightarrow \\text{C}_2\\text{H}_4\\text{Br}_2(l)",
       "addition", "laboratory", ["C", "H", "Br"],
       [_reactant("C2H4", 1, "g"), _reactant("Br2", 1, "aq")],
       [_reactant("C2H4Br2", 1, "l")],
       _thermo(None, None, None, True),
       _cond(),
       False,
       "Decolorization of bromine water indicates C=C double bond. Classic qualitative test."),

    _r("C-laboratory-005", "Aldol condensation",
       "2 CH3CHO(l) -> CH3CH(OH)CH2CHO(l)",
       "2\\,\\text{CH}_3\\text{CHO}(l) \\rightarrow \\text{CH}_3\\text{CH}(\\text{OH})\\text{CH}_2\\text{CHO}(l)",
       "condensation", "laboratory", ["C", "H", "O"],
       [_reactant("CH3CHO", 2, "l")],
       [_reactant("CH3CH(OH)CH2CHO", 1, "l")],
       _thermo(None, None, None, None),
       _cond(catalyst="NaOH (dilute)", other="low temperature favors aldol; heat gives crotonaldehyde"),
       True,
       "Aldol addition forms beta-hydroxy aldehyde. Cornerstone of organic synthesis."),

    _r("C-laboratory-006", "Cannizzaro reaction (benzaldehyde)",
       "2 C6H5CHO(l) + NaOH(aq) -> C6H5COONa(aq) + C6H5CH2OH(l)",
       "2\\,\\text{C}_6\\text{H}_5\\text{CHO}(l) + \\text{NaOH}(aq) \\rightarrow \\text{C}_6\\text{H}_5\\text{COONa}(aq) + \\text{C}_6\\text{H}_5\\text{CH}_2\\text{OH}(l)",
       "disproportionation", "laboratory", ["C", "H", "O", "Na"],
       [_reactant("C6H5CHO", 2, "l"), _reactant("NaOH", 1, "aq")],
       [_reactant("C6H5COONa", 1, "aq"), _reactant("C6H5CH2OH", 1, "l")],
       _thermo(None, None, None, True),
       _cond(other="concentrated NaOH; aldehydes without alpha-H"),
       False,
       "Non-enolizable aldehyde disproportionation to acid + alcohol. Named for Stanislao Cannizzaro."),

    # --- Electrochemistry ---
    _r("Zn-laboratory-001", "Daniell cell",
       "Zn(s) + Cu2+(aq) -> Zn2+(aq) + Cu(s)",
       "\\text{Zn}(s) + \\text{Cu}^{2+}(aq) \\rightarrow \\text{Zn}^{2+}(aq) + \\text{Cu}(s)",
       "galvanic cell", "laboratory", ["Zn", "Cu"],
       [_reactant("Zn", 1, "s"), _reactant("CuSO4", 1, "aq")],
       [_reactant("ZnSO4", 1, "aq"), _reactant("Cu", 1, "s")],
       _thermo(-212.6, -212.6, 0.05, True),
       _cond(other="E_cell = +1.10 V; salt bridge required"),
       False,
       "Classic galvanic cell. Zn anode oxidizes, Cu2+ deposits at cathode. E = 1.10V."),

    _r("Cu-laboratory-001", "Copper electroplating",
       "Cu2+(aq) + 2 e- -> Cu(s)",
       "\\text{Cu}^{2+}(aq) + 2\\,e^- \\rightarrow \\text{Cu}(s)",
       "electrolysis", "laboratory", ["Cu"],
       [_reactant("CuSO4", 1, "aq")],
       [_reactant("Cu", 1, "s")],
       _thermo(None, None, None, None),
       _cond(other="CuSO4 electrolyte; object as cathode; Cu anode; ~0.5 A/dm2"),
       False,
       "Electrodeposition of copper. Common introductory electrochemistry experiment."),

    # --- Analytical ---
    _r("Fe-laboratory-001", "Iron(III) thiocyanate test",
       "Fe3+(aq) + 3 SCN-(aq) -> Fe(SCN)3(aq)",
       "\\text{Fe}^{3+}(aq) + 3\\,\\text{SCN}^-(aq) \\rightarrow \\text{Fe}(\\text{SCN})_3(aq)",
       "complexation", "laboratory", ["Fe", "S", "C", "N"],
       [_reactant("FeCl3", 1, "aq"), _reactant("KSCN", 3, "aq")],
       [_reactant("Fe(SCN)3", 1, "aq")],
       _thermo(None, None, None, None),
       _cond(),
       False,
       "Blood-red color confirms Fe3+ ions. Very sensitive test. Used in equilibrium demonstrations."),

    _r("Cu-laboratory-002", "Copper(II) ammonia complex (deep blue)",
       "Cu2+(aq) + 4 NH3(aq) -> [Cu(NH3)4]2+(aq)",
       "\\text{Cu}^{2+}(aq) + 4\\,\\text{NH}_3(aq) \\rightarrow [\\text{Cu}(\\text{NH}_3)_4]^{2+}(aq)",
       "complexation", "laboratory", ["Cu", "N", "H"],
       [_reactant("CuSO4", 1, "aq"), _reactant("NH3", 4, "aq")],
       [_reactant("[Cu(NH3)4]SO4", 1, "aq")],
       _thermo(None, None, None, True),
       _cond(other="add excess ammonia to Cu2+ solution"),
       False,
       "Tetraamminecopper(II) complex. Deep royal blue color. Confirms Cu2+ presence."),

    _r("Fe-laboratory-002", "Prussian blue formation",
       "4 Fe3+(aq) + 3 [Fe(CN)6]4-(aq) -> Fe4[Fe(CN)6]3(s)",
       "4\\,\\text{Fe}^{3+}(aq) + 3\\,[\\text{Fe}(\\text{CN})_6]^{4-}(aq) \\rightarrow \\text{Fe}_4[\\text{Fe}(\\text{CN})_6]_3(s)",
       "precipitation", "laboratory", ["Fe", "C", "N"],
       [_reactant("FeCl3", 4, "aq"), _reactant("K4[Fe(CN)6]", 3, "aq")],
       [_reactant("Fe4[Fe(CN)6]3", 1, "s")],
       _thermo(None, None, None, True),
       _cond(),
       False,
       "Prussian blue: intense blue pigment. Historic iron detection test. Also used in art."),

    _r("Na-laboratory-001", "Flame test (sodium)",
       "Na(s) -> Na*(g) -> Na(g) + hv (589 nm)",
       "\\text{Na}(s) \\rightarrow \\text{Na}^*(g) \\rightarrow \\text{Na}(g) + h\\nu\\,(589\\,\\text{nm})",
       "emission", "laboratory", ["Na"],
       [_reactant("NaCl", 1, "s")],
       [_reactant("Na", 1, "g")],
       _thermo(None, None, None, None),
       _cond(other="nichrome wire dipped in conc. HCl then sample; Bunsen burner"),
       False,
       "Intense persistent yellow flame. Sodium D-line at 589 nm. Masks other flame colors."),

    _r("K-laboratory-001", "Flame test (potassium)",
       "K(s) -> K*(g) -> K(g) + hv (766 nm)",
       "\\text{K}(s) \\rightarrow \\text{K}^*(g) \\rightarrow \\text{K}(g) + h\\nu\\,(766\\,\\text{nm})",
       "emission", "laboratory", ["K"],
       [_reactant("KCl", 1, "s")],
       [_reactant("K", 1, "g")],
       _thermo(None, None, None, None),
       _cond(other="view through cobalt blue glass to filter out sodium yellow"),
       False,
       "Lilac/violet flame for potassium. Best observed through blue cobalt glass."),

    _r("Li-laboratory-001", "Flame test (lithium)",
       "Li(s) -> Li*(g) -> Li(g) + hv (671 nm)",
       "\\text{Li}(s) \\rightarrow \\text{Li}^*(g) \\rightarrow \\text{Li}(g) + h\\nu\\,(671\\,\\text{nm})",
       "emission", "laboratory", ["Li"],
       [_reactant("LiCl", 1, "s")],
       [_reactant("Li", 1, "g")],
       _thermo(None, None, None, None),
       _cond(other="nichrome wire; Bunsen burner"),
       False,
       "Crimson red flame for lithium at 671 nm."),

    _r("Sr-laboratory-001", "Flame test (strontium)",
       "Sr(s) -> Sr*(g) -> Sr(g) + hv (606 nm)",
       "\\text{Sr}(s) \\rightarrow \\text{Sr}^*(g) \\rightarrow \\text{Sr}(g) + h\\nu\\,(606\\,\\text{nm})",
       "emission", "laboratory", ["Sr"],
       [_reactant("SrCl2", 1, "s")],
       [_reactant("Sr", 1, "g")],
       _thermo(None, None, None, None),
       _cond(other="nichrome wire; Bunsen burner"),
       False,
       "Bright red flame for strontium. Used in fireworks for red color."),

    _r("Ba-laboratory-002", "Flame test (barium)",
       "Ba(s) -> Ba*(g) -> Ba(g) + hv (524 nm)",
       "\\text{Ba}(s) \\rightarrow \\text{Ba}^*(g) \\rightarrow \\text{Ba}(g) + h\\nu\\,(524\\,\\text{nm})",
       "emission", "laboratory", ["Ba"],
       [_reactant("BaCl2", 1, "s")],
       [_reactant("Ba", 1, "g")],
       _thermo(None, None, None, None),
       _cond(other="nichrome wire; Bunsen burner"),
       False,
       "Apple green flame for barium. Used in fireworks for green color."),

    # --- More Organic ---
    _r("C-laboratory-007", "Tollens' test (silver mirror)",
       "RCHO(aq) + 2 [Ag(NH3)2]+(aq) + 2 OH-(aq) -> RCOO-(aq) + 2 Ag(s) + 4 NH3(aq) + H2O(l)",
       "\\text{RCHO}(aq) + 2\\,[\\text{Ag}(\\text{NH}_3)_2]^+(aq) + 2\\,\\text{OH}^-(aq) \\rightarrow \\text{RCOO}^-(aq) + 2\\,\\text{Ag}(s) + 4\\,\\text{NH}_3(aq) + \\text{H}_2\\text{O}(l)",
       "oxidation", "laboratory", ["C", "H", "O", "Ag", "N"],
       [_reactant("RCHO", 1, "aq"), _reactant("[Ag(NH3)2]OH", 2, "aq")],
       [_reactant("RCOONa", 1, "aq"), _reactant("Ag", 2, "s"), _reactant("NH3", 4, "aq")],
       _thermo(None, None, None, True),
       _cond(other="warm gently; silver deposits as mirror on clean glass"),
       False,
       "Silver mirror test distinguishes aldehydes from ketones. Ag+ reduced to metallic silver."),

    _r("C-laboratory-008", "Fehling's test (reducing sugars)",
       "RCHO(aq) + 2 Cu2+(aq) + 5 OH-(aq) -> RCOO-(aq) + Cu2O(s) + 3 H2O(l)",
       "\\text{RCHO}(aq) + 2\\,\\text{Cu}^{2+}(aq) + 5\\,\\text{OH}^-(aq) \\rightarrow \\text{RCOO}^-(aq) + \\text{Cu}_2\\text{O}(s) + 3\\,\\text{H}_2\\text{O}(l)",
       "oxidation", "laboratory", ["C", "H", "O", "Cu"],
       [_reactant("RCHO", 1, "aq"), _reactant("CuSO4", 2, "aq"), _reactant("NaOH", 5, "aq")],
       [_reactant("RCOONa", 1, "aq"), _reactant("Cu2O", 1, "s"), _reactant("H2O", 3, "l")],
       _thermo(None, None, None, True),
       _cond(other="heat to boiling; blue Cu2+ to brick-red Cu2O precipitate"),
       False,
       "Brick-red Cu2O precipitate indicates reducing sugar or aldehyde. Blue to red color change."),

    # --- Decomposition ---
    _r("K-laboratory-002", "Potassium chlorate decomposition",
       "2 KClO3(s) -> 2 KCl(s) + 3 O2(g)",
       "2\\,\\text{KClO}_3(s) \\rightarrow 2\\,\\text{KCl}(s) + 3\\,\\text{O}_2(g)",
       "decomposition", "laboratory", ["K", "Cl", "O"],
       [_reactant("KClO3", 2, "s")],
       [_reactant("KCl", 2, "s"), _reactant("O2", 3, "g")],
       _thermo(-89.4, -225.4, 456, True),
       _cond(temp_k=673, catalyst="MnO2 lowers temperature to ~523K"),
       False,
       "Thermal decomposition accelerated by MnO2 catalyst. Historically used for O2 preparation."),

    _r("N-laboratory-002", "Ammonium dichromate volcano",
       "(NH4)2Cr2O7(s) -> N2(g) + Cr2O3(s) + 4 H2O(g)",
       "(\\text{NH}_4)_2\\text{Cr}_2\\text{O}_7(s) \\rightarrow \\text{N}_2(g) + \\text{Cr}_2\\text{O}_3(s) + 4\\,\\text{H}_2\\text{O}(g)",
       "decomposition", "laboratory", ["N", "H", "Cr", "O"],
       [_reactant("(NH4)2Cr2O7", 1, "s")],
       [_reactant("N2", 1, "g"), _reactant("Cr2O3", 1, "s"), _reactant("H2O", 4, "g")],
       _thermo(-429, None, None, True),
       _cond(other="ignite with match; spectacular volcanic eruption; toxic Cr(III) product"),
       False,
       "Chemical volcano demonstration. Self-sustaining decomposition once ignited."),
]


# =============================================================================
# BIOLOGICAL REACTIONS
# =============================================================================

BIOLOGICAL_REACTIONS = [
    _r("C-biological-001", "Photosynthesis (overall)",
       "6 CO2(g) + 6 H2O(l) -> C6H12O6(s) + 6 O2(g)",
       "6\\,\\text{CO}_2(g) + 6\\,\\text{H}_2\\text{O}(l) \\rightarrow \\text{C}_6\\text{H}_{12}\\text{O}_6(s) + 6\\,\\text{O}_2(g)",
       "photosynthesis", "biological", ["C", "H", "O"],
       [_reactant("CO2", 6, "g"), _reactant("H2O", 6, "l")],
       [_reactant("C6H12O6", 1, "s"), _reactant("O2", 6, "g")],
       _thermo(2803, 2879, -256, False),
       _cond(other="chlorophyll; light energy (sunlight); thylakoid membranes"),
       False,
       "Converts light energy to chemical energy. Produces all atmospheric O2. Calvin cycle fixes CO2."),

    _r("C-biological-002", "Photosynthesis light reaction (water splitting)",
       "2 H2O(l) -> O2(g) + 4 H+(aq) + 4 e-",
       "2\\,\\text{H}_2\\text{O}(l) \\rightarrow \\text{O}_2(g) + 4\\,\\text{H}^+(aq) + 4\\,e^-",
       "oxidation", "biological", ["H", "O"],
       [_reactant("H2O", 2, "l")],
       [_reactant("O2", 1, "g")],
       _thermo(None, None, None, False),
       _cond(other="photosystem II; Mn4CaO5 cluster (OEC); 680 nm light"),
       False,
       "Water oxidation at PSII oxygen-evolving complex. Source of all biological O2."),

    _r("C-biological-003", "Cellular respiration (overall)",
       "C6H12O6(s) + 6 O2(g) -> 6 CO2(g) + 6 H2O(l)",
       "\\text{C}_6\\text{H}_{12}\\text{O}_6(s) + 6\\,\\text{O}_2(g) \\rightarrow 6\\,\\text{CO}_2(g) + 6\\,\\text{H}_2\\text{O}(l)",
       "respiration", "biological", ["C", "H", "O"],
       [_reactant("C6H12O6", 1, "s"), _reactant("O2", 6, "g")],
       [_reactant("CO2", 6, "g"), _reactant("H2O", 6, "l")],
       _thermo(-2803, -2879, 256, True),
       _cond(other="glycolysis + Krebs cycle + oxidative phosphorylation; ~30-32 ATP per glucose"),
       False,
       "Complete oxidation of glucose. Exact reverse of photosynthesis. Powers all aerobic life."),

    _r("C-biological-004", "Glycolysis (net reaction)",
       "C6H12O6(s) + 2 NAD+(aq) + 2 ADP(aq) + 2 Pi(aq) -> 2 CH3COCOO-(aq) + 2 NADH(aq) + 2 ATP(aq) + 2 H2O(l)",
       "\\text{C}_6\\text{H}_{12}\\text{O}_6 + 2\\,\\text{NAD}^+ + 2\\,\\text{ADP} + 2\\,\\text{P}_i \\rightarrow 2\\,\\text{pyruvate} + 2\\,\\text{NADH} + 2\\,\\text{ATP} + 2\\,\\text{H}_2\\text{O}",
       "glycolysis", "biological", ["C", "H", "O", "N", "P"],
       [_reactant("C6H12O6", 1, "s"), _reactant("NAD+", 2, "aq"),
        _reactant("ADP", 2, "aq"), _reactant("Pi", 2, "aq")],
       [_reactant("pyruvate", 2, "aq"), _reactant("NADH", 2, "aq"),
        _reactant("ATP", 2, "aq"), _reactant("H2O", 2, "l")],
       _thermo(-74, None, None, True),
       _cond(other="cytoplasm; 10 enzyme-catalyzed steps"),
       False,
       "Glucose to 2 pyruvate. Net 2 ATP + 2 NADH. Universal pathway in all domains of life."),

    _r("P-biological-001", "ATP hydrolysis",
       "ATP(aq) + H2O(l) -> ADP(aq) + Pi(aq) + H+(aq)",
       "\\text{ATP}(aq) + \\text{H}_2\\text{O}(l) \\rightarrow \\text{ADP}(aq) + \\text{P}_i(aq) + \\text{H}^+(aq)",
       "hydrolysis", "biological", ["C", "H", "O", "N", "P"],
       [_reactant("ATP", 1, "aq"), _reactant("H2O", 1, "l")],
       [_reactant("ADP", 1, "aq"), _reactant("Pi", 1, "aq")],
       _thermo(-30.5, -30.5, 0, True),
       _cond(other="Mg2+ cofactor; standard biochemical conditions (pH 7, 1mM)"),
       False,
       "Universal energy currency of life. Drives endergonic reactions by coupling. ~100 kg ATP recycled/day in humans."),

    _r("P-biological-002", "ATP synthesis (oxidative phosphorylation)",
       "ADP(aq) + Pi(aq) + H+(aq) -> ATP(aq) + H2O(l)",
       "\\text{ADP}(aq) + \\text{P}_i(aq) + \\text{H}^+(aq) \\rightarrow \\text{ATP}(aq) + \\text{H}_2\\text{O}(l)",
       "phosphorylation", "biological", ["C", "H", "O", "N", "P"],
       [_reactant("ADP", 1, "aq"), _reactant("Pi", 1, "aq")],
       [_reactant("ATP", 1, "aq"), _reactant("H2O", 1, "l")],
       _thermo(30.5, 30.5, 0, False),
       _cond(other="ATP synthase (Complex V); proton motive force across inner mitochondrial membrane"),
       False,
       "Chemiosmotic synthesis of ATP. F1Fo-ATP synthase rotary motor. ~26-28 ATP from one glucose."),

    _r("N-biological-001", "Biological nitrogen fixation",
       "N2(g) + 8 H+(aq) + 8 e- + 16 ATP(aq) -> 2 NH3(aq) + H2(g) + 16 ADP(aq) + 16 Pi(aq)",
       "\\text{N}_2 + 8\\,\\text{H}^+ + 8\\,e^- + 16\\,\\text{ATP} \\rightarrow 2\\,\\text{NH}_3 + \\text{H}_2 + 16\\,\\text{ADP} + 16\\,\\text{P}_i",
       "nitrogen fixation", "biological", ["N", "H", "Fe", "Mo"],
       [_reactant("N2", 1, "g"), _reactant("ATP", 16, "aq")],
       [_reactant("NH3", 2, "aq"), _reactant("H2", 1, "g"), _reactant("ADP", 16, "aq")],
       _thermo(None, None, None, False),
       _cond(catalyst="nitrogenase (Fe-Mo cofactor)", other="obligate anaerobic conditions for enzyme; Rhizobium in legume root nodules"),
       False,
       "Nitrogenase reduces N2 to NH3 at ambient conditions. Requires 16 ATP per N2. FeMo cofactor active site."),

    _r("Fe-biological-001", "Hemoglobin oxygen binding",
       "Hb(aq) + 4 O2(g) <=> Hb(O2)4(aq)",
       "\\text{Hb}(aq) + 4\\,\\text{O}_2(g) \\rightleftharpoons \\text{Hb}(\\text{O}_2)_4(aq)",
       "binding", "biological", ["Fe", "O", "C", "N", "H"],
       [_reactant("Hb", 1, "aq"), _reactant("O2", 4, "g")],
       [_reactant("Hb(O2)4", 1, "aq")],
       _thermo(-67, None, None, True),
       _cond(other="cooperative binding (sigmoidal curve); Fe2+ in heme porphyrin; P50 = 26 mmHg"),
       True,
       "Cooperative O2 binding to hemoglobin's four heme groups. Fe2+ does NOT oxidize to Fe3+."),

    _r("C-biological-005", "Ethanol fermentation",
       "C6H12O6(s) -> 2 C2H5OH(l) + 2 CO2(g)",
       "\\text{C}_6\\text{H}_{12}\\text{O}_6(s) \\rightarrow 2\\,\\text{C}_2\\text{H}_5\\text{OH}(l) + 2\\,\\text{CO}_2(g)",
       "fermentation", "biological", ["C", "H", "O"],
       [_reactant("C6H12O6", 1, "s")],
       [_reactant("C2H5OH", 2, "l"), _reactant("CO2", 2, "g")],
       _thermo(-69, None, None, True),
       _cond(catalyst="zymase (yeast enzymes)", other="anaerobic conditions; Saccharomyces cerevisiae"),
       False,
       "Anaerobic glucose metabolism by yeast. Produces beer, wine, bread (CO2 leavening)."),

    _r("C-biological-006", "Lactic acid fermentation",
       "C6H12O6(s) -> 2 CH3CH(OH)COOH(aq)",
       "\\text{C}_6\\text{H}_{12}\\text{O}_6(s) \\rightarrow 2\\,\\text{CH}_3\\text{CH}(\\text{OH})\\text{COOH}(aq)",
       "fermentation", "biological", ["C", "H", "O"],
       [_reactant("C6H12O6", 1, "s")],
       [_reactant("CH3CH(OH)COOH", 2, "aq")],
       _thermo(-136, None, None, True),
       _cond(catalyst="lactate dehydrogenase", other="anaerobic; muscle cells during intense exercise"),
       False,
       "Anaerobic pathway in muscle. Also used by Lactobacillus in yogurt/sauerkraut production."),

    _r("C-biological-007", "Krebs cycle (overall per acetyl-CoA)",
       "Acetyl-CoA + 3 NAD+ + FAD + GDP + Pi + 2 H2O -> CoA-SH + 2 CO2 + 3 NADH + FADH2 + GTP",
       "\\text{Acetyl-CoA} + 3\\,\\text{NAD}^+ + \\text{FAD} + \\text{GDP} + \\text{P}_i + 2\\,\\text{H}_2\\text{O} \\rightarrow \\text{CoA-SH} + 2\\,\\text{CO}_2 + 3\\,\\text{NADH} + \\text{FADH}_2 + \\text{GTP}",
       "oxidative cycle", "biological", ["C", "H", "O", "N", "S", "P"],
       [_reactant("Acetyl-CoA", 1, "aq"), _reactant("NAD+", 3, "aq"),
        _reactant("FAD", 1, "aq"), _reactant("GDP", 1, "aq"),
        _reactant("Pi", 1, "aq"), _reactant("H2O", 2, "l")],
       [_reactant("CoA-SH", 1, "aq"), _reactant("CO2", 2, "g"),
        _reactant("NADH", 3, "aq"), _reactant("FADH2", 1, "aq"), _reactant("GTP", 1, "aq")],
       _thermo(None, None, None, True),
       _cond(other="mitochondrial matrix; 8 enzyme-catalyzed steps; regulated by [ATP]/[ADP]"),
       False,
       "Central metabolic hub. Oxidizes acetyl groups to CO2, generating electron carriers for ETC."),

    _r("Fe-biological-002", "Cytochrome c oxidase (Complex IV)",
       "4 cytochrome c (Fe2+) + O2 + 8 H+(matrix) -> 4 cytochrome c (Fe3+) + 2 H2O + 4 H+(IMS)",
       "4\\,\\text{cyt c}(\\text{Fe}^{2+}) + \\text{O}_2 + 8\\,\\text{H}^+_{\\text{matrix}} \\rightarrow 4\\,\\text{cyt c}(\\text{Fe}^{3+}) + 2\\,\\text{H}_2\\text{O} + 4\\,\\text{H}^+_{\\text{IMS}}",
       "electron transfer", "biological", ["Fe", "Cu", "O", "H"],
       [_reactant("cytochrome_c_Fe2+", 4, "aq"), _reactant("O2", 1, "g")],
       [_reactant("cytochrome_c_Fe3+", 4, "aq"), _reactant("H2O", 2, "l")],
       _thermo(None, None, None, True),
       _cond(other="inner mitochondrial membrane; contains CuA, CuB, heme a, heme a3"),
       False,
       "Terminal oxidase of electron transport chain. Reduces O2 to H2O. Pumps 4 H+ per O2."),

    _r("N-biological-002", "Urease catalysis",
       "NH2CONH2(aq) + H2O(l) -> 2 NH3(aq) + CO2(g)",
       "\\text{NH}_2\\text{CONH}_2(aq) + \\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{NH}_3(aq) + \\text{CO}_2(g)",
       "enzyme catalysis", "biological", ["N", "H", "C", "O", "Ni"],
       [_reactant("NH2CONH2", 1, "aq"), _reactant("H2O", 1, "l")],
       [_reactant("NH3", 2, "aq"), _reactant("CO2", 1, "g")],
       _thermo(-40, None, None, True),
       _cond(catalyst="urease (Ni2+ active site)", other="rate enhancement: 10^14 over uncatalyzed"),
       False,
       "Urease was the first enzyme crystallized (Sumner, 1926). Contains dinuclear Ni center."),

    _r("Zn-biological-001", "Carbonic anhydrase",
       "CO2(g) + H2O(l) <=> HCO3-(aq) + H+(aq)",
       "\\text{CO}_2(g) + \\text{H}_2\\text{O}(l) \\rightleftharpoons \\text{HCO}_3^-(aq) + \\text{H}^+(aq)",
       "enzyme catalysis", "biological", ["C", "O", "H", "Zn"],
       [_reactant("CO2", 1, "g"), _reactant("H2O", 1, "l")],
       [_reactant("HCO3-", 1, "aq"), _reactant("H+", 1, "aq")],
       _thermo(None, None, None, None),
       _cond(catalyst="carbonic anhydrase (Zn2+ active site)", other="kcat ~10^6 s-1; one of fastest enzymes"),
       True,
       "Zinc metalloenzyme. Catalyzes CO2 hydration ~10^6 times faster than uncatalyzed. Critical for CO2 transport and pH regulation."),

    _r("Mg-biological-001", "Chlorophyll light absorption",
       "Chl + hv -> Chl*",
       "\\text{Chl} + h\\nu \\rightarrow \\text{Chl}^*",
       "photoexcitation", "biological", ["Mg", "C", "N", "H", "O"],
       [_reactant("Chl", 1, "s")],
       [_reactant("Chl*", 1, "s")],
       _thermo(None, None, None, False),
       _cond(other="absorbs red (680 nm) and blue (430 nm); Mg2+ coordinated by porphyrin ring"),
       False,
       "Chlorophyll absorbs light via Mg-porphyrin chromophore. Excited electron enters photosystem."),

    _r("S-biological-001", "Disulfide bond formation (protein folding)",
       "2 R-SH -> R-S-S-R + 2 H+ + 2 e-",
       "2\\,\\text{R-SH} \\rightarrow \\text{R-S-S-R} + 2\\,\\text{H}^+ + 2\\,e^-",
       "oxidation", "biological", ["S", "C", "H", "O", "N"],
       [_reactant("R-SH", 2, "aq")],
       [_reactant("R-S-S-R", 1, "aq")],
       _thermo(None, None, None, None),
       _cond(other="endoplasmic reticulum; protein disulfide isomerase (PDI)"),
       True,
       "Cysteine disulfide bridges stabilize protein tertiary structure. Oxidizing environment of ER."),

    _r("C-biological-008", "Beta-oxidation of fatty acids (per cycle)",
       "Acyl-CoA(Cn) + FAD + NAD+ + CoA-SH + H2O -> Acyl-CoA(Cn-2) + FADH2 + NADH + Acetyl-CoA",
       "\\text{Acyl-CoA}(C_n) + \\text{FAD} + \\text{NAD}^+ + \\text{CoA-SH} + \\text{H}_2\\text{O} \\rightarrow \\text{Acyl-CoA}(C_{n-2}) + \\text{FADH}_2 + \\text{NADH} + \\text{Acetyl-CoA}",
       "oxidation", "biological", ["C", "H", "O", "S", "N"],
       [_reactant("Acyl-CoA", 1, "aq"), _reactant("FAD", 1, "aq"),
        _reactant("NAD+", 1, "aq"), _reactant("CoA-SH", 1, "aq"), _reactant("H2O", 1, "l")],
       [_reactant("Acyl-CoA(n-2)", 1, "aq"), _reactant("FADH2", 1, "aq"),
        _reactant("NADH", 1, "aq"), _reactant("Acetyl-CoA", 1, "aq")],
       _thermo(None, None, None, True),
       _cond(other="mitochondrial matrix; 4 enzymes per cycle; each cycle removes 2 carbons"),
       False,
       "Fatty acid catabolism. Each cycle yields 1 FADH2, 1 NADH, 1 acetyl-CoA (~14 ATP per cycle)."),

    _r("Ca-biological-001", "Hydroxyapatite formation (bone mineralization)",
       "10 Ca2+(aq) + 6 PO4 3-(aq) + 2 OH-(aq) -> Ca10(PO4)6(OH)2(s)",
       "10\\,\\text{Ca}^{2+}(aq) + 6\\,\\text{PO}_4^{3-}(aq) + 2\\,\\text{OH}^-(aq) \\rightarrow \\text{Ca}_{10}(\\text{PO}_4)_6(\\text{OH})_2(s)",
       "mineralization", "biological", ["Ca", "P", "O", "H"],
       [_reactant("Ca2+", 10, "aq"), _reactant("PO4_3-", 6, "aq"), _reactant("OH-", 2, "aq")],
       [_reactant("Ca10(PO4)6(OH)2", 1, "s")],
       _thermo(None, None, None, True),
       _cond(other="osteoblast-mediated; alkaline phosphatase; collagen matrix template"),
       False,
       "Hydroxyapatite is the mineral component of bone and teeth. ~70% of bone mass."),
]


# =============================================================================
# ENVIRONMENTAL REACTIONS
# =============================================================================

ENVIRONMENTAL_REACTIONS = [
    _r("O-environmental-001", "Ozone formation (stratospheric)",
       "3 O2(g) -> 2 O3(g)",
       "3\\,\\text{O}_2(g) \\rightarrow 2\\,\\text{O}_3(g)",
       "photochemical", "environmental", ["O"],
       [_reactant("O2", 3, "g")],
       [_reactant("O3", 2, "g")],
       _thermo(285.4, 326.4, -137.6, False),
       _cond(other="UV radiation (<240 nm); Chapman cycle; stratosphere 20-50 km"),
       False,
       "Stratospheric ozone formation via UV photolysis of O2 then O + O2 -> O3."),

    _r("O-environmental-002", "Ozone depletion by CFCs",
       "O3(g) + Cl(g) -> O2(g) + ClO(g)",
       "\\text{O}_3(g) + \\text{Cl}(g) \\rightarrow \\text{O}_2(g) + \\text{ClO}(g)",
       "catalytic destruction", "environmental", ["O", "Cl"],
       [_reactant("O3", 1, "g"), _reactant("Cl", 1, "g")],
       [_reactant("O2", 1, "g"), _reactant("ClO", 1, "g")],
       _thermo(None, None, None, True),
       _cond(other="Cl radical from CFC photolysis; one Cl atom destroys ~100,000 O3 molecules"),
       False,
       "Catalytic ozone destruction cycle. Cl regenerated: ClO + O -> Cl + O2. Montreal Protocol banned CFCs."),

    _r("N-environmental-001", "Acid rain from NO2",
       "4 NO2(g) + O2(g) + 2 H2O(l) -> 4 HNO3(aq)",
       "4\\,\\text{NO}_2(g) + \\text{O}_2(g) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow 4\\,\\text{HNO}_3(aq)",
       "acid rain", "environmental", ["N", "O", "H"],
       [_reactant("NO2", 4, "g"), _reactant("O2", 1, "g"), _reactant("H2O", 2, "l")],
       [_reactant("HNO3", 4, "aq")],
       _thermo(None, None, None, True),
       _cond(other="atmospheric; NOx from combustion engines and power plants"),
       False,
       "NOx forms nitric acid in atmosphere. Contributes to acid rain (pH < 5.6)."),

    _r("S-environmental-001", "Acid rain from SO2",
       "2 SO2(g) + O2(g) + 2 H2O(l) -> 2 H2SO4(aq)",
       "2\\,\\text{SO}_2(g) + \\text{O}_2(g) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{H}_2\\text{SO}_4(aq)",
       "acid rain", "environmental", ["S", "O", "H"],
       [_reactant("SO2", 2, "g"), _reactant("O2", 1, "g"), _reactant("H2O", 2, "l")],
       [_reactant("H2SO4", 2, "aq")],
       _thermo(None, None, None, True),
       _cond(other="atmospheric oxidation; SO2 from coal/ore smelting"),
       False,
       "Sulfuric acid formation in atmosphere. Major cause of acid rain. Addressed by Clean Air Act."),

    _r("C-environmental-001", "Ocean acidification",
       "CO2(g) + H2O(l) <=> H2CO3(aq) <=> HCO3-(aq) + H+(aq)",
       "\\text{CO}_2(g) + \\text{H}_2\\text{O}(l) \\rightleftharpoons \\text{H}_2\\text{CO}_3(aq) \\rightleftharpoons \\text{HCO}_3^-(aq) + \\text{H}^+(aq)",
       "acidification", "environmental", ["C", "O", "H"],
       [_reactant("CO2", 1, "g"), _reactant("H2O", 1, "l")],
       [_reactant("HCO3-", 1, "aq"), _reactant("H+", 1, "aq")],
       _thermo(None, None, None, None),
       _cond(other="oceans absorb ~30% of anthropogenic CO2; pH dropped from 8.2 to 8.1 since pre-industrial"),
       True,
       "CO2 dissolution lowers ocean pH. Threatens coral reefs and shellfish (CaCO3 dissolution)."),

    _r("Ca-environmental-001", "Chemical weathering of limestone",
       "CaCO3(s) + CO2(g) + H2O(l) -> Ca(HCO3)2(aq)",
       "\\text{CaCO}_3(s) + \\text{CO}_2(g) + \\text{H}_2\\text{O}(l) \\rightarrow \\text{Ca}(\\text{HCO}_3)_2(aq)",
       "weathering", "environmental", ["Ca", "C", "O", "H"],
       [_reactant("CaCO3", 1, "s"), _reactant("CO2", 1, "g"), _reactant("H2O", 1, "l")],
       [_reactant("Ca(HCO3)2", 1, "aq")],
       _thermo(None, None, None, False),
       _cond(other="rain dissolves atmospheric CO2 forming weak carbonic acid"),
       True,
       "Creates caves, karst topography. Major CO2 sink on geological timescale. Causes hard water."),

    _r("Si-environmental-001", "Silicate weathering (CO2 sink)",
       "CaSiO3(s) + 2 CO2(g) + 3 H2O(l) -> Ca2+(aq) + 2 HCO3-(aq) + H4SiO4(aq)",
       "\\text{CaSiO}_3(s) + 2\\,\\text{CO}_2(g) + 3\\,\\text{H}_2\\text{O}(l) \\rightarrow \\text{Ca}^{2+}(aq) + 2\\,\\text{HCO}_3^-(aq) + \\text{H}_4\\text{SiO}_4(aq)",
       "weathering", "environmental", ["Ca", "Si", "C", "O", "H"],
       [_reactant("CaSiO3", 1, "s"), _reactant("CO2", 2, "g"), _reactant("H2O", 3, "l")],
       [_reactant("Ca2+", 1, "aq"), _reactant("HCO3-", 2, "aq"), _reactant("H4SiO4", 1, "aq")],
       _thermo(None, None, None, False),
       _cond(other="geological timescale; thermostat of Earth's climate"),
       False,
       "Silicate weathering is the primary geological CO2 sink. Urey reaction controls climate over millions of years."),

    _r("Fe-environmental-001", "Iron rusting",
       "4 Fe(s) + 3 O2(g) + 6 H2O(l) -> 4 Fe(OH)3(s)",
       "4\\,\\text{Fe}(s) + 3\\,\\text{O}_2(g) + 6\\,\\text{H}_2\\text{O}(l) \\rightarrow 4\\,\\text{Fe}(\\text{OH})_3(s)",
       "corrosion", "environmental", ["Fe", "O", "H"],
       [_reactant("Fe", 4, "s"), _reactant("O2", 3, "g"), _reactant("H2O", 6, "l")],
       [_reactant("Fe(OH)3", 4, "s")],
       _thermo(-1648, None, None, True),
       _cond(other="electrochemical process; accelerated by salt, acid, scratches in protective coating"),
       False,
       "Iron corrosion is an electrochemical process. Costs ~3.4% of global GDP annually."),

    _r("Fe-environmental-002", "Rust dehydration",
       "2 Fe(OH)3(s) -> Fe2O3 * H2O(s) + 2 H2O(l)",
       "2\\,\\text{Fe}(\\text{OH})_3(s) \\rightarrow \\text{Fe}_2\\text{O}_3 \\cdot \\text{H}_2\\text{O}(s) + 2\\,\\text{H}_2\\text{O}(l)",
       "dehydration", "environmental", ["Fe", "O", "H"],
       [_reactant("Fe(OH)3", 2, "s")],
       [_reactant("Fe2O3*H2O", 1, "s"), _reactant("H2O", 2, "l")],
       _thermo(None, None, None, True),
       _cond(other="ambient conditions; Fe(OH)3 dehydrates to reddish-brown rust"),
       False,
       "Hydrated iron(III) oxide (rust). Flaky, non-protective. Unlike Al2O3, does not passivate."),

    _r("N-environmental-002", "Nitrogen fixation by lightning",
       "N2(g) + O2(g) -> 2 NO(g)",
       "\\text{N}_2(g) + \\text{O}_2(g) \\rightarrow 2\\,\\text{NO}(g)",
       "fixation", "environmental", ["N", "O"],
       [_reactant("N2", 1, "g"), _reactant("O2", 1, "g")],
       [_reactant("NO", 2, "g")],
       _thermo(180.6, 173.2, 24.8, False),
       _cond(temp_k=3000, other="lightning provides extreme local temperature; ~5-8 Tg N/year globally"),
       False,
       "Lightning fixes ~5-8 Tg N/year. High activation energy overcome by >30,000 K plasma temperature."),

    _r("N-environmental-003", "Denitrification",
       "2 NO3-(aq) + 10 e- + 12 H+(aq) -> N2(g) + 6 H2O(l)",
       "2\\,\\text{NO}_3^-(aq) + 10\\,e^- + 12\\,\\text{H}^+(aq) \\rightarrow \\text{N}_2(g) + 6\\,\\text{H}_2\\text{O}(l)",
       "denitrification", "environmental", ["N", "O", "H"],
       [_reactant("NO3-", 2, "aq")],
       [_reactant("N2", 1, "g"), _reactant("H2O", 6, "l")],
       _thermo(None, None, None, True),
       _cond(other="anaerobic bacteria (Pseudomonas); wetlands; wastewater treatment"),
       False,
       "Microbial conversion of nitrate back to N2. Closes the nitrogen cycle. Important in wastewater treatment."),

    _r("N-environmental-004", "Nitrification (ammonia to nitrate)",
       "NH3(aq) + 2 O2(g) -> NO3-(aq) + H2O(l) + H+(aq)",
       "\\text{NH}_3(aq) + 2\\,\\text{O}_2(g) \\rightarrow \\text{NO}_3^-(aq) + \\text{H}_2\\text{O}(l) + \\text{H}^+(aq)",
       "nitrification", "environmental", ["N", "H", "O"],
       [_reactant("NH3", 1, "aq"), _reactant("O2", 2, "g")],
       [_reactant("NO3-", 1, "aq"), _reactant("H2O", 1, "l")],
       _thermo(-349, None, None, True),
       _cond(other="two-step: Nitrosomonas (NH3->NO2-) then Nitrobacter (NO2-->NO3-); aerobic soil"),
       False,
       "Two-step microbial oxidation. Converts ammonium fertilizer to plant-available nitrate."),

    _r("C-environmental-002", "Methanogenesis",
       "CO2(g) + 4 H2(g) -> CH4(g) + 2 H2O(l)",
       "\\text{CO}_2(g) + 4\\,\\text{H}_2(g) \\rightarrow \\text{CH}_4(g) + 2\\,\\text{H}_2\\text{O}(l)",
       "methanogenesis", "environmental", ["C", "H", "O"],
       [_reactant("CO2", 1, "g"), _reactant("H2", 4, "g")],
       [_reactant("CH4", 1, "g"), _reactant("H2O", 2, "l")],
       _thermo(-130.7, -130.4, -1.0, True),
       _cond(other="methanogenic archaea; strictly anaerobic; wetlands, ruminant guts, landfills"),
       False,
       "Archaeal methanogenesis. Major biogenic CH4 source. Contributes to greenhouse effect."),

    _r("C-environmental-003", "Methane oxidation (atmospheric sink)",
       "CH4(g) + 2 O2(g) -> CO2(g) + 2 H2O(g)",
       "\\text{CH}_4(g) + 2\\,\\text{O}_2(g) \\rightarrow \\text{CO}_2(g) + 2\\,\\text{H}_2\\text{O}(g)",
       "oxidation", "environmental", ["C", "H", "O"],
       [_reactant("CH4", 1, "g"), _reactant("O2", 2, "g")],
       [_reactant("CO2", 1, "g"), _reactant("H2O", 2, "g")],
       _thermo(-802, None, None, True),
       _cond(other="hydroxyl radical (OH*) initiated; tropospheric lifetime ~12 years"),
       False,
       "Atmospheric methane removal by OH radical. Primary sink for CH4 greenhouse gas."),

    _r("S-environmental-002", "Sulfate reduction (deep ocean)",
       "SO4 2-(aq) + 2 CH2O(s) -> H2S(g) + 2 HCO3-(aq)",
       "\\text{SO}_4^{2-}(aq) + 2\\,\\text{CH}_2\\text{O}(s) \\rightarrow \\text{H}_2\\text{S}(g) + 2\\,\\text{HCO}_3^-(aq)",
       "sulfate reduction", "environmental", ["S", "O", "C", "H"],
       [_reactant("SO4_2-", 1, "aq"), _reactant("CH2O", 2, "s")],
       [_reactant("H2S", 1, "g"), _reactant("HCO3-", 2, "aq")],
       _thermo(None, None, None, True),
       _cond(other="sulfate-reducing bacteria (Desulfovibrio); anaerobic marine sediments"),
       False,
       "Microbial sulfate reduction. Produces H2S (rotten egg smell). Forms black FeS in sediments."),

    _r("O-environmental-003", "Tropospheric ozone formation (smog)",
       "NO2(g) + hv -> NO(g) + O(g); O(g) + O2(g) -> O3(g)",
       "\\text{NO}_2(g) + h\\nu \\rightarrow \\text{NO}(g) + \\text{O}(g);\\quad \\text{O}(g) + \\text{O}_2(g) \\rightarrow \\text{O}_3(g)",
       "photochemical", "environmental", ["N", "O"],
       [_reactant("NO2", 1, "g"), _reactant("O2", 1, "g")],
       [_reactant("NO", 1, "g"), _reactant("O3", 1, "g")],
       _thermo(None, None, None, False),
       _cond(other="UV light; VOCs sustain cycle by converting NO back to NO2 without consuming O3"),
       False,
       "Ground-level ozone is a harmful pollutant (photochemical smog). Good up high, bad nearby."),

    _r("Pb-environmental-001", "Lead pipe corrosion (historical)",
       "Pb(s) + 1/2 O2(g) + H2O(l) -> Pb(OH)2(s)",
       "\\text{Pb}(s) + \\tfrac{1}{2}\\,\\text{O}_2(g) + \\text{H}_2\\text{O}(l) \\rightarrow \\text{Pb}(\\text{OH})_2(s)",
       "corrosion", "environmental", ["Pb", "O", "H"],
       [_reactant("Pb", 1, "s"), _reactant("O2", 0.5, "g"), _reactant("H2O", 1, "l")],
       [_reactant("Pb(OH)2", 1, "s")],
       _thermo(None, None, None, True),
       _cond(other="soft water dissolves Pb(OH)2; hard water forms protective PbCO3 scale"),
       False,
       "Lead leaching from pipes. Flint water crisis (2014). Roman plumbing (plumbum = lead)."),

    _r("Hg-environmental-001", "Mercury methylation (bioaccumulation)",
       "Hg2+(aq) + CH3-(aq) -> CH3Hg+(aq)",
       "\\text{Hg}^{2+}(aq) + \\text{CH}_3^-(aq) \\rightarrow \\text{CH}_3\\text{Hg}^+(aq)",
       "methylation", "environmental", ["Hg", "C", "H"],
       [_reactant("Hg2+", 1, "aq"), _reactant("CH3-", 1, "aq")],
       [_reactant("CH3Hg+", 1, "aq")],
       _thermo(None, None, None, None),
       _cond(other="anaerobic sulfate-reducing bacteria; wetlands and sediments"),
       False,
       "Methylmercury bioaccumulates in food chains. Minamata disease. Extremely neurotoxic."),
]


# =============================================================================
# NOTABLE REACTIONS
# =============================================================================

NOTABLE_REACTIONS = [
    _r("Fe-notable-001", "Thermite reaction",
       "2 Al(s) + Fe2O3(s) -> Al2O3(s) + 2 Fe(l)",
       "2\\,\\text{Al}(s) + \\text{Fe}_2\\text{O}_3(s) \\rightarrow \\text{Al}_2\\text{O}_3(s) + 2\\,\\text{Fe}(l)",
       "aluminothermic", "notable", ["Al", "Fe", "O"],
       [_reactant("Al", 2, "s"), _reactant("Fe2O3", 1, "s")],
       [_reactant("Al2O3", 1, "s"), _reactant("Fe", 2, "l")],
       _thermo(-851.5, -818.4, -111, True),
       _cond(other="ignited with Mg ribbon; reaches ~2500C; used for rail welding"),
       False,
       "Extremely exothermic aluminothermic reaction. Produces molten iron at >2500C. Used to weld railroad rails."),

    _r("Na-notable-001", "Sodium in water",
       "2 Na(s) + 2 H2O(l) -> 2 NaOH(aq) + H2(g)",
       "2\\,\\text{Na}(s) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{NaOH}(aq) + \\text{H}_2(g)",
       "single displacement", "notable", ["Na", "H", "O"],
       [_reactant("Na", 2, "s"), _reactant("H2O", 2, "l")],
       [_reactant("NaOH", 2, "aq"), _reactant("H2", 1, "g")],
       _thermo(-368, None, None, True),
       _cond(other="violent; Na melts and skates on water surface; H2 may ignite with yellow flame"),
       False,
       "Alkali metal reacts violently with water. Na melts into a ball, dashes across surface. H2 may ignite."),

    _r("K-notable-001", "Potassium in water",
       "2 K(s) + 2 H2O(l) -> 2 KOH(aq) + H2(g)",
       "2\\,\\text{K}(s) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{KOH}(aq) + \\text{H}_2(g)",
       "single displacement", "notable", ["K", "H", "O"],
       [_reactant("K", 2, "s"), _reactant("H2O", 2, "l")],
       [_reactant("KOH", 2, "aq"), _reactant("H2", 1, "g")],
       _thermo(-392, None, None, True),
       _cond(other="more violent than Na; lilac flame from K vapor; may explode"),
       False,
       "Even more violent than sodium. H2 always ignites with lilac flame. Can explode."),

    _r("Li-notable-001", "Lithium in water",
       "2 Li(s) + 2 H2O(l) -> 2 LiOH(aq) + H2(g)",
       "2\\,\\text{Li}(s) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{LiOH}(aq) + \\text{H}_2(g)",
       "single displacement", "notable", ["Li", "H", "O"],
       [_reactant("Li", 2, "s"), _reactant("H2O", 2, "l")],
       [_reactant("LiOH", 2, "aq"), _reactant("H2", 1, "g")],
       _thermo(-446, None, None, True),
       _cond(other="gentle fizzing; lithium floats and slowly dissolves; crimson specks"),
       False,
       "Least violent alkali metal + water. Gentle fizzing. Li has highest ionization energy of group 1."),

    _r("Rb-notable-001", "Rubidium in water",
       "2 Rb(s) + 2 H2O(l) -> 2 RbOH(aq) + H2(g)",
       "2\\,\\text{Rb}(s) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{RbOH}(aq) + \\text{H}_2(g)",
       "single displacement", "notable", ["Rb", "H", "O"],
       [_reactant("Rb", 2, "s"), _reactant("H2O", 2, "l")],
       [_reactant("RbOH", 2, "aq"), _reactant("H2", 1, "g")],
       _thermo(None, None, None, True),
       _cond(other="extremely violent; spontaneous ignition of H2; handled only in inert atmosphere"),
       False,
       "Extremely violent reaction. Ignites immediately on contact with water. Must handle under argon."),

    _r("Cs-notable-001", "Caesium in water",
       "2 Cs(s) + 2 H2O(l) -> 2 CsOH(aq) + H2(g)",
       "2\\,\\text{Cs}(s) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow 2\\,\\text{CsOH}(aq) + \\text{H}_2(g)",
       "single displacement", "notable", ["Cs", "H", "O"],
       [_reactant("Cs", 2, "s"), _reactant("H2O", 2, "l")],
       [_reactant("CsOH", 2, "aq"), _reactant("H2", 1, "g")],
       _thermo(None, None, None, True),
       _cond(other="explosive; shatters glass container; Coulomb explosion (2015 Nature Chemistry)"),
       False,
       "Most violent alkali metal + water. Explodes instantly. 2015 study showed Coulomb explosion mechanism."),

    _r("Xe-notable-001", "Xenon hexafluoroplatinate (first noble gas compound)",
       "Xe(g) + PtF6(g) -> [XePtF6](s)",
       "\\text{Xe}(g) + \\text{PtF}_6(g) \\rightarrow [\\text{XePtF}_6](s)",
       "synthesis", "notable", ["Xe", "Pt", "F"],
       [_reactant("Xe", 1, "g"), _reactant("PtF6", 1, "g")],
       [_reactant("XePtF6", 1, "s")],
       _thermo(None, None, None, True),
       _cond(other="Neil Bartlett, 1962; orange-yellow solid"),
       False,
       "First noble gas compound. Neil Bartlett (1962) shattered the 'noble gases are inert' dogma."),

    _r("Xe-notable-002", "Xenon difluoride synthesis",
       "Xe(g) + F2(g) -> XeF2(s)",
       "\\text{Xe}(g) + \\text{F}_2(g) \\rightarrow \\text{XeF}_2(s)",
       "synthesis", "notable", ["Xe", "F"],
       [_reactant("Xe", 1, "g"), _reactant("F2", 1, "g")],
       [_reactant("XeF2", 1, "s")],
       _thermo(-108, None, None, True),
       _cond(other="UV irradiation or high pressure; colorless crystals; used in semiconductor etching"),
       False,
       "Simplest noble gas compound made synthetically. Linear molecule. Used in semiconductor fabrication."),

    _r("U-notable-001", "Uranium-235 nuclear fission",
       "U-235 + n -> Ba-141 + Kr-92 + 3 n",
       "{}^{235}\\text{U} + n \\rightarrow {}^{141}\\text{Ba} + {}^{92}\\text{Kr} + 3\\,n",
       "nuclear fission", "notable", ["U", "Ba", "Kr"],
       [_reactant("U-235", 1, "s")],
       [_reactant("Ba-141", 1, "s"), _reactant("Kr-92", 1, "s")],
       _thermo(None, None, None, True),
       _cond(other="thermal neutron capture; ~200 MeV per fission; chain reaction if k >= 1"),
       False,
       "Nuclear fission of U-235. ~200 MeV per event. Powers nuclear reactors and weapons. Discovered 1938 by Hahn & Strassmann."),

    _r("H-notable-001", "Hydrogen fusion (proton-proton chain)",
       "4 H-1 -> He-4 + 2 e+ + 2 ve + energy",
       "4\\,{}^1\\text{H} \\rightarrow {}^4\\text{He} + 2\\,e^+ + 2\\,\\nu_e + \\text{energy}",
       "nuclear fusion", "notable", ["H", "He"],
       [_reactant("H-1", 4, "g")],
       [_reactant("He-4", 1, "g")],
       _thermo(None, None, None, True),
       _cond(temp_k=15000000, other="stellar core; requires ~15 million K to overcome Coulomb barrier"),
       False,
       "Powers the Sun and main-sequence stars. Converts 4 protons to helium-4. E = mc^2: 0.7% mass converted to energy."),

    _r("H-notable-002", "Deuterium-tritium fusion (ITER target)",
       "D + T -> He-4 + n + 17.6 MeV",
       "{}^2\\text{H} + {}^3\\text{H} \\rightarrow {}^4\\text{He} + n + 17.6\\,\\text{MeV}",
       "nuclear fusion", "notable", ["H", "He"],
       [_reactant("D", 1, "g"), _reactant("T", 1, "g")],
       [_reactant("He-4", 1, "g")],
       _thermo(None, None, None, True),
       _cond(temp_k=100000000, other="tokamak magnetic confinement; ITER target: Q >= 10"),
       False,
       "Easiest fusion reaction to achieve. 17.6 MeV per event. ITER aims for net energy gain."),

    _r("N-notable-001", "Nitrogen triiodide detonation",
       "2 NI3(s) -> N2(g) + 3 I2(g)",
       "2\\,\\text{NI}_3(s) \\rightarrow \\text{N}_2(g) + 3\\,\\text{I}_2(g)",
       "decomposition", "notable", ["N", "I"],
       [_reactant("NI3", 2, "s")],
       [_reactant("N2", 1, "g"), _reactant("I2", 3, "g")],
       _thermo(-290, None, None, True),
       _cond(other="contact-sensitive when dry; detonates from slightest touch; purple I2 cloud"),
       False,
       "Extremely sensitive contact explosive. Detonates when dry from a feather touch. Purple iodine cloud produced."),

    _r("Mg-notable-001", "Magnesium ribbon burning",
       "2 Mg(s) + O2(g) -> 2 MgO(s)",
       "2\\,\\text{Mg}(s) + \\text{O}_2(g) \\rightarrow 2\\,\\text{MgO}(s)",
       "combustion", "notable", ["Mg", "O"],
       [_reactant("Mg", 2, "s"), _reactant("O2", 1, "g")],
       [_reactant("MgO", 2, "s")],
       _thermo(-1204, -1138, -216.6, True),
       _cond(other="brilliant white light; do not look directly; burns at ~3100C; cannot be extinguished with water"),
       False,
       "Brilliant white flame used in fireworks and flares. Burns hot enough to ignite in CO2 and N2."),

    _r("Mg-notable-002", "Magnesium in dry ice",
       "2 Mg(s) + CO2(g) -> 2 MgO(s) + C(s)",
       "2\\,\\text{Mg}(s) + \\text{CO}_2(g) \\rightarrow 2\\,\\text{MgO}(s) + \\text{C}(s)",
       "combustion", "notable", ["Mg", "C", "O"],
       [_reactant("Mg", 2, "s"), _reactant("CO2", 1, "g")],
       [_reactant("MgO", 2, "s"), _reactant("C", 1, "s")],
       _thermo(None, None, None, True),
       _cond(other="Mg ribbon burns inside CO2; demonstrates that CO2 fire extinguisher fails for Mg fires"),
       False,
       "Magnesium burns in CO2, reducing it to carbon. Demonstrates that CO2 extinguishers are useless for metal fires."),

    _r("P-notable-001", "White phosphorus spontaneous combustion",
       "P4(s) + 5 O2(g) -> P4O10(s)",
       "\\text{P}_4(s) + 5\\,\\text{O}_2(g) \\rightarrow \\text{P}_4\\text{O}_{10}(s)",
       "combustion", "notable", ["P", "O"],
       [_reactant("P4", 1, "s"), _reactant("O2", 5, "g")],
       [_reactant("P4O10", 1, "s")],
       _thermo(-2984, None, None, True),
       _cond(other="pyrophoric; ignites spontaneously in air above ~30C; stored under water"),
       False,
       "White phosphorus is pyrophoric. Ignites spontaneously in air. Must be stored under water. Extremely toxic."),

    _r("F-notable-001", "Fluorine + water",
       "2 F2(g) + 2 H2O(l) -> 4 HF(aq) + O2(g)",
       "2\\,\\text{F}_2(g) + 2\\,\\text{H}_2\\text{O}(l) \\rightarrow 4\\,\\text{HF}(aq) + \\text{O}_2(g)",
       "oxidation", "notable", ["F", "H", "O"],
       [_reactant("F2", 2, "g"), _reactant("H2O", 2, "l")],
       [_reactant("HF", 4, "aq"), _reactant("O2", 1, "g")],
       _thermo(-758, None, None, True),
       _cond(other="fluorine oxidizes water (unique among halogens); extremely dangerous"),
       False,
       "Fluorine is the only element that oxidizes water at room temperature. Most reactive element."),

    _r("Au-notable-001", "Gold dissolution in aqua regia",
       "Au(s) + 3 HCl(aq) + HNO3(aq) -> HAuCl4(aq) + NO(g) + 2 H2O(l)",
       "\\text{Au}(s) + 3\\,\\text{HCl}(aq) + \\text{HNO}_3(aq) \\rightarrow \\text{HAuCl}_4(aq) + \\text{NO}(g) + 2\\,\\text{H}_2\\text{O}(l)",
       "dissolution", "notable", ["Au", "H", "Cl", "N", "O"],
       [_reactant("Au", 1, "s"), _reactant("HCl", 3, "aq"), _reactant("HNO3", 1, "aq")],
       [_reactant("HAuCl4", 1, "aq"), _reactant("NO", 1, "g"), _reactant("H2O", 2, "l")],
       _thermo(None, None, None, True),
       _cond(other="3:1 HCl:HNO3 mixture; chloride complexation drives dissolution"),
       False,
       "Aqua regia dissolves gold. Neither acid alone works. Cl- complexation (AuCl4-) shifts equilibrium. de Hevesy saved Nobel medals this way in WWII."),

    _r("Pt-notable-001", "Platinum dissolution in aqua regia",
       "3 Pt(s) + 4 HNO3(aq) + 12 HCl(aq) -> 3 H2PtCl4(aq) + 4 NO(g) + 8 H2O(l)",
       "3\\,\\text{Pt}(s) + 4\\,\\text{HNO}_3(aq) + 12\\,\\text{HCl}(aq) \\rightarrow 3\\,\\text{H}_2\\text{PtCl}_4(aq) + 4\\,\\text{NO}(g) + 8\\,\\text{H}_2\\text{O}(l)",
       "dissolution", "notable", ["Pt", "H", "N", "O", "Cl"],
       [_reactant("Pt", 3, "s"), _reactant("HNO3", 4, "aq"), _reactant("HCl", 12, "aq")],
       [_reactant("H2PtCl4", 3, "aq"), _reactant("NO", 4, "g"), _reactant("H2O", 8, "l")],
       _thermo(None, None, None, True),
       _cond(other="aqua regia ('royal water') named for dissolving noble metals"),
       False,
       "Aqua regia dissolves platinum. Named 'royal water' by alchemists. Only method for dissolving Pt."),

    _r("Hg-notable-001", "Mercury(II) oxide decomposition (Priestley's experiment)",
       "2 HgO(s) -> 2 Hg(l) + O2(g)",
       "2\\,\\text{HgO}(s) \\rightarrow 2\\,\\text{Hg}(l) + \\text{O}_2(g)",
       "decomposition", "notable", ["Hg", "O"],
       [_reactant("HgO", 2, "s")],
       [_reactant("Hg", 2, "l"), _reactant("O2", 1, "g")],
       _thermo(181.6, 118.2, 213, False),
       _cond(temp_k=773, other="Priestley (1774) and Scheele independently discovered O2 this way"),
       False,
       "How oxygen was discovered. Priestley heated HgO with a lens and collected 'dephlogisticated air'."),

    _r("Cu-notable-001", "Copper sulfate + iron (tree of Diana)",
       "Fe(s) + CuSO4(aq) -> FeSO4(aq) + Cu(s)",
       "\\text{Fe}(s) + \\text{CuSO}_4(aq) \\rightarrow \\text{FeSO}_4(aq) + \\text{Cu}(s)",
       "single displacement", "notable", ["Fe", "Cu", "S", "O"],
       [_reactant("Fe", 1, "s"), _reactant("CuSO4", 1, "aq")],
       [_reactant("FeSO4", 1, "aq"), _reactant("Cu", 1, "s")],
       _thermo(-153, None, None, True),
       _cond(other="iron nail in CuSO4 solution; copper plates out as reddish deposit"),
       False,
       "Iron displaces copper from solution. Classic activity series demonstration. Reddish copper deposits on iron."),

    _r("C-notable-001", "Elephant toothpaste",
       "2 H2O2(aq) -> 2 H2O(l) + O2(g)",
       "2\\,\\text{H}_2\\text{O}_2(aq) \\rightarrow 2\\,\\text{H}_2\\text{O}(l) + \\text{O}_2(g)",
       "decomposition", "notable", ["H", "O"],
       [_reactant("H2O2", 2, "aq")],
       [_reactant("H2O", 2, "l"), _reactant("O2", 1, "g")],
       _thermo(-196, None, None, True),
       _cond(catalyst="KI or MnO2", other="30% H2O2 + dish soap + KI; massive foam eruption"),
       False,
       "Rapid catalytic decomposition of concentrated H2O2. O2 trapped in soap creates massive foam column."),

    _r("C-notable-002", "Barking dog reaction",
       "2 CS2(l) + 6 N2O(g) -> 2 CO2(g) + 6 N2(g) + 4 SO2(g)",
       "2\\,\\text{CS}_2(l) + 6\\,\\text{N}_2\\text{O}(g) \\rightarrow 2\\,\\text{CO}_2(g) + 6\\,\\text{N}_2(g) + 4\\,\\text{SO}_2(g)",
       "combustion", "notable", ["C", "S", "N", "O"],
       [_reactant("CS2", 2, "l"), _reactant("N2O", 6, "g")],
       [_reactant("CO2", 2, "g"), _reactant("N2", 6, "g"), _reactant("SO2", 4, "g")],
       _thermo(None, None, None, True),
       _cond(other="ignite CS2/N2O mixture in long tube; blue flame + barking sound"),
       False,
       "CS2 vapor and N2O combust with distinctive barking sound and blue flash. Popular demonstration."),

    _r("Mn-notable-001", "Briggs-Rauscher oscillating reaction",
       "IO3-(aq) + 2 H2O2(aq) + CH2(COOH)2(aq) + H+(aq) -> ICH(COOH)2(aq) + 2 O2(g) + 3 H2O(l)",
       "\\text{IO}_3^-(aq) + 2\\,\\text{H}_2\\text{O}_2(aq) + \\text{CH}_2(\\text{COOH})_2(aq) + \\text{H}^+(aq) \\rightarrow \\text{ICH}(\\text{COOH})_2(aq) + 2\\,\\text{O}_2(g) + 3\\,\\text{H}_2\\text{O}(l)",
       "oscillating", "notable", ["I", "O", "H", "C", "Mn"],
       [_reactant("KIO3", 1, "aq"), _reactant("H2O2", 2, "aq"),
        _reactant("CH2(COOH)2", 1, "aq"), _reactant("H2SO4", 1, "aq")],
       [_reactant("ICH(COOH)2", 1, "aq"), _reactant("O2", 2, "g"), _reactant("H2O", 3, "l")],
       _thermo(None, None, None, True),
       _cond(catalyst="MnSO4", other="oscillates amber -> blue -> colorless for ~10 cycles"),
       False,
       "Chemical oscillator. Color cycles amber -> blue -> colorless repeatedly. Non-equilibrium thermodynamics."),

    _r("N-notable-002", "Ammonium nitrate decomposition (explosive)",
       "2 NH4NO3(s) -> 2 N2(g) + O2(g) + 4 H2O(g)",
       "2\\,\\text{NH}_4\\text{NO}_3(s) \\rightarrow 2\\,\\text{N}_2(g) + \\text{O}_2(g) + 4\\,\\text{H}_2\\text{O}(g)",
       "decomposition", "notable", ["N", "H", "O"],
       [_reactant("NH4NO3", 2, "s")],
       [_reactant("N2", 2, "g"), _reactant("O2", 1, "g"), _reactant("H2O", 4, "g")],
       _thermo(-236, None, None, True),
       _cond(other="detonation >270C; Beirut 2020 (2750 tonnes); Texas City 1947; Oppau 1921"),
       False,
       "Ammonium nitrate detonation. Responsible for numerous industrial disasters. Also a common fertilizer."),

    _r("W-notable-001", "Tungsten highest melting point metal",
       "WO3(s) + 3 H2(g) -> W(s) + 3 H2O(g)",
       "\\text{WO}_3(s) + 3\\,\\text{H}_2(g) \\rightarrow \\text{W}(s) + 3\\,\\text{H}_2\\text{O}(g)",
       "reduction", "notable", ["W", "O", "H"],
       [_reactant("WO3", 1, "s"), _reactant("H2", 3, "g")],
       [_reactant("W", 1, "s"), _reactant("H2O", 3, "g")],
       _thermo(117.3, None, None, False),
       _cond(temp_k=1073, other="hydrogen reduction of tungsten trioxide; used for light bulb filaments"),
       False,
       "Production of tungsten metal. Highest melting point of any element (3695K). Used in light bulb filaments."),

    _r("B-notable-001", "Boron green flame",
       "4 B(s) + 3 O2(g) -> 2 B2O3(s)",
       "4\\,\\text{B}(s) + 3\\,\\text{O}_2(g) \\rightarrow 2\\,\\text{B}_2\\text{O}_3(s)",
       "combustion", "notable", ["B", "O"],
       [_reactant("B", 4, "s"), _reactant("O2", 3, "g")],
       [_reactant("B2O3", 2, "s")],
       _thermo(-2546, None, None, True),
       _cond(other="brilliant green flame; boron compounds used in fireworks for green color"),
       False,
       "Boron burns with a distinctive bright green flame. Boric acid used in firework coloring."),

    _r("Ga-notable-001", "Gallium melts in hand",
       "Ga(s) -> Ga(l)",
       "\\text{Ga}(s) \\rightarrow \\text{Ga}(l)",
       "phase change", "notable", ["Ga"],
       [_reactant("Ga", 1, "s")],
       [_reactant("Ga", 1, "l")],
       _thermo(5.59, None, None, False),
       _cond(temp_k=302.91, other="melting point 29.76C; melts from body heat; 'disappearing spoon' trick"),
       False,
       "Gallium melts just above room temperature (29.76C). 'Disappearing spoon' demonstration."),

    _r("Pu-notable-001", "Plutonium-239 fission",
       "Pu-239 + n -> Ba-146 + Sr-91 + 3 n",
       "{}^{239}\\text{Pu} + n \\rightarrow {}^{146}\\text{Ba} + {}^{91}\\text{Sr} + 3\\,n",
       "nuclear fission", "notable", ["Pu", "Ba", "Sr"],
       [_reactant("Pu-239", 1, "s")],
       [_reactant("Ba-146", 1, "s"), _reactant("Sr-91", 1, "s")],
       _thermo(None, None, None, True),
       _cond(other="thermal neutron capture; ~210 MeV per fission; bred from U-238 in reactors"),
       False,
       "Plutonium fission. Pu-239 bred from U-238 by neutron capture. Used in breeder reactors and weapons."),

    _r("Bi-notable-001", "Bismuth subsalicylate (Pepto-Bismol pink)",
       "Bi2O3(s) + 6 HCl(aq) -> 2 BiCl3(aq) + 3 H2O(l)",
       "\\text{Bi}_2\\text{O}_3(s) + 6\\,\\text{HCl}(aq) \\rightarrow 2\\,\\text{BiCl}_3(aq) + 3\\,\\text{H}_2\\text{O}(l)",
       "acid dissolution", "notable", ["Bi", "O", "H", "Cl"],
       [_reactant("Bi2O3", 1, "s"), _reactant("HCl", 6, "aq")],
       [_reactant("BiCl3", 2, "aq"), _reactant("H2O", 3, "l")],
       _thermo(None, None, None, True),
       _cond(),
       False,
       "Bismuth is the heaviest non-radioactive element used in medicine. Pepto-Bismol's pink active ingredient."),

    _r("Os-notable-001", "Osmium tetroxide formation",
       "Os(s) + 2 O2(g) -> OsO4(s)",
       "\\text{Os}(s) + 2\\,\\text{O}_2(g) \\rightarrow \\text{OsO}_4(s)",
       "oxidation", "notable", ["Os", "O"],
       [_reactant("Os", 1, "s"), _reactant("O2", 2, "g")],
       [_reactant("OsO4", 1, "s")],
       _thermo(None, None, None, True),
       _cond(other="OsO4 is volatile, extremely toxic; powerful oxidizer used in organic chemistry for dihydroxylation"),
       False,
       "Osmium forms volatile, highly toxic OsO4. Used in organic synthesis (Sharpless dihydroxylation, Nobel 2001). Densest element."),
]


# =============================================================================
# GENERATION LOGIC
# =============================================================================

ALL_CATEGORIES = {
    "industrial": INDUSTRIAL_REACTIONS,
    "laboratory": LABORATORY_REACTIONS,
    "biological": BIOLOGICAL_REACTIONS,
    "environmental": ENVIRONMENTAL_REACTIONS,
    "notable": NOTABLE_REACTIONS,
}


def build_element_file_map():
    """Build a map of element symbol -> file path from elements/ directory."""
    file_map = {}
    if not ELEMENTS_DIR.exists():
        print(f"ERROR: Elements directory not found: {ELEMENTS_DIR}", file=sys.stderr)
        sys.exit(1)
    for filepath in sorted(ELEMENTS_DIR.glob("*.json")):
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                data = json.load(f)
            symbol = data.get("symbol")
            if symbol:
                file_map[symbol] = filepath
        except (json.JSONDecodeError, KeyError) as e:
            print(f"WARNING: Could not parse {filepath.name}: {e}", file=sys.stderr)
    return file_map


def build_element_reactions_map(all_reactions):
    """Map each element symbol to a list of simplified reaction dicts for injection."""
    element_map = {}  # symbol -> list of simplified reaction dicts
    for rxn in all_reactions:
        simplified = {
            "id": rxn["id"],
            "name": rxn["name"],
            "equation": rxn["equation"],
            "type": rxn["type"],
            "category": rxn["category"],
            "delta_h_kj": rxn["thermodynamics"]["delta_h_kj"],
            "conditions": _summarize_conditions(rxn["conditions"]),
            "description": rxn["description"],
            "reversible": rxn["reversible"],
        }
        for symbol in rxn["elements_involved"]:
            if symbol not in element_map:
                element_map[symbol] = []
            element_map[symbol].append(simplified)
    return element_map


def _summarize_conditions(conditions):
    """Create a human-readable condition string from conditions dict."""
    parts = []
    if conditions.get("temperature_k") is not None:
        parts.append(f"{conditions['temperature_k']} K")
    if conditions.get("pressure_atm") is not None:
        parts.append(f"{conditions['pressure_atm']} atm")
    if conditions.get("catalyst"):
        parts.append(conditions["catalyst"])
    if conditions.get("other"):
        parts.append(conditions["other"])
    return "; ".join(parts) if parts else None


def build_index(all_reactions):
    """Build master index mapping reaction IDs to categories and elements."""
    index = {}
    for rxn in all_reactions:
        index[rxn["id"]] = {
            "name": rxn["name"],
            "category": rxn["category"],
            "type": rxn["type"],
            "elements_involved": rxn["elements_involved"],
            "equation": rxn["equation"],
        }
    return index


def validate_reactions(all_reactions):
    """Validate reaction data for common issues."""
    seen_ids = set()
    errors = []
    for rxn in all_reactions:
        rxn_id = rxn["id"]
        # Check for duplicate IDs
        if rxn_id in seen_ids:
            errors.append(f"Duplicate reaction ID: {rxn_id}")
        seen_ids.add(rxn_id)

        # Check elements_involved are valid symbols
        for sym in rxn["elements_involved"]:
            if sym not in SYMBOL_TO_ELEMENT:
                errors.append(f"{rxn_id}: Unknown element symbol '{sym}'")

        # Check ID format matches pattern: Symbol-category-NNN
        parts = rxn_id.split("-")
        if len(parts) != 3:
            errors.append(f"{rxn_id}: ID should have format Symbol-category-NNN")
        else:
            symbol_part, cat_part, num_part = parts
            if symbol_part not in SYMBOL_TO_ELEMENT:
                errors.append(f"{rxn_id}: ID symbol '{symbol_part}' not a valid element")
            if cat_part != rxn["category"]:
                errors.append(f"{rxn_id}: ID category '{cat_part}' doesn't match category '{rxn['category']}'")

    return errors


def write_json(filepath, data):
    """Write JSON to file with consistent formatting."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")


def main():
    """Generate all reaction files and inject into element files."""
    print("=" * 70)
    print("Chemical Reactions Generator")
    print("=" * 70)

    # Collect all reactions
    all_reactions = []
    for category, reactions in ALL_CATEGORIES.items():
        all_reactions.extend(reactions)
    print(f"\nTotal reactions: {len(all_reactions)}")
    for cat, rxns in ALL_CATEGORIES.items():
        print(f"  {cat:15s}: {len(rxns):3d} reactions")

    # Validate
    print("\nValidating reactions...")
    errors = validate_reactions(all_reactions)
    if errors:
        print(f"\nFOUND {len(errors)} VALIDATION ERROR(S):", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        sys.exit(1)
    print("  All reactions valid.")

    # Write category files
    print("\nWriting category files...")
    for category, reactions in ALL_CATEGORIES.items():
        filepath = REACTIONS_DIR / f"{category}.json"
        write_json(filepath, reactions)
        print(f"  {filepath.name:25s} ({len(reactions)} reactions)")

    # Write master index
    print("\nWriting index file...")
    index = build_index(all_reactions)
    index_path = REACTIONS_DIR / "index.json"
    write_json(index_path, index)
    print(f"  {index_path.name:25s} ({len(index)} entries)")

    # Count unique elements covered
    all_elements = set()
    for rxn in all_reactions:
        all_elements.update(rxn["elements_involved"])
    print(f"\nUnique elements covered by reactions: {len(all_elements)} / 118")

    # Build element -> reactions map and inject into element files
    print("\nInjecting reactions into element files...")
    element_file_map = build_element_file_map()
    element_reactions_map = build_element_reactions_map(all_reactions)

    injected_count = 0
    skipped_count = 0
    for symbol in sorted(element_reactions_map.keys()):
        if symbol not in element_file_map:
            print(f"  WARNING: No element file found for {symbol}", file=sys.stderr)
            skipped_count += 1
            continue
        filepath = element_file_map[symbol]
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)
        # Sort reactions by category then ID for consistent ordering
        reactions = sorted(element_reactions_map[symbol],
                           key=lambda r: (r["category"], r["id"]))
        data["reactions"] = reactions
        write_json(filepath, data)
        injected_count += 1

    # Also clear reactions for elements that have no reactions
    for symbol, filepath in element_file_map.items():
        if symbol not in element_reactions_map:
            with open(filepath, "r", encoding="utf-8") as f:
                data = json.load(f)
            if data.get("reactions") != []:
                data["reactions"] = []
                write_json(filepath, data)

    print(f"  Updated {injected_count} element files with reactions")
    if skipped_count:
        print(f"  Skipped {skipped_count} elements (no file found)")

    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Category files written:    {len(ALL_CATEGORIES)}")
    print(f"  Index entries:             {len(index)}")
    print(f"  Total unique reactions:    {len(all_reactions)}")
    print(f"  Elements with reactions:   {injected_count}")
    print(f"  Elements without reactions: {118 - injected_count}")

    # Print element coverage details
    print(f"\nElements covered:")
    for symbol in sorted(all_elements):
        z, name = SYMBOL_TO_ELEMENT[symbol]
        count = len(element_reactions_map.get(symbol, []))
        print(f"  {z:3d} {symbol:3s} ({name:15s}): {count} reactions")

    print("\nDone.")


if __name__ == "__main__":
    main()
