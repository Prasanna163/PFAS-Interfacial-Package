from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

# Prefer project-level data but keep tool-local fallback for portability.
RAW_DIR = PROJECT_ROOT / "Raw data"
if not RAW_DIR.exists():
    RAW_DIR = SCRIPT_DIR / "Raw data"
SOURCE_CSV = RAW_DIR / "pmc6374777_table1.csv"
OUTPUT_CSV = RAW_DIR / "unified_pfas_ki_smiles.csv"

EXCLUDED_SECTIONS = {"Hydrocarbons", "OIL-Water"}

FORMULA_BY_ACRONYM = {
    "PFAA": "CF3CO2Na",
    "PFPrA": "C2F5CO2Na",
    "PFBA": "C3F7CO2Na",
    "PFPeA": "C4F9CO2Na",
    "PFHxA": "C5F11CO2Na",
    "PFHpA": "C6F13CO2Na",
    "PFOA": "C7F15CO2Na",
    "PFNA": "C8F17CO2Na",
    "PFDA": "C9F19CO2Na",
    "PFUnA": "C10F21CO2Na",
    "Iso-PFOA": "(CF3)2CF(CF2)4CO2Na",
    "Iso-PFDA": "(CF3)2CF(CF2)6CO2Na",
    "PFBS": "C4F9SO3K",
    "PFHxS": "C6F13SO3K",
    "PFHpS": "C7F15SO3Na",
    "PFOS": "C8F17SO3K",
    "PFNS": "C9F19SO3K",
    "9H-PFNA": "HCF2(CF2)7CO2Na",
    "7H-PFHpA": "HCF2(CF2)5CO2NH4",
    "SPBS": "C9F17OC6H4SO3Na",
    "FC-53": "CF3(CF2)5O(CF2)2SO3K",
    "TDFHD": "CF3(CF2)3CF(CF3)(CH2)10CO2Na",
    "F9-CTAB": "CF3(CF2)3(CH2)11N(CH3)3Br",
    "F12-CTAB": "HCF2(CF2)5(CH2)9N(CH3)3Br",
    "F17-CTAB": "CF3(CF2)7(CH2)6N(CH3)3Br",
    "UDFOS": "CH3CH2CHF(CF2)5SO3Na",
    "NFHES": "CF3CF2O(CF2)2(CH2)2SO3Na",
    "UDFHES": "CF3(CF2)2O(CF2)2(CH2)2SO3Na",
    "TDFP": "CF3(CF2)2C(CF3)2CH2CO2Na",
    "TDHP": "CF3(CF2)2C(CF3)2(CH2)2CO2Na",
    "HDFPEC": "CF3(CF2)2OCF(CF3)CF2OCF(CF3)CO2Na",
    "TDFPBP": "CF3(CF2)2C(CF3)2CH2C6H4PO3Li2",
    "TDFTDE": "CF3(CF2)5C2H4SC2H4(CH2CH2O)2OH",
    "TDFTTE": "CF3(CF2)5C2H4SC2H4(CH2CH2O)3OH",
    "TDFTPE": "CF3(CF2)5C2H4SC2H4(CH2CH2O)5OH",
    "PFOA-amide": "C7F15CONHC3H6OH",
    "NFTME": "CF3(CF2)3CH2O(CH2CH2O)3CH3",
    "TDFTME": "CF3(CF2)5CH2O(CH2CH2O)3CH3",
    "HOFTME": "HCF2(CF2)3CH2O(CH2CH2O)3CH3",
    "HDDFTME": "HCF2(CF2)5CH2O(CH2CH2O)3CH3",
    "8:1 FTOH": "CF3(CF2)7CH2OH",
    "FC8diol": "(CF2)6(CH2)2(OH)2",
}


def canonicalize(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


def perfluoro_chain(total_carbons: int, terminal: str = "CF3") -> str:
    units = ["C(F)(F)"] * (total_carbons - 1)
    units.append("C(F)(F)F" if terminal == "CF3" else "C(F)F")
    return "".join(units)


def carboxylate(counterion: str, tail_carbons: int, terminal: str = "CF3") -> str:
    return canonicalize(f"{counterion}.[O-]C(=O){perfluoro_chain(tail_carbons, terminal)}")


def sulfonate(counterion: str, tail_carbons: int) -> str:
    return canonicalize(f"{counterion}.[O-]S(=O)(=O){perfluoro_chain(tail_carbons)}")


SMILES_BY_ACRONYM = {
    "PFAA": carboxylate("[Na+]", 1),
    "PFPrA": carboxylate("[Na+]", 2),
    "PFBA": carboxylate("[Na+]", 3),
    "PFPeA": carboxylate("[Na+]", 4),
    "PFHxA": carboxylate("[Na+]", 5),
    "PFHpA": carboxylate("[Na+]", 6),
    "PFOA": carboxylate("[Na+]", 7),
    "PFNA": carboxylate("[Na+]", 8),
    "PFDA": carboxylate("[Na+]", 9),
    "PFUnA": carboxylate("[Na+]", 10),
    "Iso-PFOA": canonicalize(
        "[Na+].[O-]C(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "Iso-PFDA": canonicalize(
        "[Na+].[O-]C(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "PFBS": sulfonate("[K+]", 4),
    "PFHxS": sulfonate("[K+]", 6),
    "PFHpS": sulfonate("[Na+]", 7),
    "PFOS": sulfonate("[K+]", 8),
    "PFNS": sulfonate("[K+]", 9),
    "9H-PFNA": carboxylate("[Na+]", 8, terminal="CF2H"),
    "7H-PFHpA": carboxylate("[NH4+]", 6, terminal="CF2H"),
    "SPBS": canonicalize(
        "[Na+].[O-]S(=O)(=O)c1ccc(OC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)cc1"
    ),
    "FC-53": canonicalize(
        "[K+].[O-]S(=O)(=O)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "TDFHD": canonicalize(
        "[Na+].[O-]C(=O)CCCCCCCCCC(C(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "F9-CTAB": canonicalize(
        "[Br-].C[N+](C)(C)CCCCCCCCCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "F12-CTAB": canonicalize(
        "[Br-].C[N+](C)(C)CCCCCCCCCC(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "F17-CTAB": canonicalize(
        "[Br-].C[N+](C)(C)CCCCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "UDFOS": canonicalize(
        "[Na+].[O-]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)CC"
    ),
    "NFHES": canonicalize(
        "[Na+].[O-]S(=O)(=O)CCOC(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "UDFHES": canonicalize(
        "[Na+].[O-]S(=O)(=O)CCOC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "TDFP": canonicalize(
        "[Na+].[O-]C(=O)CC(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "TDHP": canonicalize(
        "[Na+].[O-]C(=O)CCC(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "HDFPEC": canonicalize(
        "[Na+].[O-]C(=O)C(F)(C(F)(F)F)OC(F)(F)C(F)(C(F)(F)F)OC(F)(F)C(F)(F)C(F)(F)F"
    ),
    "TDFPBP": canonicalize(
        "[Li+].[Li+].[O-]P(=O)([O-])c1ccc(CC(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F)cc1"
    ),
    "TDFTDE": canonicalize(
        "OCCOCCOCCSCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "TDFTTE": canonicalize(
        "OCCOCCOCCOCCSCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "TDFTPE": canonicalize(
        "OCCOCCOCCOCCOCCOCCSCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "PFOA-amide": canonicalize(
        "CC(O)CNC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "NFTME": canonicalize(
        "COCCOCCOCCOCCOCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "TDFTME": canonicalize(
        "COCCOCCOCCOCCOCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "HOFTME": canonicalize(
        "COCCOCCOCCOCCOCC(F)(F)C(F)(F)C(F)(F)C(F)F"
    ),
    "HDDFTME": canonicalize(
        "COCCOCCOCCOCCOCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)F"
    ),
    "8:1 FTOH": canonicalize(
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ),
    "FC8diol": canonicalize(
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CO"
    ),
}


def build_dataset() -> pd.DataFrame:
    df = pd.read_csv(SOURCE_CSV)
    df = df.loc[~df["section"].isin(EXCLUDED_SECTIONS)].copy()
    df = df.drop_duplicates(subset=["acronym"], keep="first")

    missing_formulas = sorted(set(df["acronym"]) - set(FORMULA_BY_ACRONYM))
    missing_smiles = sorted(set(df["acronym"]) - set(SMILES_BY_ACRONYM))
    if missing_formulas or missing_smiles:
        raise KeyError(
            f"Missing mappings. formulas={missing_formulas}, smiles={missing_smiles}"
        )

    df["source_formula"] = df["acronym"].map(FORMULA_BY_ACRONYM)
    df["smiles"] = df["acronym"].map(SMILES_BY_ACRONYM)
    df["log10_Ki_cm"] = np.log10(df["Ki_cm"])
    df["source_file"] = SOURCE_CSV.name
    df["source_reference"] = "Brusseau et al. 2019 Table 1"
    df["smiles_representation"] = "source-reported ionic or neutral species"

    return df[
        [
            "acronym",
            "name",
            "section",
            "Ki_cm",
            "log10_Ki_cm",
            "source_formula",
            "smiles",
            "source_reference",
            "source_file",
            "smiles_representation",
        ]
    ].sort_values(["section", "acronym"], ignore_index=True)


def main() -> None:
    df = build_dataset()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
