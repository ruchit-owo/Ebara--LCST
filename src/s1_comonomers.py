from pathlib import Path
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "lcst.csv"
OUT = ROOT / "data" / "processed" / "comonomers.csv"
FIG = ROOT / "results" / "figures" / "comonomer_structures.png"

COLUMNS = [
    "entry_id", "role", "source", "source_location", "read_method", "base_monomer",
    "comonomer_name", "comonomer_smiles", "comonomer_mol_frac", "frac_type",
    "Mw", "Mw_method", "polymer_conc", "conc_units", "solvent", "pH", "salt",
    "lcst_method", "heating_rate", "lcst_C", "notes",
]


def main():
    df = pd.read_csv(RAW)
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        sys.exit(f"missing columns: {missing}")
    if df["entry_id"].duplicated().any():
        sys.exit(f"duplicate entry_id: {df.loc[df.entry_id.duplicated(), 'entry_id'].tolist()}")

    co = df[df["comonomer_name"].notna() & (df["comonomer_name"] != "none")]
    names = co.groupby("comonomer_name")["comonomer_smiles"].nunique()
    if (names > 1).any():
        sys.exit(f"name with more than one SMILES: {names[names > 1].index.tolist()}")

    rows, bad = [], []
    for name, g in co.groupby("comonomer_name"):
        smi = g["comonomer_smiles"].iloc[0]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            bad.append(name)
            continue
        rows.append({
            "comonomer_name": name,
            "comonomer_smiles": smi,
            "canonical_smiles": Chem.MolToSmiles(mol),
            "n_rows": len(g),
            "n_sources": g["source"].nunique(),
            "frac_min": g["comonomer_mol_frac"].min(),
            "frac_max": g["comonomer_mol_frac"].max(),
        })
    if bad:
        sys.exit(f"SMILES failed to parse: {bad}")

    out = pd.DataFrame(rows).sort_values("n_rows", ascending=False)
    dup = out["canonical_smiles"].duplicated(keep=False)
    if dup.any():
        sys.exit(f"different names, same structure: {out.loc[dup, 'comonomer_name'].tolist()}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)

    mols = [Chem.MolFromSmiles(s) for s in out["comonomer_smiles"]]
    img = Draw.MolsToGridImage(mols, molsPerRow=6, subImgSize=(260, 220),
                               legends=[f"{n} ({k})" for n, k in zip(out.comonomer_name, out.n_rows)])
    FIG.parent.mkdir(parents=True, exist_ok=True)
    img.save(FIG)

    print(f"rows: {len(df)}")
    print(f"homopolymer rows: {len(df) - len(co)}")
    print(f"distinct comonomers: {len(out)}")
    print(f"sources: {df['source'].nunique()}")
    print("\nrows per lcst_method:")
    print(df["lcst_method"].value_counts(dropna=False).to_string())
    print("\nrows per frac_type:")
    print(co["frac_type"].value_counts(dropna=False).to_string())
    print("\nrows per comonomer:")
    print(out[["comonomer_name", "n_rows", "n_sources", "frac_min", "frac_max"]].to_string(index=False))


if __name__ == "__main__":
    main()
