from pathlib import Path
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "lcst.csv"
COMONOMERS = ROOT / "data" / "processed" / "comonomers.csv"
PROC = ROOT / "data" / "processed"
ALVA_IN = PROC / "alvadesc_input.smi"
ALVA_EXPORT = ROOT / "data" / "raw" / "alvadesc_export.txt"
RDKIT_OUT = PROC / "descriptors_rdkit.csv"
OUT = PROC / "model_table.csv"

NIPAM_SMILES = "C=CC(=O)NC(C)C"

CONC_TO_WT = {"wt%": 1.0, "w/v%": 1.0, "mg/mL": 0.1}


def rdkit_desc(smi):
    m = Chem.MolFromSmiles(smi)
    return {
        "logP": Crippen.MolLogP(m),
        "MR": Crippen.MolMR(m),
        "TPSA": rdMolDescriptors.CalcTPSA(m),
        "HBD": Lipinski.NumHDonors(m),
        "HBA": Lipinski.NumHAcceptors(m),
        "MW": Descriptors.MolWt(m),
        "RotB": Lipinski.NumRotatableBonds(m),
        "HeavyAtoms": m.GetNumHeavyAtoms(),
        "AromRings": rdMolDescriptors.CalcNumAromaticRings(m),
        "FracCSP3": rdMolDescriptors.CalcFractionCSP3(m),
        "charged_at_pH7": int(m.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))
                              or m.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H0;!$(NC=O);!$(N-a);!$(n)](C)(C)C"))),
    }


def main():
    df = pd.read_csv(RAW)
    co = pd.read_csv(COMONOMERS)

    names = ["NIPAM"] + co["comonomer_name"].tolist()
    smiles = [NIPAM_SMILES] + co["comonomer_smiles"].tolist()
    ALVA_IN.write_text("".join(f"{s}\t{n}\n" for n, s in zip(names, smiles)))

    desc = pd.DataFrame([{"comonomer_name": n, **rdkit_desc(s)} for n, s in zip(names, smiles)])
    desc.to_csv(RDKIT_OUT, index=False)
    nipam = desc.iloc[0]
    cdesc = desc.iloc[1:].copy()
    dcols = [c for c in cdesc.columns if c != "comonomer_name"]
    cdesc = cdesc.rename(columns={c: f"co_{c}" for c in dcols})

    if ALVA_EXPORT.exists():
        alva = pd.read_csv(ALVA_EXPORT, sep="\t")
        alva = alva.rename(columns={alva.columns[0]: "comonomer_name"})
        alva = alva[alva["comonomer_name"] != "NIPAM"]
        alva = alva.rename(columns={c: f"alva_{c}" for c in alva.columns if c != "comonomer_name"})
        cdesc = cdesc.merge(alva, on="comonomer_name", how="left")
        print(f"alvaDesc descriptors joined: {alva.shape[1] - 1}")
    else:
        print("no alvaDesc export found; RDKit set only")

    n0 = len(df)
    t = df.merge(cdesc, on="comonomer_name", how="left")
    if len(t) != n0:
        sys.exit(f"row count changed on merge: {n0} -> {len(t)}")

    f = t["comonomer_mol_frac"].astype(float)
    t["frac"] = f
    t["frac_x_logP"] = (f * t["co_logP"]).fillna(0.0)
    t["frac_x_dlogP"] = (f * (t["co_logP"] - nipam["logP"])).fillna(0.0)
    t["frac_x_dTPSA"] = (f * (t["co_TPSA"] - nipam["TPSA"])).fillna(0.0)
    t["frac_x_dHBD"] = (f * (t["co_HBD"] - nipam["HBD"])).fillna(0.0)
    t["frac_x_charged"] = (f * t["co_charged_at_pH7"]).fillna(0.0)

    t["conc_wt"] = t["polymer_conc"] * t["conc_units"].map(CONC_TO_WT)
    t["log10_M"] = np.log10(pd.to_numeric(t["Mw"], errors="coerce"))
    t["is_buffer"] = t["solvent"].isin(["PBS", "phosphate_buffer"]).astype(int)
    t["is_secondary"] = t["notes"].fillna("").str.startswith("secondary").astype(int)

    co_rows = t[t["comonomer_name"] != "none"]
    miss = co_rows["co_logP"].isna().sum()
    if miss:
        sys.exit(f"{miss} comonomer rows without descriptors")

    t.to_csv(OUT, index=False)
    print(f"rows in: {n0}, rows out: {len(t)}")
    print(f"train rows: {(t.role == 'train').sum()}, ebara_check rows: {(t.role == 'ebara_check').sum()}")
    print(f"missing conc_wt: {t.conc_wt.isna().sum()}, missing log10_M: {t.log10_M.isna().sum()}")
    print("\nRDKit descriptors (NIPAM first):")
    print(desc.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
