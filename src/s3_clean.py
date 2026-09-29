from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
RDKIT = PROC / "descriptors_rdkit.csv"
ALVA = PROC / "descriptors_alvadesc.csv"
MODEL_TABLE = PROC / "model_table.csv"
OUT = PROC / "features_clean.csv"
SUMMARY = RES / "feature_cleaning.csv"
DROPLOG = RES / "feature_drop_log.csv"

MIN_OFF_MODE = 4
CORR_MAX = 0.95


def main():
    rd = pd.read_csv(RDKIT).set_index("comonomer_name").add_prefix("rd_")
    X = rd
    if ALVA.exists():
        al = pd.read_csv(ALVA).set_index("comonomer_name").add_prefix("al_")
        X = rd.join(al, how="left").copy()

    mt = pd.read_csv(MODEL_TABLE)
    used = sorted(set(mt.loc[mt.comonomer_name != "none", "comonomer_name"]))
    C = X.loc[used]
    n = len(C)

    log, steps = [], []
    cols = list(C.columns)

    def step(name, keep, reason):
        nonlocal cols
        dropped = [c for c in cols if c not in keep]
        log.extend({"feature": c, "step": name, "reason": reason} for c in dropped)
        steps.append({"step": name, "n_before": len(cols), "n_dropped": len(dropped),
                      "n_after": len(keep), "rule": reason})
        cols = keep

    step("missing", [c for c in cols if C[c].notna().all()],
         "value missing for at least one comonomer")
    step("constant", [c for c in cols if C[c].nunique() > 1],
         "same value for every comonomer")
    near = []
    for c in cols:
        mode_n = C[c].round(8).value_counts().iloc[0]
        if n - mode_n >= MIN_OFF_MODE:
            near.append(c)
    step("near_constant", near,
         f"fewer than {MIN_OFF_MODE} of {n} comonomers differ from the most common value")

    Z = C[cols].astype(float)
    corr = Z.corr().abs().values
    keep_idx = []
    for i in range(len(cols)):
        if all(corr[i, j] <= CORR_MAX for j in keep_idx):
            keep_idx.append(i)
    step("correlated", [cols[i] for i in keep_idx],
         f"|r| > {CORR_MAX} with a feature kept earlier (RDKit columns first, then alvaDesc export order)")

    out = X.loc[["NIPAM"] + used, cols].reset_index()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    RES.mkdir(exist_ok=True)
    pd.DataFrame(steps).to_csv(SUMMARY, index=False)
    pd.DataFrame(log).to_csv(DROPLOG, index=False)

    print(f"comonomers used for filtering: {n}")
    print(pd.DataFrame(steps)[["step", "n_before", "n_dropped", "n_after"]].to_string(index=False))
    kept_rd = [c for c in cols if c.startswith("rd_")]
    print(f"\nsurviving: {len(cols)} ({len(kept_rd)} RDKit, {len(cols) - len(kept_rd)} alvaDesc)")
    print("RDKit kept:", ", ".join(c[3:] for c in kept_rd))
    print(f"features / distinct comonomers = {len(cols) / n:.1f}; further selection happens inside CV folds in Stage 4")


if __name__ == "__main__":
    main()
