from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Lipinski, rdMolDescriptors
from sklearn.linear_model import LinearRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s4_model import ridge, CHARGED, SEED
from s5_trends import phys_rows, pnipam_reference, PHYS, INK, INK2, BLUE, ORANGE, AQUA

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
FIG = RES / "figures"

N_BOOT = 500
TARGET = (38.0, 40.0)
NEW = {
    "HPMA": ("C=C(C)C(=O)NCC(C)O", "N-(2-hydroxypropyl)methacrylamide"),
    "NAM": ("C=CC(=O)N1CCOCC1", "N-acryloylmorpholine"),
    "NVP": ("C=CN1CCCC1=O", "N-vinylpyrrolidone"),
    "MEA": ("C=CC(=O)OCCOC", "2-methoxyethyl acrylate"),
}
HEADLINE = "NAM"


def desc_of(smi):
    m = Chem.MolFromSmiles(smi)
    return {"logP": Crippen.MolLogP(m), "TPSA": rdMolDescriptors.CalcTPSA(m),
            "HBD": Lipinski.NumHDonors(m), "charged_at_pH7": 0}


def load_table():
    mt = pd.read_csv(PROC / "model_table.csv")
    mt = mt[~mt.comonomer_name.isin(CHARGED)].reset_index(drop=True)
    mt["frac"] = mt.comonomer_mol_frac.astype(float)
    return mt


def boot_predict(train, X, resid_pool, rng):
    groups = np.array(sorted(train.comonomer_name.unique()))
    out = np.empty((N_BOOT, len(X)))
    for b in range(N_BOOT):
        pick = rng.choice(groups, len(groups), replace=True)
        B = pd.concat([train[train.comonomer_name == g] for g in pick])
        out[b] = ridge().fit(B[PHYS], B.lcst_C).predict(X) + rng.choice(resid_pool)
    return out


def logo_residuals(train):
    res = []
    for g in train.comonomer_name.unique():
        tr, te = train[train.comonomer_name != g], train[train.comonomer_name == g]
        res.extend(te.lcst_C - ridge().fit(tr[PHYS], tr.lcst_C).predict(te[PHYS]))
    return np.array(res)


def applicability(train_desc, name, d, fps, fp_new):
    lo, hi = train_desc[["logP", "TPSA", "HBD"]].min(), train_desc[["logP", "TPSA", "HBD"]].max()
    inside = {k: bool(lo[k] <= d[k] <= hi[k]) for k in ["logP", "TPSA", "HBD"]}
    Z = train_desc[["logP", "TPSA", "HBD"]].astype(float)
    mu, sd = Z.mean(), Z.std()
    Xs = ((Z - mu) / sd).values
    x = ((pd.Series({k: d[k] for k in ["logP", "TPSA", "HBD"]}) - mu) / sd).values
    Xa = np.column_stack([np.ones(len(Xs)), Xs])
    xa = np.concatenate([[1.0], x])
    h = float(xa @ np.linalg.pinv(Xa.T @ Xa) @ xa)
    h_star = 3 * Xa.shape[1] / len(Xa)
    sims = {n: DataStructs.TanimotoSimilarity(fp_new, f) for n, f in fps.items()}
    nn = max(sims, key=sims.get)
    return {"comonomer": name, "logP": d["logP"], "TPSA": d["TPSA"], "HBD": d["HBD"],
            "logP_in_range": inside["logP"], "TPSA_in_range": inside["TPSA"], "HBD_in_range": inside["HBD"],
            "leverage": h, "leverage_threshold": h_star,
            "nearest_train_comonomer": nn, "tanimoto_nn": sims[nn],
            "domain": "inside" if all(inside.values()) and h <= h_star else "extrapolation"}


def main():
    mt = load_table()
    train_all = mt[mt.role == "train"].reset_index(drop=True)
    check = mt[mt.role == "ebara_check"].reset_index(drop=True)
    assert not set(check.source) & set(train_all.source), "Ebara source found in training data"
    rng = np.random.default_rng(SEED)
    desc = pd.read_csv(PROC / "descriptors_rdkit.csv").set_index("comonomer_name")
    nipam = desc.loc["NIPAM"]

    # A. Ebara check: HMAAm fully held out (headline) and with non-Ebara HMAAm rows kept (reference)
    rows = []
    for label, train in [("HMAAm_held_out", train_all[train_all.comonomer_name != "HMAAm"]),
                         ("HMAAm_in_training", train_all)]:
        m = ridge().fit(train[PHYS], train.lcst_C)
        p = m.predict(check[PHYS])
        boot = boot_predict(train, check[PHYS], logo_residuals(train), rng)
        lf = LinearRegression().fit(train[["frac"]], train.lcst_C).predict(check[["frac"]])
        for i, r in check.iterrows():
            rows.append({"training": label, "entry_id": r.entry_id, "source": r.source, "HMAAm_frac": r.frac,
                         "frac_type": r.frac_type, "solvent": r.solvent, "measured": r.lcst_C,
                         "pred_ridge_physics": p[i], "pi05": np.quantile(boot[:, i], 0.05),
                         "pi95": np.quantile(boot[:, i], 0.95), "pred_linear_frac": lf[i], "pred_const32": 32.0})
    eb = pd.DataFrame(rows)
    eb["err_ridge"] = eb.pred_ridge_physics - eb.measured
    eb["in_90pi"] = (eb.measured >= eb.pi05) & (eb.measured <= eb.pi95)
    eb.round(2).to_csv(RES / "ebara_check.csv", index=False)
    summ = eb.groupby("training").apply(lambda g: pd.Series({
        "MAE_ridge": g.err_ridge.abs().mean(), "MAE_linear_frac": (g.pred_linear_frac - g.measured).abs().mean(),
        "MAE_const32": (32 - g.measured).abs().mean(), "n_in_90pi": int(g.in_90pi.sum()), "n": len(g)}),
        include_groups=False)

    hm = mt[mt.comonomer_name == "HMAAm"].copy()
    ref = pnipam_reference(hm, mt[mt.comonomer_name == "none"])
    hm["dT"] = hm.lcst_C - ref
    src_slope = hm.groupby(["role", "source"]).apply(
        lambda g: pd.Series({"n": len(g), "slope_per_10mol": 0.1 * (g.dT * g.frac).sum() / (g.frac ** 2).sum()}),
        include_groups=False).reset_index()
    src_slope.round(2).to_csv(RES / "hmaam_source_slopes.csv", index=False)

    # B. New copolymers: prediction, uncertainty, applicability domain
    train = train_all
    m = ridge().fit(train[PHYS], train.lcst_C)
    resid = logo_residuals(train)
    fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    train_names = sorted(set(train.comonomer_name) - {"none"})
    fps = {n: fpgen.GetFingerprint(Chem.MolFromSmiles(s)) for n, s in
           pd.read_csv(PROC / "comonomers.csv").set_index("comonomer_name").comonomer_smiles.items() if n in train_names}
    tdesc = desc.loc[train_names]
    frac_max_train = train.frac.max()
    grid = np.round(np.arange(0, 0.51, 0.01), 2)
    new_rows, ad_rows, curves, contrib = [], [], {}, []
    for name, (smi, full_name) in NEW.items():
        d = desc_of(smi)
        dd = pd.DataFrame([d], index=[name])
        X = phys_rows(pd.concat([desc, dd]), nipam, name, grid, conc=0.1, buffer=0)
        pred = m.predict(X)
        boot = boot_predict(train, X, resid, rng)
        lo, hi = np.quantile(boot, 0.05, axis=0), np.quantile(boot, 0.95, axis=0)
        curves[name] = (grid, pred, lo, hi)
        hit = grid[(pred >= TARGET[0]) & (pred <= TARGET[1])]
        ad = applicability(tdesc, name, d, fps, fpgen.GetFingerprint(Chem.MolFromSmiles(smi)))
        ad_rows.append({**ad, "full_name": full_name, "smiles": smi})
        for f in [0.05, 0.10, 0.20, 0.30]:
            i = int(np.where(np.isclose(grid, f))[0][0])
            new_rows.append({"comonomer": name, "frac": f, "pred": pred[i], "pi05": lo[i], "pi95": hi[i]})
        sc = m[1].scale_
        coefs = pd.Series(m[-1].coef_ / sc, index=PHYS)
        x20 = X.loc[np.isclose(grid, 0.20)].iloc[0]
        contrib.append({"comonomer": name, "frac": 0.20,
                        "dlogP_term_C": coefs["frac_x_dlogP"] * x20["frac_x_dlogP"],
                        "dTPSA_term_C": coefs["frac_x_dTPSA"] * x20["frac_x_dTPSA"],
                        "dHBD_term_C": coefs["frac_x_dHBD"] * x20["frac_x_dHBD"],
                        "frac_term_C": coefs["frac"] * x20["frac"]})
        new_rows.append({"comonomer": name, "frac": f"{hit.min():.2f}-{hit.max():.2f}" if len(hit) else "none <= 0.50",
                         "pred": np.nan, "pi05": np.nan, "pi95": np.nan, "note": f"fraction range predicted to give {TARGET[0]:.0f}-{TARGET[1]:.0f} C"})
    newp = pd.DataFrame(new_rows)
    newp.round(2).to_csv(RES / "new_predictions.csv", index=False)
    ct = pd.DataFrame(contrib)
    for c in ["HMAAm", "DMAAm", "HIPAAm"]:
        x20 = phys_rows(desc, nipam, c, [0.20], conc=0.1, buffer=0).iloc[0]
        coefs = pd.Series(m[-1].coef_ / m[1].scale_, index=PHYS)
        ct.loc[len(ct)] = {"comonomer": c + " (training)", "frac": 0.20,
                           "dlogP_term_C": coefs["frac_x_dlogP"] * x20["frac_x_dlogP"],
                           "dTPSA_term_C": coefs["frac_x_dTPSA"] * x20["frac_x_dTPSA"],
                           "dHBD_term_C": coefs["frac_x_dHBD"] * x20["frac_x_dHBD"],
                           "frac_term_C": coefs["frac"] * x20["frac"]}
    ct.round(2).to_csv(RES / "prediction_term_breakdown.csv", index=False)
    adt = pd.DataFrame(ad_rows)
    adt.round(3).to_csv(RES / "applicability_domain.csv", index=False)

    # Figures
    FIG.mkdir(parents=True, exist_ok=True)
    fig, a = plt.subplots(figsize=(5.8, 4.2))
    tr0 = train_all[train_all.comonomer_name != "HMAAm"]
    m0 = ridge().fit(tr0[PHYS], tr0.lcst_C)
    g = np.linspace(0, 0.32, 60)
    band = boot_predict(tr0, phys_rows(desc, nipam, "HMAAm", g, conc=0.1, buffer=0), logo_residuals(tr0), rng)
    a.fill_between(g * 100, np.quantile(band, 0.05, axis=0), np.quantile(band, 0.95, axis=0), color=ORANGE, alpha=0.15, lw=0)
    a.plot(g * 100, m0.predict(phys_rows(desc, nipam, "HMAAm", g, conc=0.1, buffer=0)), color=ORANGE, lw=2,
           label="predicted, HMAAm never seen (water, 0.1 wt%)")
    oth = hm[hm.role == "train"]
    a.scatter(oth.frac * 100, oth.lcst_C, s=30, color=AQUA, edgecolor="white", linewidth=1.2, zorder=3, label="other groups (training data)")
    a.scatter(check.frac * 100, check.lcst_C, s=44, marker="s", color=BLUE, edgecolor="white", linewidth=1.2, zorder=4, label="Ebara group (check set)")
    a.set_xlabel("HMAAm (mol%)"); a.set_ylabel("LCST (C)")
    a.set_title("NIPAAm-co-HMAAm: Ebara check", fontsize=9, color=INK, loc="left")
    a.legend(loc="upper left", fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "ebara_check.png", dpi=200); plt.close(fig)

    fig, a = plt.subplots(figsize=(5.8, 4.2))
    gx, gp, glo, ghi = curves[HEADLINE]
    a.axhspan(*TARGET, color=INK2, alpha=0.08, lw=0)
    a.text(1, TARGET[1] + 0.4, "target 38-40 C", fontsize=7, color=INK2)
    a.fill_between(gx * 100, glo, ghi, color=BLUE, alpha=0.15, lw=0, label="90% interval")
    a.plot(gx * 100, gp, color=BLUE, lw=2, label="predicted")
    a.set_xlabel(f"{HEADLINE} (mol%)"); a.set_ylabel("Predicted LCST (C), water 0.1 wt%")
    dom = adt.set_index("comonomer").loc[HEADLINE]
    a.set_title(f"NIPAM-co-{HEADLINE} (not in data): domain = {dom.domain}", fontsize=9, color=INK, loc="left")
    a.legend(loc="upper left", fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "new_prediction.png", dpi=200); plt.close(fig)

    pd.set_option("display.width", 200)
    print("A. Ebara check")
    print(eb[["training", "entry_id", "source", "HMAAm_frac", "solvent", "measured", "pred_ridge_physics", "pi05", "pi95",
              "pred_linear_frac"]].round(1).to_string(index=False))
    print("\n", summ.round(2).to_string())
    print("\nHMAAm shift per 10 mol% by source:")
    print(src_slope.round(2).to_string(index=False))
    print("\nB. New copolymers (water, 0.1 wt%)")
    print(newp.round(1).to_string(index=False))
    print("\nTerm breakdown at 20 mol% (C relative to PNIPAM):")
    print(ct.round(1).to_string(index=False))
    print("\nApplicability domain (training comonomers, logP/TPSA/HBD):")
    print(adt[["comonomer", "logP", "TPSA", "HBD", "logP_in_range", "TPSA_in_range", "HBD_in_range", "leverage",
               "leverage_threshold", "nearest_train_comonomer", "tanimoto_nn", "domain"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
