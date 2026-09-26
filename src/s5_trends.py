from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s4_model import load, ridge, CHARGED, COND, SEED

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures"
RDKIT = ROOT / "data" / "processed" / "descriptors_rdkit.csv"

PHYS = ["frac", "frac_x_dlogP", "frac_x_dTPSA", "frac_x_dHBD", "frac_x_charged"] + COND
FIX_FRAC = 0.10
N_BOOT = 200
LABEL_OFFSET = {"HEAAm": (4, -9), "NtBAm": (-4, -11), "NEAm": (3, -9), "NnPAm": (3, 4), "AProOMe": (-10, -11), "HMAAm": (4, 4)}
SERIES = ["DMAAm", "HMAAm", "NtBAm", "AAm", "HIPAAm", "NEAm"]

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False})


def phys_rows(desc, nipam, name, fracs, conc, buffer=0):
    d = desc.loc[name]
    f = np.asarray(fracs, float)
    return pd.DataFrame({
        "frac": f,
        "frac_x_dlogP": f * (d.logP - nipam.logP),
        "frac_x_dTPSA": f * (d.TPSA - nipam.TPSA),
        "frac_x_dHBD": f * (d.HBD - nipam.HBD),
        "frac_x_charged": f * d.charged_at_pH7,
        "conc_wt": conc,
        "is_buffer": buffer,
    })


def measured_slope(g, pnipam_ref):
    g = g.copy()
    ref = g.source.map(pnipam_ref).fillna(32.0)
    low = g[g.frac <= 0.35]
    use = low if len(low) else g
    d = use.lcst_C - ref.loc[use.index]
    return float((d * use.frac).sum() / (use.frac ** 2).sum())


def main():
    mt, _, _ = load()
    D = mt[~mt.comonomer_name.isin(CHARGED)].reset_index(drop=True)
    desc = pd.read_csv(RDKIT).set_index("comonomer_name")
    nipam = desc.loc["NIPAM"]
    conc_med = float(D.conc_wt.median())
    comonomers = sorted(set(D.comonomer_name) - {"none"})
    pnipam_ref = D[D.comonomer_name == "none"].groupby("source").lcst_C.mean()

    full = ridge().fit(D[PHYS], D.lcst_C)

    # Trend 1: predicted LCST at fixed fraction vs comonomer logP (full fit and comonomer held out)
    t1 = []
    for c in comonomers:
        tr = D[D.comonomer_name != c]
        m = ridge().fit(tr[PHYS], tr.lcst_C)
        X = phys_rows(desc, nipam, c, [FIX_FRAC], conc_med)
        g = D[D.comonomer_name == c]
        t1.append({"comonomer_name": c, "logP": desc.loc[c, "logP"], "TPSA": desc.loc[c, "TPSA"],
                   "HBD": desc.loc[c, "HBD"], "n_rows": len(g),
                   "pred_full_at_0.10": float(full.predict(X)[0]),
                   "pred_heldout_at_0.10": float(m.predict(X)[0]),
                   "measured_slope_per_10mol": 0.1 * measured_slope(g, pnipam_ref)})
    t1 = pd.DataFrame(t1)
    t1.round(3).to_csv(RES / "trend1_fixed_fraction.csv", index=False)
    rho = {
        "pred_full_vs_logP": spearmanr(t1["pred_full_at_0.10"], t1.logP),
        "pred_heldout_vs_logP": spearmanr(t1["pred_heldout_at_0.10"], t1.logP),
        "measured_slope_vs_logP": spearmanr(t1.measured_slope_per_10mol, t1.logP),
        "pred_heldout_vs_measured_slope": spearmanr(t1["pred_heldout_at_0.10"], t1.measured_slope_per_10mol),
    }

    # Trend 2: composition series, comonomer held out
    t2, curves = [], {}
    for c in SERIES:
        tr = D[D.comonomer_name != c]
        m = ridge().fit(tr[PHYS], tr.lcst_C)
        g = D[D.comonomer_name == c]
        p = m.predict(g[PHYS])
        grid = np.linspace(0, g.frac.max() * 1.05, 50)
        curves[c] = (grid, m.predict(phys_rows(desc, nipam, c, grid, conc_med)))
        ref = g.source.map(pnipam_ref).fillna(32.0)
        s_meas = float(((g.lcst_C - ref) * g.frac).sum() / (g.frac ** 2).sum())
        s_pred = float(((p - ref) * g.frac).sum() / (g.frac ** 2).sum())
        t2.append({"comonomer_name": c, "n_rows": len(g), "n_sources": g.source.nunique(),
                   "frac_range": f"{g.frac.min():.2f}-{g.frac.max():.2f}",
                   "spearman_meas_vs_pred": spearmanr(g.lcst_C, p).statistic if g.frac.nunique() > 2 else np.nan,
                   "MAE": float(np.abs(p - g.lcst_C).mean()),
                   "slope_measured_per_10mol": 0.1 * s_meas, "slope_predicted_per_10mol": 0.1 * s_pred})
    t2 = pd.DataFrame(t2)
    t2.round(3).to_csv(RES / "trend2_series.csv", index=False)

    # Feature reliance: standardized coefficients + sign stability over comonomer bootstrap
    coef = pd.Series(full[-1].coef_, index=PHYS)
    rng = np.random.default_rng(SEED)
    groups = np.array(sorted(D.comonomer_name.unique()))
    boot = []
    for _ in range(N_BOOT):
        pick = rng.choice(groups, len(groups), replace=True)
        B = pd.concat([D[D.comonomer_name == gname] for gname in pick])
        boot.append(ridge().fit(B[PHYS], B.lcst_C)[-1].coef_)
    boot = np.array(boot)
    ct = pd.DataFrame({"feature": PHYS, "std_coef": coef.values,
                       "boot_q05": np.quantile(boot, 0.05, axis=0), "boot_q95": np.quantile(boot, 0.95, axis=0),
                       "frac_boot_same_sign": (np.sign(boot) == np.sign(coef.values)).mean(axis=0)})
    ct = ct[ct.feature != "frac_x_charged"]
    ct.round(3).to_csv(RES / "coefficients_ridge_physics.csv", index=False)

    # Sensitivity: drop feature groups, grouped 5-fold MAE (5 shuffled repeats)
    def gmae(cols, sub):
        e = []
        for rep in range(5):
            pr = np.zeros(len(sub))
            for tr, te in GroupKFold(5, shuffle=True, random_state=SEED + rep).split(sub, groups=sub.group):
                pr[te] = ridge().fit(sub.iloc[tr][cols], sub.iloc[tr].lcst_C).predict(sub.iloc[te][cols])
            e.append(np.abs(pr - sub.lcst_C).mean())
        return float(np.mean(e))
    no_aa = D[~D.comonomer_name.isin({"APheOMe", "ATrpOMe", "ANalOMe", "AProOMe"})].reset_index(drop=True)
    sens = []
    for label, cols in [("all physics features", PHYS),
                        ("without dTPSA", [c for c in PHYS if c != "frac_x_dTPSA"]),
                        ("without dTPSA, dHBD", [c for c in PHYS if c not in ("frac_x_dTPSA", "frac_x_dHBD")])]:
        sens.append({"features": label, "group_MAE_core": gmae(cols, D), "group_MAE_no_amino_acid_esters": gmae(cols, no_aa)})
    sens = pd.DataFrame(sens)
    sens.round(2).to_csv(RES / "sensitivity_physics_features.csv", index=False)
    r_lt = float(D.frac_x_dlogP.corr(D.frac_x_dTPSA))

    # Figures
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 4), sharex=True)
    for a, col, lab, color in [(ax[0], "pred_heldout_at_0.10", "Predicted LCST at 10 mol% (C)\n(comonomer held out)", BLUE),
                               (ax[1], "measured_slope_per_10mol", "Measured LCST shift per 10 mol% (C)", ORANGE)]:
        a.scatter(t1.logP, t1[col], s=36, color=color, edgecolor="white", linewidth=1.5, zorder=3)
        for _, r in t1.iterrows():
            a.annotate(r.comonomer_name, (r.logP, r[col]), fontsize=7, color=INK2,
                       xytext=LABEL_OFFSET.get(r.comonomer_name, (3, 3)), textcoords="offset points")
        a.set_xlabel("Comonomer Crippen logP")
        a.set_ylabel(lab)
    ax[0].set_title(f"Model: rho = {rho['pred_heldout_vs_logP'].statistic:.2f}", fontsize=9, color=INK, loc="left")
    ax[1].set_title(f"Data: rho = {rho['measured_slope_vs_logP'].statistic:.2f}", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "trend1_logp.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(10, 6))
    for a, c in zip(axes.flat, SERIES):
        g = D[D.comonomer_name == c]
        a.scatter(g.frac * 100, g.lcst_C, s=30, color=BLUE, edgecolor="white", linewidth=1.2, zorder=3, label="measured")
        gx, gy = curves[c]
        a.plot(gx * 100, gy, color=ORANGE, lw=2, label="predicted (comonomer held out)")
        a.set_title(f"{c}  (n={len(g)}, {g.source.nunique()} source{'s' if g.source.nunique() > 1 else ''})",
                    fontsize=9, color=INK, loc="left")
        a.set_xlabel("Comonomer (mol%)")
        a.set_ylabel("LCST (C)")
    axes.flat[0].legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "trend2_series.png", dpi=200)
    plt.close(fig)

    fig, a = plt.subplots(figsize=(5.5, 3))
    y = np.arange(len(ct))
    a.barh(y, ct.std_coef, color=BLUE, height=0.55, zorder=3)
    a.errorbar(ct.std_coef, y, xerr=[ct.std_coef - ct.boot_q05, ct.boot_q95 - ct.std_coef],
               fmt="none", ecolor=INK2, elinewidth=1, capsize=3, zorder=4)
    a.set_yticks(y, ct.feature)
    a.axvline(0, color=INK2, lw=0.8)
    a.set_xlabel("Standardized coefficient (C per SD)")
    a.set_title("Ridge, physics features; bars = 90% comonomer bootstrap", fontsize=9, color=INK, loc="left")
    a.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(FIG / "coefficients_ridge_physics.png", dpi=200)
    plt.close(fig)

    lp = pd.read_csv(RES / "predictions_logo.csv")
    lp = lp[lp.subset == "core_uncharged"]
    fig, a = plt.subplots(figsize=(4.6, 4.4))
    lim = [lp.lcst_C.min() - 3, lp.lcst_C.max() + 3]
    a.plot(lim, lim, color=INK2, lw=0.8, zorder=1)
    for mname, color, lab in [("baseline_linear_frac", ORANGE, "linear in fraction"), ("ridge_physics", BLUE, "ridge, physics features")]:
        q = lp[lp.model == mname]
        mae = np.abs(q.pred - q.lcst_C).mean()
        a.scatter(q.lcst_C, q.pred, s=22, color=color, edgecolor="white", linewidth=1, zorder=3, label=f"{lab} (MAE {mae:.1f})")
    a.set_xlim(lim); a.set_ylim(lim)
    a.set_xlabel("Measured LCST (C)"); a.set_ylabel("Predicted LCST (C), comonomer held out")
    a.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "parity_logo.png", dpi=200)
    plt.close(fig)

    print("Trend 1 (Spearman, n = %d comonomers, fraction fixed at %.2f):" % (len(t1), FIX_FRAC))
    for k, v in rho.items():
        print(f"  {k:32s} rho = {v.statistic:6.2f}  p = {v.pvalue:.3f}")
    print("\n", t1.sort_values("logP").round(2).to_string(index=False))
    print("\nTrend 2 (comonomer held out):")
    print(t2.round(2).to_string(index=False))
    print("\nRidge (physics) standardized coefficients:")
    print(ct.round(2).to_string(index=False))
    print(f"\ncorrelation frac_x_dlogP vs frac_x_dTPSA across rows: {r_lt:.2f}")
    print(sens.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
