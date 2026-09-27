from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s5_trends import INK, INK2, GRID, BLUE, ORANGE

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures"

LABELS = {
    "baseline_const32": "32 C for everything",
    "baseline_linear_frac": "linear in mole fraction",
    "ridge_frac_logP": "ridge: fraction x logP",
    "ridge_physics": "ridge: fraction x (logP, TPSA, HBD)",
    "ridge_rdkit": "ridge: 8 RDKit descriptors",
    "ridge_alva_top10": "ridge: top-10 alvaDesc",
    "rf_rdkit": "random forest: RDKit",
    "gbr_rdkit": "gradient boosting: RDKit",
    "rf_alva_top10": "random forest: top-10 alvaDesc",
}


def main():
    m = pd.read_csv(RES / "metrics.csv")
    m = m[m.subset == "core_uncharged"].pivot(index="model", columns="scheme", values="MAE")
    m = m.sort_values("group_5fold", ascending=False)
    y = np.arange(len(m))

    fig, a = plt.subplots(figsize=(7, 3.8))
    for i, (name, r) in enumerate(m.iterrows()):
        a.plot([r.random_5fold, r.group_5fold], [i, i], color=GRID, lw=3, zorder=1, solid_capstyle="round")
    a.scatter(m.random_5fold, y, s=48, color=BLUE, edgecolor="white", linewidth=1.5, zorder=3, label="random 5-fold")
    a.scatter(m.group_5fold, y, s=48, color=ORANGE, edgecolor="white", linewidth=1.5, zorder=3, label="grouped 5-fold (comonomer held out)")
    for i, (name, r) in enumerate(m.iterrows()):
        a.annotate(f"{r.group_5fold:.1f}", (r.group_5fold, i), xytext=(7, -3), textcoords="offset points", fontsize=7, color=INK2)
    a.set_yticks(y, [LABELS[n] for n in m.index])
    a.set_xlabel("Mean absolute error (C), uncharged comonomers")
    a.set_xlim(0, m.group_5fold.max() * 1.15)
    a.grid(axis="y", visible=False)
    a.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "results_summary.png", dpi=200)
    plt.close(fig)

    d = pd.read_csv(ROOT / "data" / "raw" / "lcst.csv")
    t = d[(d.role == "train") & (d.comonomer_name != "none")]
    g = t.groupby("comonomer_name").agg(rows=("entry_id", "size"), sources=("source", "nunique"),
                                        lo=("lcst_C", "min"), hi=("lcst_C", "max")).sort_values("rows")
    fig, ax = plt.subplots(1, 2, figsize=(8, 5), sharey=True, gridspec_kw={"width_ratios": [1, 1.4]})
    yy = np.arange(len(g))
    ax[0].barh(yy, g.rows, color=BLUE, height=0.6, zorder=3)
    for i, (n, r) in enumerate(g.iterrows()):
        ax[0].annotate(f"{int(r.rows)}" + (f"  ({int(r.sources)} sources)" if r.sources > 1 else ""), (r.rows, i),
                       xytext=(3, -3), textcoords="offset points", fontsize=7, color=INK2)
    ax[0].set_yticks(yy, g.index)
    ax[0].set_xlabel("Training rows")
    ax[0].set_xlim(0, g.rows.max() * 1.6)
    ax[0].grid(axis="y", visible=False)
    for i, (n, r) in enumerate(g.iterrows()):
        ax[1].plot([r.lo, r.hi], [i, i], color=BLUE, lw=2, solid_capstyle="round", zorder=3)
        ax[1].scatter([r.lo, r.hi], [i, i], s=14, color=BLUE, zorder=4)
    ax[1].axvline(32, color=INK2, lw=0.8, ls="--")
    ax[1].text(32.5, len(g) - 0.6, "PNIPAM ~32 C", fontsize=7, color=INK2)
    ax[1].set_xlabel("Measured LCST range (C)")
    ax[1].grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(FIG / "dataset_overview.png", dpi=200)
    plt.close(fig)
    print("wrote results_summary.png, dataset_overview.png")


if __name__ == "__main__":
    main()
