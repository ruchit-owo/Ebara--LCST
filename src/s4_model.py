from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold, KFold, LeaveOneGroupOut
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
MODEL_TABLE = PROC / "model_table.csv"
FEATURES = PROC / "features_clean.csv"

SEED = 0
N_REPEATS = 5
K_ALVA = 10
ALPHAS = np.logspace(-3, 3, 13)
CHARGED = {"AAc", "DMAPAAm", "DMAPMAm"}
COND = ["conc_wt", "is_buffer"]


def load():
    mt = pd.read_csv(MODEL_TABLE)
    mt = mt[mt.role == "train"].reset_index(drop=True)
    feat = pd.read_csv(FEATURES).set_index("comonomer_name")
    nipam = feat.loc["NIPAM"]
    raw = feat.reindex(mt.comonomer_name.replace("none", "NIPAM")).reset_index(drop=True)
    delta = (raw - nipam).mul(mt.frac.values, axis=0)
    delta.columns = ["fx_" + c for c in delta.columns]
    mt["group"] = mt.comonomer_name
    return mt, raw, delta


def ridge():
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), RidgeCV(alphas=ALPHAS))


def forest():
    return make_pipeline(SimpleImputer(strategy="median"),
                         RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=SEED, n_jobs=-1))


def boosting():
    return make_pipeline(SimpleImputer(strategy="median"),
                         GradientBoostingRegressor(n_estimators=300, max_depth=2, learning_rate=0.05,
                                                   subsample=0.8, random_state=SEED))


def top_k(X, y, k):
    Xc = X - X.mean()
    yc = y - y.mean()
    r = (Xc.mul(yc, axis=0).sum() / np.sqrt((Xc ** 2).sum() * (yc ** 2).sum())).abs().fillna(0)
    return r.sort_values(ascending=False).index[:k].tolist()


RD_FX = ["fx_rd_logP", "fx_rd_MR", "fx_rd_TPSA", "fx_rd_HBD", "fx_rd_HBA", "fx_rd_RotB", "fx_rd_AromRings", "fx_rd_FracCSP3"]
RD_RAW = [c[3:] for c in RD_FX]


def make_models(mt, raw, delta):
    D = pd.concat([mt, raw, delta], axis=1)
    alva_fx = [c for c in delta.columns if c.startswith("fx_al_")]

    def fixed(cols, est):
        def fit_predict(tr, te):
            m = est()
            m.fit(D.loc[tr, cols], D.loc[tr, "lcst_C"])
            return m.predict(D.loc[te, cols])
        return fit_predict

    def selected(pool, k, est, raw_too=False):
        def fit_predict(tr, te):
            sel = top_k(D.loc[tr, pool], D.loc[tr, "lcst_C"], k)
            cols = ["frac"] + COND + (sel if not raw_too else [c[3:] for c in sel])
            m = est()
            m.fit(D.loc[tr, cols], D.loc[tr, "lcst_C"])
            return m.predict(D.loc[te, cols])
        return fit_predict

    def const(tr, te):
        return np.full(len(te), 32.0)

    return D, {
        "baseline_const32": const,
        "baseline_linear_frac": fixed(["frac"], LinearRegression),
        "ridge_frac_logP": fixed(["frac", "frac_x_dlogP"], ridge),
        "ridge_physics": fixed(["frac", "frac_x_dlogP", "frac_x_dTPSA", "frac_x_dHBD", "frac_x_charged"] + COND, ridge),
        "ridge_rdkit": fixed(["frac"] + RD_FX + COND, ridge),
        f"ridge_alva_top{K_ALVA}": selected(alva_fx, K_ALVA, ridge),
        "rf_rdkit": fixed(["frac"] + RD_RAW + COND, forest),
        "gbr_rdkit": fixed(["frac"] + RD_RAW + COND, boosting),
        f"rf_alva_top{K_ALVA}": selected(alva_fx, K_ALVA, forest, raw_too=True),
    }


def splits(scheme, D, rep):
    idx = np.arange(len(D))
    if scheme == "random_5fold":
        return KFold(5, shuffle=True, random_state=SEED + rep).split(idx)
    if scheme == "group_5fold":
        return GroupKFold(5, shuffle=True, random_state=SEED + rep).split(idx, groups=D.group)
    return LeaveOneGroupOut().split(idx, groups=D.group)


def metrics(y, p):
    return mean_absolute_error(y, p), np.sqrt(mean_squared_error(y, p)), r2_score(y, p)


def main():
    mt, raw, delta = load()
    rows, logo_preds = [], []
    subsets = {"core_uncharged": ~mt.comonomer_name.isin(CHARGED), "all": pd.Series(True, index=mt.index)}
    for sub_name, mask in subsets.items():
        D, models = make_models(mt[mask].reset_index(drop=True), raw[mask].reset_index(drop=True),
                                delta[mask].reset_index(drop=True))
        y = D.lcst_C.values
        for scheme in ["random_5fold", "group_5fold", "logo"]:
            reps = 1 if scheme == "logo" else N_REPEATS
            for name, fp in models.items():
                scores = []
                for rep in range(reps):
                    pred = np.zeros(len(D))
                    for tr, te in splits(scheme, D, rep):
                        pred[te] = fp(tr, te)
                    scores.append(metrics(y, pred))
                    if scheme == "logo":
                        logo_preds.append(pd.DataFrame({"subset": sub_name, "model": name, "entry_id": D.entry_id,
                                                        "comonomer_name": D.comonomer_name, "frac": D.frac,
                                                        "lcst_C": y, "pred": pred}))
                s = np.array(scores)
                rows.append({"subset": sub_name, "n_rows": len(D), "n_groups": D.group.nunique(),
                             "scheme": scheme, "model": name,
                             "MAE": s[:, 0].mean(), "MAE_sd": s[:, 0].std(),
                             "RMSE": s[:, 1].mean(), "R2": s[:, 2].mean()})
                print(f"{sub_name:15s} {scheme:13s} {name:22s} MAE {s[:,0].mean():5.2f}  RMSE {s[:,1].mean():5.2f}  R2 {s[:,2].mean():6.2f}")

    RES.mkdir(exist_ok=True)
    out = pd.DataFrame(rows)
    out.round(3).to_csv(RES / "metrics.csv", index=False)
    pd.concat(logo_preds).round(4).to_csv(RES / "predictions_logo.csv", index=False)

    piv = out[out.subset == "core_uncharged"].pivot(index="model", columns="scheme", values="MAE")
    piv["gap_group_minus_random"] = piv["group_5fold"] - piv["random_5fold"]
    print("\ncore_uncharged, MAE (C):")
    print(piv[["random_5fold", "group_5fold", "logo", "gap_group_minus_random"]].round(2).sort_values("group_5fold").to_string())


if __name__ == "__main__":
    main()
