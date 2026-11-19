"""Week 4: can models trained on these cycles beat the published PhenoAge at predicting death?

Same 9 markers + age + sex for every model. One train/test split (70/30, stratified by cycle and
death); all tuning happens by 5-fold cross-validation inside the training set; the test set is
touched once. Writes results/week4_*.csv and figures/week4_*.png.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index, k_fold_cross_validation
from sklearn.model_selection import train_test_split

SEED = 42
T, E = "permth_exm", "mortstat"
d = pd.read_csv("data/pooled.csv")
d["female"] = (d.sex == 2).astype(int)
for c in ("creat_umol", "glucose_mmol", "alp", "wbc"):  # right-skewed markers
    d[f"log_{c}"] = np.log(d[c])
FEATS = ["age", "female", "albumin_gl", "log_creat_umol", "log_glucose_mmol", "lncrp",
         "lymph", "mcv", "rdw", "log_alp", "log_wbc"]
LABELS = {"age": "Age", "female": "Female", "albumin_gl": "Albumin", "log_creat_umol": "Creatinine",
          "log_glucose_mmol": "Glucose", "lncrp": "CRP", "lymph": "Lymphocyte %", "mcv": "MCV",
          "rdw": "RDW", "log_alp": "Alkaline phosphatase", "log_wbc": "White cells"}

tr, te = train_test_split(d, test_size=0.3, random_state=SEED, stratify=d.cycle + "_" + d[E].astype(str))
assert not set(tr.seqn.astype(str) + tr.cycle) & set(te.seqn.astype(str) + te.cycle)  # no person in both
mu, sd = tr[FEATS].mean(), tr[FEATS].std()
Ztr, Zte = (tr[FEATS] - mu) / sd, (te[FEATS] - mu) / sd


def cidx(df, score):
    return concordance_index(df[T], -np.asarray(score), df[E])


scores = {}
scores["Age only"] = te.age.values
cox_as = CoxPHFitter().fit(tr[[T, E, "age", "female"]], T, E)
scores["Age + sex (Cox)"] = cox_as.predict_log_partial_hazard(te[["age", "female"]]).values
scores["PhenoAge (published, not refitted)"] = te.phenoage.values

# Elastic-net Cox: penalty chosen by 5-fold CV C-index on the training set
tr_z = pd.concat([tr[[T, E]].reset_index(drop=True), Ztr.reset_index(drop=True)], axis=1)
cv = {}
for pen in (0.0003, 0.001, 0.003, 0.01, 0.03, 0.1):
    cv[pen] = float(np.mean(k_fold_cross_validation(CoxPHFitter(penalizer=pen, l1_ratio=0.5), tr_z, T, E,
                                                    k=5, scoring_method="concordance_index", seed=SEED)))
best_pen = max(cv, key=cv.get)
enet = CoxPHFitter(penalizer=best_pen, l1_ratio=0.5).fit(tr_z, T, E)
scores["Elastic-net Cox (trained here)"] = enet.predict_log_partial_hazard(Zte.reset_index(drop=True)).values

# Gradient boosting with a Cox objective: depth and number of trees chosen by 5-fold CV
lab = np.where(tr[E] == 1, tr[T] + 0.5, -(tr[T] + 0.5))  # xgboost: negative time = censored
dtr = xgb.DMatrix(tr[FEATS], label=lab, feature_names=FEATS)
best = None
for depth in (2, 3, 4):
    p = dict(objective="survival:cox", eval_metric="cox-nloglik", eta=0.03, max_depth=depth,
             subsample=0.8, colsample_bytree=0.8, min_child_weight=10, seed=SEED)
    r = xgb.cv(p, dtr, num_boost_round=3000, nfold=5, early_stopping_rounds=100, seed=SEED)
    cand = (r["test-cox-nloglik-mean"].min(), depth, len(r), p)
    best = cand if best is None or cand[0] < best[0] else best
_, depth, n_rounds, params = best
booster = xgb.train(params, dtr, num_boost_round=n_rounds)
dte = xgb.DMatrix(te[FEATS], feature_names=FEATS)
scores["Gradient boosting (XGBoost Cox)"] = booster.predict(dte, output_margin=True)

# Test-set C-index with bootstrap 95% CIs, and differences vs PhenoAge
rng = np.random.default_rng(SEED)
names = list(scores)
point = {k: cidx(te, v) for k, v in scores.items()}
boot = {k: [] for k in names}
idx = np.arange(len(te))
for _ in range(500):
    b = rng.choice(idx, len(idx), replace=True)
    sub = te.iloc[b]
    for k in names:
        boot[k].append(cidx(sub, scores[k][b]))
ref = np.array(boot["PhenoAge (published, not refitted)"])
rows = []
for k in names:
    bk = np.array(boot[k])
    diff = bk - ref
    rows.append(dict(model=k, c_index=point[k], lo=np.percentile(bk, 2.5), hi=np.percentile(bk, 97.5),
                     diff_vs_phenoage=point[k] - point["PhenoAge (published, not refitted)"],
                     diff_lo=np.percentile(diff, 2.5), diff_hi=np.percentile(diff, 97.5)))
res = pd.DataFrame(rows)
res.to_csv("results/week4_cindex.csv", index=False)

# Explanations: SHAP values from xgboost (exact TreeSHAP), elastic-net coefficients
contrib = booster.predict(dte, pred_contribs=True)[:, :-1]
shap = pd.DataFrame(contrib, columns=FEATS)
imp = shap.abs().mean().sort_values(ascending=False)
imp.rename("mean_abs_shap").to_frame().assign(feature=lambda x: x.index.map(LABELS)) \
   .to_csv("results/week4_shap.csv")
coef = enet.params_.rename("coef_per_sd").to_frame()
coef["hr_per_sd"] = np.exp(coef.coef_per_sd)
coef.assign(feature=coef.index.map(LABELS)).sort_values("coef_per_sd").to_csv("results/week4_enet_coefs.csv")
meta = dict(n_train=len(tr), n_test=len(te), deaths_train=int(tr[E].sum()), deaths_test=int(te[E].sum()),
            enet_cv=cv, enet_penalizer=best_pen, xgb_depth=depth, xgb_rounds=n_rounds)
json.dump(meta, open("results/week4_meta.json", "w"), indent=2)

# Figure 1: test C-index
fig, ax = plt.subplots(figsize=(7, 3.6))
y = np.arange(len(res))[::-1]
ax.errorbar(res.c_index, y, xerr=[res.c_index - res.lo, res.hi - res.c_index], fmt="o", capsize=3)
ax.set_yticks(y, res.model)
ax.set_xlabel("C-index on held-out test set (95% bootstrap CI)")
ax.set_title(f"Predicting death: test n={len(te):,}, {meta['deaths_test']} deaths")
fig.tight_layout()
fig.savefig("figures/week4_cindex.png", dpi=90)

# Figure 2: SHAP importance + beeswarm-style spread for the gradient-boosting model
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw=dict(width_ratios=[1, 1.4]))
order = imp.index[::-1]
axes[0].barh([LABELS[f] for f in order], imp[order])
axes[0].set(xlabel="Mean |SHAP| (log hazard)", title="Which markers drive the prediction")
for i, f in enumerate(order):
    v = te[f].values
    col = (v - np.percentile(v, 5)) / (np.percentile(v, 95) - np.percentile(v, 5) + 1e-9)
    jitter = rng.uniform(-0.3, 0.3, len(v))
    axes[1].scatter(shap[f], i + jitter, c=np.clip(col, 0, 1), cmap="coolwarm", s=2, alpha=0.4)
axes[1].set_yticks(range(len(order)), [LABELS[f] for f in order])
axes[1].axvline(0, color="k", lw=0.5)
axes[1].set(xlabel="SHAP value (red = high marker value)", title="Direction of each effect")
fig.tight_layout()
fig.savefig("figures/week4_shap.png", dpi=90)

pd.set_option("display.width", 160)
print(json.dumps(meta))
print(res.round(4).to_string(index=False))
print(imp.round(3).to_string())
print(coef.round(3).to_string())
