"""Week 3, part 2: Cox models of PhenoAge acceleration vs all-cause mortality in the pooled cycles.

Reads data/pooled.csv; writes results/week3_cox.csv and figures/km_accel_quartiles.png.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import concordance_index

d = pd.read_csv("data/pooled.csv")
# The acceleration has a long right tail (a few people with dialysis-level creatinine or extreme WBC
# reach PhenoAge > 150). Main models use it winsorised at the 0.5th/99.5th percentiles.
lo_q, hi_q = d.accel.quantile([0.005, 0.995])
d["accel_w"] = d.accel.clip(lo_q, hi_q)
sd = d.accel_w.std()
d["accel_sd"] = d.accel_w / sd
d["accel_raw_sd"] = d.accel / d.accel.std()
d["female"] = (d.sex == 2).astype(int)
T, E = "permth_exm", "mortstat"
rows = []


def fit(name, data, covars, strata=None, weights=None):
    cols = [T, E] + covars + (strata or []) + ([weights] if weights else [])
    cph = CoxPHFitter().fit(data[cols], duration_col=T, event_col=E, strata=strata,
                            weights_col=weights, robust=bool(weights))
    s = cph.summary.loc["accel_sd"]
    rows.append(dict(model=name, n=len(data), deaths=int(data[E].sum()), hr_per_sd=s["exp(coef)"],
                     lo=s["exp(coef) lower 95%"], hi=s["exp(coef) upper 95%"], p=s["p"]))
    return cph


base = ["accel_sd", "age", "female"]
m1 = fit("main: age + sex, stratified by cycle", d, base, strata=["cycle"])
fit("survey weights (robust SE; no strata/PSU)", d, base, strata=["cycle"], weights="weight")
raw = d.drop(columns="accel_sd").rename(columns={"accel_raw_sd": "accel_sd"})
fit("not winsorised (raw acceleration)", raw, base, strata=["cycle"])
fit("drop raw |accel| > 5 SD", raw[raw.accel_sd.abs() <= 5], base, strata=["cycle"])
fit("drop deaths in first 2 years", d[~((d[E] == 1) & (d[T] < 24))], base, strata=["cycle"])
sub = [(d.age < 60, "age < 60", base), (d.age >= 60, "age >= 60", base),
       (d.female == 0, "men", ["accel_sd", "age"]), (d.female == 1, "women", ["accel_sd", "age"])]
for g, lab, cov in sub:
    fit(f"subgroup: {lab}", d[g], cov, strata=["cycle"])
for c, g in d.groupby("cycle"):
    fit(f"cycle {c}", g, base)

res = pd.DataFrame(rows)
res["hr_per_5y"] = np.exp(np.log(res.hr_per_sd) / sd * 5)
ph = proportional_hazard_test(m1, d[[T, E, "accel_sd", "age", "female", "cycle"]], time_transform="rank")
res.attrs["ph_p_accel"] = float(ph.summary.loc["accel_sd", "p"])
res.to_csv("results/week3_cox.csv", index=False)

c_age = concordance_index(d[T], -d.age, d[E])
c_pheno = concordance_index(d[T], -d.phenoage, d[E])
pd.DataFrame([dict(sd_accel_years=sd, winsor_lo=lo_q, winsor_hi=hi_q, n_winsorised=int((d.accel != d.accel_w).sum()), c_index_age=c_age, c_index_phenoage=c_pheno,
                   ph_test_p_accel=res.attrs["ph_p_accel"])]).to_csv("results/week3_summary.csv", index=False)

d["q"] = pd.qcut(d.accel, 4, labels=["Q1 (youngest blood)", "Q2", "Q3", "Q4 (oldest blood)"])
fig, ax = plt.subplots(figsize=(6, 4.5))
for q, g in d.groupby("q", observed=True):
    KaplanMeierFitter().fit(g[T] / 12, g[E], label=str(q)).plot_survival_function(ax=ax, ci_show=False)
ax.set(xlabel="Years since exam", ylabel="Survival", title=f"Survival by PhenoAge acceleration quartile (n={len(d):,})")
fig.tight_layout()
fig.savefig("figures/km_accel_quartiles.png", dpi=90)

pd.set_option("display.width", 160)
print(res.round(3).to_string(index=False))
print(f"SD of acceleration: {sd:.2f} y | C-index age {c_age:.3f}, PhenoAge {c_pheno:.3f} | PH test p(accel) {res.attrs['ph_p_accel']:.3g}")
