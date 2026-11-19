"""Week 2: PhenoAge (Levine et al. 2018) on the analysis table, and a first mortality check.

Reads data/analysis.csv (made by build_analysis.py), writes data/phenoage.csv,
figures/phenoage_vs_age.png and results/week2_summary.csv.
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

GAMMA = 0.0076927
COEF = {  # Levine 2018 coefficients; inputs in the units of analysis.csv
    "albumin_gl": -0.0336, "creat_umol": 0.0095, "glucose_mmol": 0.1953, "lncrp": 0.0954,
    "lymph": -0.0120, "mcv": 0.0268, "rdw": 0.3306, "alp": 0.00188, "wbc": 0.0554, "age": 0.0804,
}
INTERCEPT = -19.9067
HORIZON = 120  # months


def phenoage(d: pd.DataFrame) -> pd.DataFrame:
    xb = INTERCEPT + sum(c * d[k] for k, c in COEF.items())
    # -ln(1 - M) is computed directly: for low-risk people M is ~1e-6 and 1 - M loses precision
    cum_hazard = np.exp(xb) * (np.exp(HORIZON * GAMMA) - 1) / GAMMA
    out = d.copy()
    out["xb"] = xb
    out["mort_score"] = -np.expm1(-cum_hazard)
    out["phenoage"] = 141.50225 + np.log(0.00553 * cum_hazard) / 0.090165
    return out


def test_known_value() -> None:
    """A 50-year-old with round-number labs should land near 50; guards against a typo in COEF."""
    row = pd.DataFrame([dict(age=50, albumin_gl=43, creat_umol=75, glucose_mmol=5.3, lncrp=np.log(0.2),
                             lymph=30, mcv=90, rdw=13.2, alp=70, wbc=6.5)])
    pa = float(phenoage(row)["phenoage"].iloc[0])
    assert 35 < pa < 55, pa


if __name__ == "__main__":
    test_known_value()
    for p in ("figures", "results"):
        pathlib.Path(p).mkdir(exist_ok=True)
    d = phenoage(pd.read_csv("data/analysis.csv"))
    assert np.isfinite(d["phenoage"]).all()
    # age acceleration = residual of PhenoAge regressed on chronological age
    slope, icpt = np.polyfit(d["age"], d["phenoage"], 1)
    d["accel"] = d["phenoage"] - (icpt + slope * d["age"])
    d.to_csv("data/phenoage.csv", index=False)

    r = np.corrcoef(d["age"], d["phenoage"])[0, 1]
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(d["age"], d["phenoage"], s=3, alpha=0.3)
    ax.plot([20, 85], [20, 85], "k--", lw=1)
    ax.set(xlabel="Chronological age (years)", ylabel="PhenoAge (years)",
           title=f"NHANES 2015-16, n={len(d)}, r={r:.2f}")
    fig.tight_layout()
    fig.savefig("figures/phenoage_vs_age.png", dpi=150)

    # mortality: C-index and age-adjusted hazard ratio per 1 SD of acceleration
    ev = d["mortstat"].astype(bool)
    t = d["permth_exm"]
    c_age = concordance_index(t, -d["age"], ev)  # lifelines wants higher score = longer survival
    c_pheno = concordance_index(t, -d["phenoage"], ev)
    d["accel_sd"] = d["accel"] / d["accel"].std()
    cph = CoxPHFitter().fit(d[["permth_exm", "mortstat", "accel_sd", "age", "sex"]],
                            duration_col="permth_exm", event_col="mortstat")
    hr = cph.summary.loc["accel_sd"]
    res = pd.DataFrame([dict(n=len(d), deaths=int(ev.sum()), corr_age_phenoage=r, slope=slope,
                             c_index_age=c_age, c_index_phenoage=c_pheno, hr_per_sd=np.exp(hr["coef"]),
                             hr_lo=hr["exp(coef) lower 95%"], hr_hi=hr["exp(coef) upper 95%"], p=hr["p"])])
    res.to_csv("results/week2_summary.csv", index=False)
    print(res.T.round(4).to_string())
