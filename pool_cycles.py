"""Week 3, part 1: pool five NHANES cycles that measured CRP and have mortality follow-up to 2019.

2011-2014 are skipped because CRP was not measured in those cycles.
Writes data/pooled.csv with PhenoAge and age acceleration for every adult.
"""
import numpy as np
import pandas as pd

from nhanes_io import read_mort, read_xpt
from phenoage import phenoage

# cycle -> (file suffix, year used in the CDC URL, CRP file, CRP column, factor to mg/dL)
CYCLES = {
    "2005-2006": ("D", 2005, "CRP_D", "lbxcrp", 1.0),
    "2007-2008": ("E", 2007, "CRP_E", "lbxcrp", 1.0),
    "2009-2010": ("F", 2009, "CRP_F", "lbxcrp", 1.0),
    "2015-2016": ("I", 2015, "HSCRP_I", "lbxhscrp", 0.1),
    "2017-2018": ("J", 2017, "HSCRP_J", "lbxhscrp", 0.1),
}
DEMO_COLS = ["ridstatr", "ridageyr", "riagendr", "wtmec2yr"]
CBC_COLS = ["lbxwbcsi", "lbxlypct", "lbxmcvsi", "lbxrdw"]
BIO_COLS = ["lbxsal", "lbxscr", "lbxsgl", "lbxsapsi"]
AGE_CAP = 80  # top-coding differs between cycles (85 in 2005-06, 80 later); use the common cap


def one_cycle(label: str, sfx: str, year: int, crp_file: str, crp_col: str, crp_to_mgdl: float) -> pd.DataFrame:
    d = read_xpt(year, f"DEMO_{sfx}", DEMO_COLS)
    for part in (read_xpt(year, f"CBC_{sfx}", CBC_COLS), read_xpt(year, f"BIOPRO_{sfx}", BIO_COLS),
                 read_xpt(year, crp_file, [crp_col]), read_mort(label)):
        d = d.merge(part, on="seqn", how="left")
    n = {"all": len(d)}
    d = d[(d.ridageyr >= 20) & (d.ridstatr == 2)]
    n["adults_exam"] = len(d)
    d = d[(d.eligstat == 1) & d.permth_exm.notna()]
    n["mort_eligible"] = len(d)
    d = d.dropna(subset=CBC_COLS + BIO_COLS + [crp_col])
    d = d[d[crp_col] > 0]
    n["complete"] = len(d)
    out = pd.DataFrame({
        "seqn": d.seqn, "cycle": label, "age": d.ridageyr.clip(upper=AGE_CAP), "sex": d.riagendr,
        "weight": d.wtmec2yr / len(CYCLES),  # NHANES rule for pooling 2-year MEC weights
        "albumin_gl": d.lbxsal * 10, "creat_umol": d.lbxscr * 88.4, "glucose_mmol": d.lbxsgl / 18.016,
        "alp": d.lbxsapsi, "wbc": d.lbxwbcsi, "lymph": d.lbxlypct, "mcv": d.lbxmcvsi, "rdw": d.lbxrdw,
        "lncrp": np.log(d[crp_col] * crp_to_mgdl),
        "mortstat": d.mortstat.astype(int), "permth_exm": d.permth_exm,
    })
    n["deaths"] = int(out.mortstat.sum())
    print(label, n)
    return out, dict(cycle=label, **n)


if __name__ == "__main__":
    parts, flow = zip(*(one_cycle(k, *v) for k, v in CYCLES.items()))
    d = phenoage(pd.concat(parts, ignore_index=True))
    assert np.isfinite(d.phenoage).all()
    assert d.albumin_gl.between(15, 70).all() and d.lncrp.between(-8, 5).all()
    slope, icpt = np.polyfit(d.age, d.phenoage, 1)
    d["accel"] = d.phenoage - (icpt + slope * d.age)
    d.to_csv("data/pooled.csv", index=False)
    flow = pd.DataFrame(flow)
    flow.loc[len(flow)] = ["total"] + flow.drop(columns="cycle").sum().tolist()
    flow.to_csv("results/week3_dataflow.csv", index=False)
    print(flow.to_string(index=False))
