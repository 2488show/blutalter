"""Run the SQL in sql/ in order, check units, and write analysis.csv."""
import pathlib

import duckdb

con = duckdb.connect("data/nhanes.duckdb")
for f in ["01_adults_raw.sql", "04_analysis_table.sql"]:
    con.execute((pathlib.Path("sql") / f).read_text())


def gdl_to_gl(x):            return x * 10
def mgdl_to_umol_creat(x):   return x * 88.4
def mgdl_to_mmol_glucose(x): return x / 18.016


assert abs(mgdl_to_umol_creat(1.0) - 88.4) < 1e-9
assert abs(mgdl_to_mmol_glucose(100) - 5.55) < 0.01
assert gdl_to_gl(4.5) == 45

a = con.execute("SELECT * FROM analysis").df()
assert a.isna().sum().sum() == 0, "analysis table must have no missing values"
assert a.albumin_gl.between(15, 70).all(), "albumin out of usual range: check units"
a.to_csv("data/analysis.csv", index=False)
print(a.shape)
