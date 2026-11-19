"""Download NHANES 2015-2016 files and load them into DuckDB.

Files that already exist in data/raw/ are not downloaded again. If the CDC
server refuses the download (HTTP 403), download the files by hand from the
NHANES 2015-2016 page and put them in data/raw/.
"""
import pathlib
import urllib.request

import duckdb
import pandas as pd

RAW = pathlib.Path("data/raw")
RAW.mkdir(parents=True, exist_ok=True)

BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/"
XPT = {"demo": "DEMO_I.xpt", "cbc": "CBC_I.xpt", "biopro": "BIOPRO_I.xpt", "hscrp": "HSCRP_I.xpt"}
MORT_NAME = "NHANES_2015_2016_MORT_2019_PUBLIC.dat"
MORT_URL = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/" + MORT_NAME


def fetch(url: str, dest: pathlib.Path) -> None:
    if dest.exists():
        return
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        f.write(r.read())
    print("downloaded", dest.name)


con = duckdb.connect("data/nhanes.duckdb")

for table, fname in XPT.items():
    fetch(BASE + fname, RAW / fname)
    df = pd.read_sas(RAW / fname, format="xport")
    df.columns = df.columns.str.lower()
    df["seqn"] = df["seqn"].astype(int)
    con.register("df_tmp", df)
    con.execute(f"CREATE OR REPLACE TABLE {table} AS SELECT * FROM df_tmp")
    con.unregister("df_tmp")

fetch(MORT_URL, RAW / MORT_NAME)
# Fixed-width file. Positions from the NCHS file layout (1-based); pandas is 0-based, end-exclusive.
mort = pd.read_fwf(
    RAW / MORT_NAME,
    colspecs=[(0, 6), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21), (42, 45), (45, 48)],
    names=["seqn", "eligstat", "mortstat", "ucod_leading", "diabetes", "hyperten", "permth_int", "permth_exm"],
    na_values=[".", ""],
)
con.register("mort_tmp", mort)
con.execute("CREATE OR REPLACE TABLE mort AS SELECT * FROM mort_tmp")

for t in ["demo", "cbc", "biopro", "hscrp", "mort"]:
    n = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
    print(f"{t:7s} {n:6d} rows")
