"""Shared download and parsing helpers for NHANES files."""
import pathlib
import ssl
import urllib.request

import pandas as pd

try:  # python.org builds on macOS ship without root certificates
    import certifi
    CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    CTX = ssl.create_default_context()

RAW = pathlib.Path("data/raw")
XPT_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{name}.xpt"
MORT_URL = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/{name}"
MORT_COLS = [(0, 6), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21), (42, 45), (45, 48)]
MORT_NAMES = ["seqn", "eligstat", "mortstat", "ucod_leading", "diabetes", "hyperten", "permth_int", "permth_exm"]


def fetch(url: str, dest: pathlib.Path) -> pathlib.Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, context=CTX, timeout=120) as r:
            dest.write_bytes(r.read())
        print("downloaded", dest.name)
    return dest


def read_xpt(year: int, name: str, cols: list[str]) -> pd.DataFrame:
    path = fetch(XPT_URL.format(year=year, name=name), RAW / f"{name}.xpt")
    df = pd.read_sas(path, format="xport")
    df.columns = df.columns.str.lower()
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{name}: missing columns {missing}")
    df = df[["seqn"] + cols].copy()
    df["seqn"] = df["seqn"].astype(int)
    return df


def read_mort(years: str) -> pd.DataFrame:
    name = f"NHANES_{years.replace('-', '_')}_MORT_2019_PUBLIC.dat"
    path = fetch(MORT_URL.format(name=name), RAW / name)
    df = pd.read_fwf(path, colspecs=MORT_COLS, names=MORT_NAMES, na_values=["."])
    df["seqn"] = df["seqn"].astype(int)
    return df
