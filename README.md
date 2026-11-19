# Blutalter: biological age from routine blood tests

How old is your blood? This project computes **PhenoAge**, a biological-age score built from nine routine blood markers plus age, for adults in the US National Health and Nutrition Examination Survey (NHANES), and tests whether being "older than your age" predicts mortality.

Status: **in progress** (week 1 of 4: data pipeline done; PhenoAge next).

## Plan

| Step | What | Status |
| --- | --- | --- |
| 1 | Data pipeline: NHANES 2015–2016 + linked mortality into DuckDB, one analysis table | done |
| 2 | PhenoAge implementation, checked against the R package BioAge | planned |
| 3 | Pool NHANES cycles; Cox model of age acceleration vs. all-cause mortality | planned |
| 4 | Machine learning (elastic net, gradient boosting) vs. PhenoAge, compared by C-index; SHAP | planned |

## Data

Public NHANES 2015–2016 files from CDC: demographics (`DEMO_I`), complete blood count (`CBC_I`), biochemistry (`BIOPRO_I`), hs-CRP (`HSCRP_I`), and the NCHS public-use linked mortality file. No data are stored in this repository; `load_nhanes.py` downloads them.

## Reproduce

```bash
conda env create -f environment.yml
conda activate blutalter
python load_nhanes.py      # download + load into data/nhanes.duckdb
python build_analysis.py   # writes data/analysis.csv
```

## Data flow: who is dropped and why

| Step | Rule | People left |
| --- | --- | --- |
| All participants | `demo` | 9,971 |
| Adults with exam | age >= 20, `ridstatr = 2` | 5,474 (of 5,719 adults) |
| Eligible for mortality follow-up | `eligstat = 1`, follow-up months present | 5,459 |
| Complete markers | all 9 markers present, hs-CRP > 0 | 5,087 (222 deaths, 4.4%) |

*Numbers are filled in from the actual run.*

## Repository layout

- `load_nhanes.py`: download and load raw files
- `sql/`: queries, numbered in run order
- `build_analysis.py`: builds the analysis table with unit checks

## Limitations (so far)

- NHANES survey weights are not yet used, so results describe the sample, not the US population.
- One 2015–2016 cycle has only a few years of follow-up and few deaths; later steps pool cycles.

## Reference

Levine ME et al. (2018). An epigenetic biomarker of aging for lifespan and healthspan. *Aging*, 10(4), 573–591.
