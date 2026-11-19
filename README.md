# Blutalter: biological age from routine blood tests

How old is your blood? This project computes **PhenoAge**, a biological-age score built from nine routine blood markers plus age, for adults in the US National Health and Nutrition Examination Survey (NHANES), and tests whether being "older than your age" predicts mortality.

Status: **in progress** (week 2 of 4 done: PhenoAge implemented and checked; pooled Cox models next).

## Plan

| Step | What | Status |
| --- | --- | --- |
| 1 | Data pipeline: NHANES 2015–2016 + linked mortality into DuckDB, one analysis table | done |
| 2 | PhenoAge implementation, checked against the R package BioAge | done |
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
python phenoage.py         # writes data/phenoage.csv, figures/, results/
Rscript r/check_bioage.R   # cross-check against the BioAge R package
```

## Data flow: who is dropped and why

| Step | Rule | People left |
| --- | --- | --- |
| All participants | `demo` | 9,971 |
| Adults with exam | age >= 20, `ridstatr = 2` | 5,474 (of 5,719 adults) |
| Eligible for mortality follow-up | `eligstat = 1`, follow-up months present | 5,459 |
| Complete markers | all 9 markers present, hs-CRP > 0 | 5,087 (222 deaths, 4.4%) |

*Numbers are from the actual run.*

## Results so far (NHANES 2015–2016, n = 5,087, 222 deaths)

![PhenoAge vs chronological age](figures/phenoage_vs_age.png)

| Check | Result |
| --- | --- |
| Correlation of PhenoAge with age | 0.90 |
| Agreement with BioAge R package (original Levine coefficients) | max difference 0.0008 years on 5,077 people |
| C-index for mortality, age alone | 0.837 |
| C-index for mortality, PhenoAge | 0.866 |
| Hazard ratio per 1 SD of age acceleration, adjusted for age and sex | 1.50 (95% CI 1.41–1.59) |

Notes:
- The BioAge formula returns infinity for 10 people whose predicted mortality rounds to 1 (PhenoAge 129–163 years). The Python version computes the cumulative hazard directly, so these rows stay finite. They are extreme values and are kept for now; step 3 will check whether they change the Cox results.
- The C-index is in-sample here. Gains from PhenoAge over age are small and will be compared properly with cross-validation in step 4.

## Repository layout

- `load_nhanes.py`: download and load raw files
- `sql/`: queries, numbered in run order
- `build_analysis.py`: builds the analysis table with unit checks
- `phenoage.py`: PhenoAge, age acceleration, C-index and Cox model
- `r/check_bioage.R`: cross-check against the BioAge R package

## Limitations (so far)

- NHANES survey weights are not yet used, so results describe the sample, not the US population.
- One 2015–2016 cycle has only a few years of follow-up and few deaths; later steps pool cycles.

## Reference

Levine ME et al. (2018). An epigenetic biomarker of aging for lifespan and healthspan. *Aging*, 10(4), 573–591.
