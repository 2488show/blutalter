# Blutalter: biological age from routine blood tests

How old is your blood? This project computes **PhenoAge**, a biological-age score built from nine routine blood markers plus age, for adults in the US National Health and Nutrition Examination Survey (NHANES), and tests whether being "older than your age" predicts mortality.

Status: **in progress** (week 3 of 4 done: five pooled NHANES cycles, Cox models; machine-learning comparison next).

## Plan

| Step | What | Status |
| --- | --- | --- |
| 1 | Data pipeline: NHANES 2015–2016 + linked mortality into DuckDB, one analysis table | done |
| 2 | PhenoAge implementation, checked against the R package BioAge | done |
| 3 | Pool NHANES cycles; Cox model of age acceleration vs. all-cause mortality | done |
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
python pool_cycles.py      # five cycles -> data/pooled.csv
python cox.py              # Cox models -> results/week3_cox.csv, figures/km_accel_quartiles.png
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

## Week 3: five pooled cycles (n = 25,300, 2,874 deaths)

Cycles 2005–2006, 2007–2008, 2009–2010, 2015–2016 and 2017–2018 (2011–2014 measured no CRP). Mortality follow-up runs to 31 Dec 2019, so the older cycles carry most of the deaths.

| Cycle | Participants | Adults with exam | Mortality-eligible | Complete markers | Deaths |
| --- | --- | --- | --- | --- | --- |
| 2005–2006 | 10,348 | 4,773 | 4,771 | 4,427 | 874 |
| 2007–2008 | 10,149 | 5,707 | 5,701 | 5,276 | 935 |
| 2009–2010 | 10,537 | 6,059 | 6,045 | 5,645 | 731 |
| 2015–2016 | 9,971 | 5,474 | 5,459 | 5,087 | 222 |
| 2017–2018 | 9,254 | 5,265 | 5,231 | 4,865 | 112 |
| **Total** | **50,259** | **27,278** | **27,207** | **25,300** | **2,874** |

**Main result.** Cox model with age and sex, stratified by cycle: each 1 SD (8.4 years) of PhenoAge acceleration carries a hazard ratio of **1.54 (95% CI 1.50–1.57)** for all-cause death, i.e. **1.29 per 5 years** of acceleration. The proportional-hazards test for acceleration gives p = 0.55 (no evidence of violation).

| Model | n | Deaths | HR per SD (95% CI) |
| --- | --- | --- | --- |
| Main: age + sex, stratified by cycle | 25,300 | 2,874 | 1.54 (1.50–1.57) |
| Survey weights (robust SE; no strata/PSU) | 25,300 | 2,874 | 1.59 (1.53–1.66) |
| Not winsorised (raw acceleration, SD 9.3 y) | 25,300 | 2,874 | 1.28 (1.26–1.29) |
| Drop raw acceleration > 5 SD | 25,204 | 2,835 | 1.61 (1.56–1.66) |
| Drop deaths in first 2 years | 24,790 | 2,364 | 1.50 (1.46–1.55) |
| Age < 60 | 16,613 | 508 | 1.66 (1.57–1.75) |
| Age ≥ 60 | 8,687 | 2,366 | 1.52 (1.48–1.57) |
| Men | 12,255 | 1,669 | 1.50 (1.45–1.55) |
| Women | 13,045 | 1,205 | 1.61 (1.55–1.68) |
| Each cycle separately | 4,427–5,645 | 112–935 | 1.50 to 1.62 |

Full table: `results/week3_cox.csv`.

![Survival by acceleration quartile](figures/km_accel_quartiles.png)

The curves are unadjusted; use the hazard ratios above for the size of the effect.

Decisions worth knowing:
- **Winsorising.** Acceleration has a long right tail: 96 people lie more than 5 SD above the mean, mostly with dialysis-level creatinine or extreme white-cell counts. Untreated, these points flatten the per-SD estimate (1.28) and made the age < 60 model fail to converge. Main models therefore clip acceleration at the 0.5th and 99.5th percentiles (−14.9 to +41.8 years; 254 people affected). Both versions are reported.
- **Age** is capped at 80 in all cycles, because 2005–2006 top-codes at 85 and later cycles at 80.
- **Dropping early deaths** barely changes the estimate, so the result is not driven only by people who were already very ill at the exam.
- C-index in the pooled data (in-sample): age 0.843, PhenoAge 0.865.

## Repository layout

- `load_nhanes.py`: download and load raw files
- `sql/`: queries, numbered in run order
- `build_analysis.py`: builds the analysis table with unit checks
- `phenoage.py`: PhenoAge, age acceleration, C-index and Cox model
- `r/check_bioage.R`: cross-check against the BioAge R package
- `nhanes_io.py`: download and parsing helpers shared by the pooling script
- `pool_cycles.py`: pools five cycles, harmonises units (CRP mg/L vs mg/dL), computes PhenoAge
- `cox.py`: Cox models, sensitivity analyses, survival figure

## Limitations (so far)

- Survey weights are used only in one sensitivity model, and without the NHANES strata and primary sampling units, so its confidence interval is too narrow for population inference. Results describe the sample.
- Serum creatinine is not recalibrated across cycles; NHANES documentation recommends adjustments for some cycles. This could shift PhenoAge slightly between cycles; the per-cycle hazard ratios are similar (1.50–1.62).
- CRP in 2005–2010 is standard CRP, in 2015–2018 high-sensitivity CRP; both are converted to mg/dL but the assays differ at low concentrations.
- PhenoAge was trained on NHANES III, not these cycles, so this is an external check of a published score, not a new model. Step 4 compares it with models trained here.

## Reference

Levine ME et al. (2018). An epigenetic biomarker of aging for lifespan and healthspan. *Aging*, 10(4), 573–591.
