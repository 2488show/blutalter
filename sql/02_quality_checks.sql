-- Non-missing count per marker
SELECT COUNT(*) AS n_adults,
       COUNT(lbxsal) AS has_albumin, COUNT(lbxscr) AS has_creat, COUNT(lbxsgl) AS has_glucose,
       COUNT(lbxsapsi) AS has_alp, COUNT(lbxwbcsi) AS has_wbc, COUNT(lbxlypct) AS has_lymph,
       COUNT(lbxmcvsi) AS has_mcv, COUNT(lbxrdw) AS has_rdw, COUNT(lbxhscrp) AS has_crp
FROM adults_raw;

-- hs-CRP missingness and 99th percentile by age group
SELECT CASE WHEN ridageyr < 40 THEN '20-39' WHEN ridageyr < 60 THEN '40-59' ELSE '60+' END AS age_grp,
       ROUND(100.0 * SUM(CASE WHEN lbxhscrp IS NULL THEN 1 ELSE 0 END) / COUNT(*), 1) AS crp_missing_pct,
       ROUND(quantile_cont(lbxhscrp, 0.99), 2) AS crp_p99
FROM adults_raw
GROUP BY 1 ORDER BY 1;
