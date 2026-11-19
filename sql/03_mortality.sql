-- Crude mortality per 1000 person-years by age group (mortstat: 1 = assumed deceased)
SELECT CASE WHEN a.ridageyr < 40 THEN '20-39' WHEN a.ridageyr < 60 THEN '40-59' ELSE '60+' END AS age_grp,
       COUNT(*) AS n,
       SUM(m.mortstat) AS deaths,
       ROUND(SUM(m.permth_exm) / 12.0, 0) AS person_years,
       ROUND(1000.0 * SUM(m.mortstat) / (SUM(m.permth_exm) / 12.0), 1) AS rate_per_1000py
FROM adults_raw a JOIN mort m ON m.seqn = a.seqn
WHERE m.eligstat = 1 AND m.permth_exm IS NOT NULL
GROUP BY 1 ORDER BY 1;
