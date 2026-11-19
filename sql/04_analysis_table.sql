-- Analysis table in PhenoAge units: albumin g/L, creatinine umol/L, glucose mmol/L, ln(CRP mg/dL)
CREATE OR REPLACE TABLE analysis AS
SELECT a.seqn, a.ridageyr AS age, a.riagendr AS sex,
       a.lbxsal * 10       AS albumin_gl,
       a.lbxscr * 88.4     AS creat_umol,
       a.lbxsgl / 18.016   AS glucose_mmol,
       a.lbxsapsi          AS alp,
       a.lbxwbcsi          AS wbc,
       a.lbxlypct          AS lymph,
       a.lbxmcvsi          AS mcv,
       a.lbxrdw            AS rdw,
       LN(a.lbxhscrp / 10) AS lncrp,
       m.mortstat, m.permth_exm
FROM adults_raw a JOIN mort m ON m.seqn = a.seqn
WHERE m.eligstat = 1 AND m.permth_exm IS NOT NULL
  AND a.lbxsal IS NOT NULL AND a.lbxscr IS NOT NULL AND a.lbxsgl IS NOT NULL AND a.lbxsapsi IS NOT NULL
  AND a.lbxwbcsi IS NOT NULL AND a.lbxlypct IS NOT NULL AND a.lbxmcvsi IS NOT NULL AND a.lbxrdw IS NOT NULL
  AND a.lbxhscrp > 0;
