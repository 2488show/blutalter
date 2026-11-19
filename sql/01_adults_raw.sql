-- Adults (>= 20 y) who completed the exam, with the 9 PhenoAge markers.
-- LEFT JOIN keeps people without lab results, so missingness stays visible.
CREATE OR REPLACE VIEW adults_raw AS
SELECT d.seqn, d.ridageyr, d.riagendr,
       c.lbxwbcsi, c.lbxlypct, c.lbxmcvsi, c.lbxrdw,
       b.lbxsal, b.lbxscr, b.lbxsgl, b.lbxsapsi,
       h.lbxhscrp
FROM demo d
LEFT JOIN cbc    c ON c.seqn = d.seqn
LEFT JOIN biopro b ON b.seqn = d.seqn
LEFT JOIN hscrp  h ON h.seqn = d.seqn
WHERE d.ridageyr >= 20 AND d.ridstatr = 2;
