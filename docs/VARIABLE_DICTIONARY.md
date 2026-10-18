# Core variables

The original code controls transformations. `config/nhanes_variable_index.csv` contains metadata lookup columns only.

|Field|Source/definition|Unit/type|
|---|---|---|
|SEQN|NHANES linkage key; no actual records bundled|integer|
|age / age10|RIDAGEYR /10|years / decades|
|sex|RIAGENDR|1 male, 2 female|
|race|RIDRETH1|original categories|
|cycle / year|survey cycle / start year|factor / integer|
|scr_mgdl|LBXSCR or 2001–2002 LBDSCR, original calibration|mg/dL|
|egfr / egfr10|CKD-EPI 2021 /10|mL/min/1.73 m²|
|ualb_mgl|URXUMA; µg/mL numerically equals mg/L|mg/L|
|ucr_mgdl|URXUCR before coordination|mg/dL|
|ucr_adjusted|CDC-coordinated urine creatinine|mg/dL|
|acr_adjusted|100*ualb_mgl/ucr_adjusted|mg/g|
|group / group20|30 or later 20 mg/L concentration threshold, ACR30|G0–G3|
|hba1c|LBXGH|%|
|bmi / bmi5|BMXBMI /5|kg/m²|
|smoking|SMQ020/SMQ040|never/former/current; unknown missing|
|htn|BPQ020|yes1/no0; otherwise missing|
|cvd|MCQ160B/C/D/E/F|any yes1; all no0; unresolved missing|
|pir|INDFMPIR|income-to-poverty ratio|
|weight|WTMEC4YR/WTMEC2YR coordination|survey weight|
|SDMVPSU / SDMVSTRA|official design variables|cluster / stratum|
|death / months|MORTSTAT / PERMTH_EXM|0/1 / months|
|time|max(months,0.5)|months|
|domain / ccB|study domain / complete Model B|indicator|
|flow|sum URXVOL1–3 / sum matching URDTIME1–3|mL/min|
|AER / CER|ualb_mgl*flow / ucr_adjusted*flow/100|µg/min / mg/min|

Hospital conversions: urine creatinine mmol/L*11.312 -> mg/dL; serum creatinine µmol/L/88.4 -> mg/dL; urine albumin mg/L. Same-specimen calculated UACR preferred. Qualified results are not arbitrarily imputed. Identity-related field names in hospital code are schema requirements, not identifier values. Instrument/reagent provenance cannot be inferred from the scripts.
