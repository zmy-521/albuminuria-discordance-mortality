# Frozen study design

- Ten nonoverlapping 1999–2018 NHANES cycles, official 2019 public mortality linkage. Adults with diabetes, eGFRcr >=60, paired urine measurements and eligible examination follow-up.
- Diabetes: self-reported diagnosis, insulin/oral diabetes medication or HbA1c >=6.5%; pregnancy exclusion follows the original code.
- CKD-EPI 2021 race-free eGFRcr, serum creatinine mg/dL. Original calibration: 1999–2000 `1.013*x+0.147`; 2005–2006 `0.978*x-0.016`; 2001–2002 LBDSCR is not additionally calibrated.
- Before 2007, urine creatinine x (mg/dL): `(1.02*sqrt(x)-0.36)^2` for x<75; `(1.05*sqrt(x)-0.74)^2` for 75<=x<250; otherwise `(1.01*sqrt(x)-0.10)^2`. UACR=`100*albumin(mg/L)/creatinine(mg/dL)`.
- Primary operational thresholds: concentration >=30 mg/L and UACR >=30 mg/g. G0 neither positive; G1 concentration only; G2 UACR only; G3 both. The concentration cutoff is not a universal guideline recommendation. The 20 mg/L analysis is a later literature-motivated sensitivity.
- Primary G2/G0; other contrasts retained. MEC time origin, all-cause death and right censoring. Original zero-month model handling remains 0.5 month.
- Weights: early four-year MEC weight*0.2; later two-year MEC weight*0.1. PSU SDMVPSU, strata SDMVSTRA, original full design before domain selection.
- Model A: group + natural spline age/10 (3 df) + sex + race/ethnicity + cycle. Model B adds BMI/5, HbA1c, eGFR/10, smoking, hypertension, previous CVD and PIR. Functions/missingness and complete-case flags are unchanged.
- Five-year source: 1999–2014, n=4315. Model B standardization target n=3789; known-outcome logistic fit n=3785. Early censoring is not observed survival. Survey JKn uncertainty and hypothetical censoring assignments retain their distinct interpretations.
- Component model excludes group; uses log2 albumin and log2 urine creatinine. The coefficient-sum Wald test uses full covariance. A2: G2 versus G3 at 30<=UACR<300 with log2 UACR adjustment. These later explanatory analyses are not initial preregistration or causal decomposition.
- Flow: 2009–2018, 3539 eligible and 3261 evaluable. Composite volume/time. Short-interval AER=albumin*flow; CER=creatinine*flow/100. Not 24-hour excretion, a diagnostic gold standard or proof of low muscle mass.
- Wanbei: first measurement-qualified same-code/same-day paired assessment, then eGFR eligibility. Classification replication only, not mortality external validation. No hospital records are distributed.
