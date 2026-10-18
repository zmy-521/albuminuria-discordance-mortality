# Albuminuria discordance and mortality

Code for a study of discordance between urinary albumin concentration and urinary albumin-to-creatinine ratio (UACR), and subsequent all-cause mortality in adults with diabetes and creatinine-based eGFR >=60 mL/min/1.73 m虏.

## Release status
The original executed Python/R scripts have been traced and hashed. Release adaptations cover portable paths, extraction of relevant blocks, exclusion of private reporting/packaging, and the previously approved figure wording. No models were fitted for release preparation. This is not an independently rerun end-to-end replication. No license has been selected; the author must decide licensing before granting reuse rights.

The 5,905-person NHANES cohort and 5,155-person Model B sample (1,115 deaths) are unchanged. There is no ML/PGSS pipeline. Authentic raw-reconstruction helpers retain some ancillary variables and D1鈥揇5 calculations because later processing depended on this historical code; these are not additional analyses or claims for this clinical paper.

## Structure
- `code/data_preparation`: raw NHANES reader and equations.
- `code/cohort_construction`: mortality linkage, urine-creatinine coordination and covariates.
- `code/survival_analysis`: survey Cox, fixed-five-year risk, selection and later sensitivities.
- `code/component_analysis`: component Wald test, A2 comparison and survey summaries.
- `code/urine_flow`: composite urine collection, AER/CER and medication mapping.
- `code/hospital_classification`: authorized LIS processing and index selection.
- `code/figures`: original plot exporters and accepted presentation adapters.
- `code/quality_control`: original checks and table-payload assembly.
- `results/figure_inputs`: approved aggregate figure cells and model coefficients/covariance only.

Read [reproducibility instructions](docs/REPRODUCIBILITY.md), [data availability](docs/DATA_AVAILABILITY.md), [study design](docs/STUDY_DESIGN.md) and [variable dictionary](docs/VARIABLE_DICTIONARY.md). The [source manifest](docs/SOURCE_MANIFEST.csv) and [release changes](docs/RELEASE_CHANGES.md) distinguish copied sources from adaptations. Full local source paths stay outside this repository.

## Setup
Use Python 3.11+ and an R installation compatible with the recorded environment. Install `requirements.txt`; consult `config/R_ENVIRONMENT.md`. Set environment variables as described in `config/config.example.yaml`, then:

```sh
python tools/initialize_workspace.py
python tools/check_release.py
```

These commands do not fit models. Obtain the public inputs and follow the explicit analytical commands in `docs/REPRODUCIBILITY.md`. Do not import analytical scripts to test them: many execute processing at top level. Raw/derived individual-level records and model objects must remain private.

Hospital data are not bundled; hospital classifications and scatterplots cannot be reproduced from this repository alone. Public NHANES derivative tables/model objects are also not bundled. Full re-estimation and replicate-weight analyses were not rerun for packaging. The existing five-figure rendering regression was rerun using private local frozen inputs; all five final PNGs matched the approved outputs pixel-for-pixel. Those individual-level plotting inputs are not bundled. See `docs/VALIDATION.md` for tests actually performed.
