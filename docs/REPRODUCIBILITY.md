# Reproducibility and execution order

## Important scope
The commands below document the real pipeline and the supplied release adaptations. They were **not run to refit models during packaging**. Use a new writable workspace; do not point it at frozen research output. No participant-level table or RDS object is included. Original output-file names are preserved so stages can interoperate.

The early raw reader used an eleven-cycle historical configuration and variable index. It includes 2017–2020/2021–2023 as historical source coverage; the mortality builder explicitly keeps only through 2015–2016, then appends original 2017–2018. These overlapping releases are never added as mortality follow-up cohorts. This historical dependency is documented rather than silently rewritten. Extra raw dependencies and ancillary helper columns do not define additional analyses in this paper.

## Environment and raw layout
Set absolute local paths of your own, for example in PowerShell:
```powershell
$env:ODKD_WORK_ROOT = (Join-Path $PWD 'work')
$env:NHANES_ROOT = (Join-Path $PWD 'data/nhanes')
$env:WANBEI_RAW_ROOT = (Join-Path $PWD 'data/wanbei_authorized')
$env:ODKD_FIGURE_ROOT = (Join-Path $PWD 'work/figures')
python -m pip install -r requirements.txt
python tools/initialize_workspace.py
python tools/check_release.py
```
R requires survey, survival and splines with the versions noted in `config/R_ENVIRONMENT.md`. Use your own library; author-specific runtime paths were removed. No renv lock is fabricated.

Obtain official data from the CDC links in DATA_AVAILABILITY. Under NHANES_ROOT use:
```
index/variables.csv                 # supplied metadata copied by setup
archive/xpt/<cycle>/<file>.xpt
archive/codebooks/<cycle>/<file>.htm
archive/mortality/NHANES_YYYY_YYYY_MORT_2019_PUBLIC.dat
```
The raw reader tries `.XPT` then `.xpt`; later scripts use lowercase `.xpt`. Use those filenames on case-sensitive systems. `config/nhanes_variable_index.csv` lists source stems/variables, not downloadable patient records. Required components include DEMO, DIQ, glycated hemoglobin, biochemistry, urine albumin/creatinine, BMI, smoking, blood pressure and medical conditions, plus historical raw-reader ancillary components. Early aliases LAB10/L10_B/L10_C, LAB18/L40_B/L40_C, LAB16/L16_B/L16_C are explicitly in source. The historical reader also needs SSCYST_A/B. Flow additionally needs UCFLOW_F–J, MCQ_F–J, RXQ_RX_F–J, codebooks and `1988-2020/RXQ_DRUG.xpt`.

Downloads are manual in this release: no authentic complete raw-download script was found for this pipeline, and a new downloader is not represented as original code. The metadata index supplies lookup information; refer to CDC file pages to obtain the actual versions, and retain hashes locally.

## Public NHANES analytical route
After obtaining all original raw dependencies:
```sh
python code/data_preparation/03_nhanes.py
python code/cohort_construction/00_link_mortality.py
python code/cohort_construction/01_urine_harmonization.py
python code/cohort_construction/02_prepare_stage1.py
Rscript code/survival_analysis/stage1/02_survival.R
Rscript code/survival_analysis/stage1/03_absolute_risk.R
Rscript code/quality_control/stage1/05_result_QC.R
python code/quality_control/stage1/04_independent_QC.py
python code/cohort_construction/03_prepare_stage2.py
Rscript code/survival_analysis/stage2/02_cox.R
Rscript code/survival_analysis/stage2/03_fixed_risk.R
Rscript code/survival_analysis/stage2/04_selection.R
Rscript code/quality_control/stage2/05_independent_QC.R
Rscript code/survival_analysis/stage3/02_supplements.R
Rscript code/quality_control/stage3/05_R_QC.R
python code/urine_flow/01_prepare.py
Rscript code/component_analysis/02_analysis.R
python code/quality_control/review/03_independent_QC.py
Rscript code/quality_control/01_baseline_descriptives.R
```
The original QC R scripts may themselves fit verification models; they are provided for a future full replication, but **were not executed during release preparation**. JKn risk calculations are not cheap smoke tests. Component analysis contains the original limited eGFR spline sensitivity and reparameterization check; no new version was introduced.

Expected stage checkpoints: Stage1 cohort 5905, complete B 5155/1115 deaths; five-year source 4315, target 3789, fit 3785; flow source 3539, evaluable 3261; A2 complete-B 1120 (G2 222/G3 898). These are expected original outputs, not synthetic fixtures or evidence that the public release has been rerun successfully.

## Authorized hospital route
Only with institutional permission:
```sh
python code/hospital_classification/01_extract_wanbei.py
python code/hospital_classification/04_hospitals.py Wanbei
python code/hospital_classification/05_index_assessments.py Wanbei
```
This requires the original LIS schema and stable identifiers in protected storage. Salt and linkage files are generated locally and must not be committed. The original helper can also process Anyi, and the two aggregate count scripts loop over both centers; Anyi upstream extraction is deliberately outside this clinical release. Do not run those two-center scripts without both authorized indexed inputs. This is an explicit limitation, not permission to create placeholder patients. Wanbei expected urine counts are 6451 ->2414 ->2261 under the original index rule.

## Tables
`build_main_tables.py` and `build_supplement_tables.py` preserve real clinical table-payload assembly. The main builder still needs the original hospital aggregate audit files; the retired ML/renal sheets and internal narrative were removed. No new statistical calculations replace missing sources. The original workbook exporter is documented separately if its nonstandard spreadsheet dependency is available. CSV audit outputs remain the numerical source even when workbook styling cannot be reproduced.

## Figure inputs and commands
The original renderers use aggregate workbooks under ODKD_FIGURE_ROOT/Inputs. Setup installs the reviewed **reduced** workbooks and coefficients. Full anonymous plotting rows are not supplied. The actual export pipeline is `01_inspect_objects.R` followed by `02_export.py` and requires the original model objects, NHANES derivatives and authorized Wanbei data. It writes private linkage files under work, never into the public repository. Its provenance list has been adapted to omit relocated historical script paths.

|Figure|Actual function|Required inputs|Final export|
|---|---|---|---|
|1|render_v3.classification_v3|01_NHANES_5905_plot.csv, 02_Wanbei_2261_plot.csv, aggregate classifications|PNG/TIFF/PDF|
|2|v3_mortality.make_mortality|Group counts, Cox evidence, Period Cox, period_interaction.json|PNG/TIFF/PDF; final text revision|
|3|render_v3.components_v3|NHANES actual Model B/A2 flags; component coefficients/covariance; Effects|PNG/TIFF/PDF|
|4|render_v2.absolute_risk_v2|Absolute risk aggregate cells|PNG/TIFF/PDF; final text revision|
|5|render_v3.flow_v3|03_Flow_plot.csv; Flow characteristics weighted quartiles|PNG/TIFF/PDF|

After authorized local exports, copy the original three plotting CSVs into ODKD_FIGURE_ROOT/Inputs. The original render_v3 imports all three at module load, so the combined five-figure driver requires all three even if selecting only some functions. This coupling is retained, not hidden:
```sh
Rscript code/figures/01_inspect_objects.R
python code/figures/02_export.py
# Copy the locally generated three plotting CSVs from work/Figure_v3_PlotData_export/package to work/figures/Inputs.
python code/figures/render_candidates.py
python code/figures/render_text_revision.py
```
Use first-pass Figures 1/3/5 and text-revised Figures 2/4. Drivers retain accepted 170 mm layout, 600 dpi PNG, 300 dpi LZW TIFF and embedded-font PDF. Old-source-image equality checks are unavailable to public users and explicitly disabled; data/axis invariance assertions remain. No model fitting occurs in figure rendering. Supplementary `base_plotting.sensitivity` and `render_v2.a2_supp` retain real original plotting code but need their additional frozen aggregate worksheet rows; they were not rerun for this release.

## Unresolved reproducibility limits
No independent end-to-end re-estimation has been performed. No original automated full download script or renv lock was found. Hospital records/identities, fitted objects and internal audit logs are withheld. Table styling may depend on a nonstandard author-environment spreadsheet library. Historical raw-reader dependencies are broader than the mortality cohort. These limitations must accompany the public release rather than be described as full replication.


### Original workbook rendering dependency
`code/quality_control/export_main_workbook.mjs` and `export_supplement_workbook.mjs` are actual author-executed Node exporters. They require `@oai/artifact-tool`, which was available in the original author environment but is not vendored or promised to be generally installable here. Only `node --check` was performed during packaging; the dependency was not redistributed. The clinical payload builders exclude retired ML sheets, so their output is a clinical subset, not a byte-identical replacement for the historic full workbooks. The original numerical cell-check scripts are included. No missing exporter has been replaced by fabricated historical code.
