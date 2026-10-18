# Release adaptations

Original hashes identify authentic sources; public origin paths are logical, not machine paths. Full local paths are in an author-only manifest outside this repository. No originals were edited.

- `code/data_preparation/common.py`: Environment-based path substitution only.
- `code/data_preparation/03_nhanes.py`: Environment-based path substitution only.
- `code/hospital_classification/01_extract_wanbei.py`: Extract Wanbei raw LIS reader; remove unrelated NHANES index exploration; paths parameterized.
- `code/hospital_classification/common.py`: Environment-based path substitution only.
- `code/hospital_classification/04_hospitals.py`: Environment-based path substitution only.
- `code/hospital_classification/core.py`: Environment-based path substitution only.
- `code/hospital_classification/05_index_assessments.py`: Environment-based path substitution only.
- `code/cohort_construction/00_link_mortality.py`: Extract original save/eg/xpt/nhanes functions; omit hospital/CHARLS/other exploratory branches; preserve NHANES construction and counts.
- `code/cohort_construction/01_urine_harmonization.py`: Keep original NHANES urine-coordination prefix and original cohort save; omit SSCARD/CHARLS/hospital branches.
- `code/cohort_construction/02_prepare_stage1.py`: Environment-based path substitution only.
- `code/survival_analysis/stage1/02_survival.R`: Environment-based path substitution only.
- `code/survival_analysis/stage1/03_absolute_risk.R`: Environment-based path substitution only.
- `code/cohort_construction/03_prepare_stage2.py`: Remove guard requiring private Stage1 release manifest; retain all cohort checks and definitions; paths parameterized.
- `code/survival_analysis/stage2/02_cox.R`: Environment-based path substitution only.
- `code/survival_analysis/stage2/03_fixed_risk.R`: Environment-based path substitution only.
- `code/survival_analysis/stage2/04_selection.R`: Environment-based path substitution only.
- `code/survival_analysis/stage3/02_supplements.R`: Environment-based path substitution only.
- `code/urine_flow/01_prepare.py`: Remove private protocol copy/hash guard; preserve raw reads, matching, flow calculations and checks; paths parameterized.
- `code/component_analysis/02_analysis.R`: Environment-based path substitution only.
- `code/quality_control/stage1/05_result_QC.R`: Environment-based path substitution only.
- `code/quality_control/stage2/05_independent_QC.R`: Environment-based path substitution only.
- `code/quality_control/stage3/05_R_QC.R`: Environment-based path substitution only.
- `code/quality_control/stage1/04_independent_QC.py`: Remove legacy-v3 private hash manifest check only; retain raw-source scalar reconstruction and assertions.
- `code/quality_control/review/03_independent_QC.py`: Remove private frozen-input manifest guard only; retain source-hash checks and statistical QC.
- `code/quality_control/01_baseline_descriptives.R`: Environment-based path substitution only.
- `code/quality_control/build_supplement_tables.py`: Retain original table-payload assembly only; omit internal narrative and file-delivery operations.
- `code/quality_control/build_main_tables.py`: Extract clinical table-payload prefix; exclude retired ML/renal sheets and internal narrative; add same payload-write operation.
- `code/figures/01_inspect_objects.R`: Environment-based path substitution only.
- `code/figures/02_export.py`: Environment-based path substitution only; source-list bookkeeping excludes old script paths after relocation.
- `code/figures/source/base_plotting.py`: Parameterize figure root; mortality private handoff replaced by identical aggregate P input; original rendering functions retained.
- `code/figures/source/render_v2.py`: Parameterize figure root; mortality private handoff replaced by identical aggregate P input; original rendering functions retained.
- `code/figures/source/render_v3.py`: Parameterize figure root; mortality private handoff replaced by identical aggregate P input; original rendering functions retained.
- `code/figures/source/v3_mortality.py`: Parameterize figure root; mortality private handoff replaced by identical aggregate P input; original rendering functions retained.
- `code/figures/render_candidates.py`: Adapt original accepted presentation driver paths; omit original-image/private-profile checks unavailable publicly; preserve all plot data/axis and text assertions.
- `code/figures/render_text_revision.py`: Adapt original accepted presentation driver paths; omit original-image/private-profile checks unavailable publicly; preserve all plot data/axis and text assertions.
- `config/raw_reconstruction_original.json`: Unchanged historical config; final mortality builder excludes post-2016 releases then appends original 2017-2018.
- `code/hospital_classification/06_group_counts.py`: Extract original hospital count block only; no ML or private laboratory reference audit.
- `code/hospital_classification/07_classification_fractions.py`: Extract hospital fraction block only; omit entire ML preprocessing and overlap analysis.
- `config/nhanes_variable_index.csv`: Metadata only: years/file/variable/label in original order.
- `results/figure_inputs/02_supplement_statistical_results.xlsx`: Approved figure aggregate cells only; retain header row4 and source row coordinates.
- `results/figure_inputs/19_final_clinical_results.xlsx`: Approved figure aggregate cells only; retain header row4 and source row coordinates.
- `results/figure_inputs/07_Frozen_component_coefficients.csv`: Unchanged aggregate coefficients/covariance; no participant records.
- `results/figure_inputs/08_Frozen_component_covariance.csv`: Unchanged aggregate coefficients/covariance; no participant records.

New tools are packaging/setup/static-test utilities, not historical analysis code. Omitting private document/hash guards does not change statistical algorithms or establish preregistration.

- `code/quality_control/export_main_workbook.mjs`: Use ODKD_WORK_ROOT; original workbook renderer and dependency retained.

- `code/quality_control/export_supplement_workbook.mjs`: Use ODKD_WORK_ROOT; original workbook renderer and dependency retained.

- `code/quality_control/check_main_workbook.py`: Extract original cell comparison/P-zero checks; omit preview contact sheets; path parameterized.

- `code/quality_control/check_supplement_workbook.py`: Extract original cell comparison/P-zero checks; omit preview contact sheets; path parameterized.
