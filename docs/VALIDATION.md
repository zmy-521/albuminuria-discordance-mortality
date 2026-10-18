# Release validation (2026-10-18)

- 29 Python files parsed/compiled without importing analytical top-level code.
- 12 R files parsed with R 4.6.1; survey 4.5, survival 3.8-6 and splines imported. Initial machine locale/library-path issues were resolved in the private test launcher by using the existing runtime alias and explicit UTF-8 locale; no analysis formulas were changed.
- 2 Node workbook files passed `node --check`; the nonstandard spreadsheet library was not tested for public availability.
- Ten Python dependencies imported successfully.
- New setup utility ran in an isolated private test workspace without fitting models.
- Existing five-figure regression test: release renderers plus local frozen inputs produced exact approved final PNG pixels for all five figures, including final Figure 2/4 wording. No new model estimates. The patient plotting CSVs and test images remain private, outside the repository.
- `.gitignore` tested in a temporary local Git repository: 14 sensitive probes ignored and six allowed probes retained; no network, login, commit or push.
- Frozen source SHA256 and released-file manifest hashes verified. All original sources unchanged.
- Recursive text/XML/extension/size scan: no unresolved personal path, token, credential, restricted record file or unexpected large-file finding. Long decimal and URL false positives were reviewed and the scanner refined; no values were removed to pass the scan.

No analytical model training, Cox/logistic refitting, new cohort selection, bootstrap, new figure design, end-to-end public-data rerun or manuscript editing was performed. The existing figure regression is a rendering test with private inputs, **not proof that all data preparation/statistical stages were reproduced**. No participant-level models/data or private logs are included.

The final author must review the repository contents and decide licensing before uploading. `.gitignore` cannot prevent manual browser uploads; upload the audited repository ZIP only, not the parent release/audit directory.
