# Benchmarking MICA and recent cryo-EM protein model-building methods on 2026 EMDB/PDB entries

Code and documents for a benchmark of MICA (Gyawali, Dhakal & Cheng, *Communications Chemistry* 2025) against
ModelAngelo, EModelX(+AF), CryoAtom and EMProt, with CryoZeta as a later, separately reported addition, on 24 cryo-EM
maps whose EMDB and PDB entries were first released in 2026 (six per resolution band from 1.5 to 4.0 Å).

- **Report:** [report source](report/REPORT_DRAFT.md), [result tables](report/RESULTS_TABLES.md), and [figures](report/figures/).
- **Dataset and benchmark outputs** (input sequences, reference models, all scored models, run records, scores):
  Zenodo DOI [10.5281/zenodo.23581358](https://doi.org/10.5281/zenodo.23581358); [dataset record and download](https://zenodo.org/records/23581358).

## Contents

| Folder | What it holds |
|---|---|
| `dataset/` | The frozen 24-case candidate file (EMDB/PDB accessions, initial release dates, scope replacements and reasons), the case list and round order (`case_rounds.csv`, `selected12.csv`), the scripts that selected and ordered the cases (`select12.py`, `make_rounds.py`), map parameters and integrity reports, and the homology and density audits. |
| `alphafold_server_job_inputs/` | The AlphaFold Server job files (sequences only) used for MICA and the scripts that wrote them. |
| `job_scripts/` | The Slurm (GAIVI) and Kaggle scripts that installed and ran each method, including CryoZeta and its stage-recording change (`cycle_predict_recording_change.diff`). |
| `scoring/` | The evaluator (`evaluate.py`), the per-method scoring driver (`score_benchmark.py`), the table generator (`make_tables.py`), the US-align re-scoring scripts and their outputs, the synthetic control checks, and the scores (`results/*.json`). |
| `report/` | Report source, generated tables, per-case results (`raw_results_24cases.csv`), figures and the build and consistency-check scripts. |

## Reproducing the tables

Run `python3 scoring/verify_release.py` from the clone to check the frozen release dates, 120 primary scores and reported means using the Python standard library. It does not re-run inference or independently query database metadata.

`scoring/make_tables.py` also reads runtime records and audit files from the original project layout. Set `MICA_PROJECT_ROOT` to that full layout (containing `benchmark/`, `scoring/` and `report/`) to regenerate all 18 tables. A code-only clone is insufficient for that operation; the release check above works directly from the clone. Other archived scripts still contain paths that must be adapted to the run environment.

To re-score models, download the dataset from Zenodo and run `scoring/score_benchmark.py METHOD canonical "<pattern>"`.
It needs US-align Version 20261012 (commit 746ad8a; see `scoring/USalign_build/BUILD_INFO.txt` — earlier builds have a
complex-mode bug that affects some scores), cctbx (`chain_comparison`) and gemmi.

## Not included

- **Cryo-EM maps:** public in EMDB under the accessions listed in `dataset/primary24_scope_screened_candidate.csv`;
  `job_scripts/kaggle/kaggle_bench12_inputs/canonicalize_map.py` reproduces the axis-order canonicalisation used.
- **The scripts of the initial 2026 EMDB query and metadata screen:** their output is the frozen candidate file in `dataset/`.
- **Phenix** (licensed; version 1.20.1-4487 was used), **AlphaFold Server outputs**, and the methods' code and weights
  (versions and commits are given in the report, section 2).

## Use of AI tools

AI assistants were used during this project to help draft and debug scripts, and to help edit the report and
documentation. The benchmark design, case selection, method runs, scoring, and the interpretation of results were
directed, checked and decided by the author, who takes full responsibility for the code, data and conclusions.

## Licence

Code: MIT (see `LICENSE`); the dataset on Zenodo is CC BY 4.0. The reference structures and maps are from the PDB and EMDB (CC0); the methods' weights carry their
own licences.
