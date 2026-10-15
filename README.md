# Benchmarking MICA and recent cryo-EM protein model-building methods on 2026 EMDB/PDB entries

Code and documents for a benchmark of MICA (Gyawali, Dhakal & Cheng, *Communications Chemistry* 2025) against
ModelAngelo, EModelX(+AF), CryoAtom and EMProt, with CryoZeta as a later, separately reported addition, on 24 cryo-EM
maps whose EMDB and PDB entries were first released in 2026 (six per resolution band from 1.5 to 4.0 Å).

- **Report:** `report/` (text source `REPORT_DRAFT.md`, all result tables `RESULTS_TABLES.md`, figures).
- **Dataset and benchmark outputs** (input sequences, reference models, all scored models, run records, scores):
  Zenodo, DOI `<ZENODO_DOI>`.

## Contents

| Folder | What it holds |
|---|---|
| `dataset/` | The frozen 24-case candidate file (EMDB/PDB accessions, initial release dates, scope replacements and reasons), the case list and round order (`case_rounds.csv`, `selected12.csv`), the scripts that selected and ordered the cases (`select12.py`, `make_rounds.py`), map parameters and integrity reports, and the homology and density audits. |
| `alphafold_server_job_inputs/` | The AlphaFold Server job files (sequences only) used for MICA and the scripts that wrote them. |
| `job_scripts/` | The Slurm (GAIVI) and Kaggle scripts that installed and ran each method, including CryoZeta and its stage-recording change (`cycle_predict_recording_change.diff`). |
| `scoring/` | The evaluator (`evaluate.py`), the per-method scoring driver (`score_benchmark.py`), the table generator (`make_tables.py`), the US-align re-scoring scripts and their outputs, the synthetic control checks, and the scores (`results/*.json`). |
| `report/` | Report source, generated tables, per-case results (`raw_results_24cases.csv`), figures and the build and consistency-check scripts. |

## Reproducing the tables

`scoring/make_tables.py` reads only `scoring/results/*.json` and the files in `dataset/` and writes the 18 result tables.
The scripts retain their original execution logic; the absolute paths near the top of each script point to the original
project layout and need to be set to your clone.

To re-score models, download the dataset from Zenodo and run `scoring/score_benchmark.py METHOD canonical "<pattern>"`.
It needs US-align Version 20261012 (commit 746ad8a; see `scoring/USalign_build/BUILD_INFO.txt` — earlier builds have a
complex-mode bug that affects some scores), cctbx (`chain_comparison`) and gemmi.

## Not included

- **Cryo-EM maps:** public in EMDB under the accessions listed in `dataset/primary24_scope_screened_candidate.csv`;
  `job_scripts/kaggle/kaggle_bench12_inputs/canonicalize_map.py` reproduces the axis-order canonicalisation used.
- **The scripts of the initial 2026 EMDB query and metadata screen:** their output is the frozen candidate file in `dataset/`.
- **Phenix** (licensed; version 1.20.1-4487 was used), **AlphaFold Server outputs**, and the methods' code and weights
  (versions and commits are given in the report, section 2).

## Licence

Code: MIT (see `LICENSE`); the dataset on Zenodo is CC BY 4.0. The reference structures and maps are from the PDB and EMDB (CC0); the methods' weights carry their
own licences.
