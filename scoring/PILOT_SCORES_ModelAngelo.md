# ModelAngelo pilot scores — 2026-10-06

**Scope warning.** These are the 4 *engineering pilots*, which the project's protocol keeps
separate from the 24-case accuracy cohort. One method only — **no comparison, no ranking**.
They test whether the prediction → frame correction → scoring pipeline works end to end.

## Evaluator
Project's own `metrics_work/evaluate.py` (copied unchanged) on this Mac (macOS arm64):
Python 3.12, NumPy 2.3.5, SciPy 1.17.0, Gemmi 0.7.5, cctbx-base 2025.11
(`chain_comparison.py` sha256 matches the recorded `9a2dd855…`), US-align @ `1fa25a95`.

**Validation:** 7 recorded synthetic controls (9QOP) rerun; **every computed metric matches**
the recorded values, including all 9 cctbx fields per control. Only labels and file paths differ.

**Deviation:** cctbx is installed `--no-deps` per the README, but needs `six` (installed
1.17.0). Without it the cctbx metrics silently fail; an earlier, looser check missed this.

## Frame correction (header-derived, no fitting)
ModelAngelo ignores the MRC ORIGIN and writes coordinates from the map's first voxel.
Each prediction was shifted by its cube map's recorded origin
(`canonical_cube.mrc.transform.json` → `output_origin_xyz_A_serialized`).
Provenance and hashes: `predictions_map_frame/*.provenance.json`.

## Results
| Pilot | Reference (observed Cα) | Pred Cα | TM-score | Aligned Cα | Seq. identity | Cα match (paper denom.) | Cα quality | Seq. match | Map-frame P / R @3 Å |
|---|---|---|---|---|---|---|---|---|---|
| EMD-65678 / 9W66 | assembly 1, 2 chains (1,162) | 1,144 | **0.982** | 1,142 | 0.987 | 99.7% | 0.982 | 98.8% | 0.997 / 0.982 |
| EMD-64908 / 9VAX | deposited (1,152) | 874 | **0.696** | 809 | 0.974 | 98.7% | 0.749 | 97.2% | 0.987 / 0.749 |
| EMD-53266 / 9QOP | deposited (1,083) | 97 | **0.031** | 34 | 0.941 | **100.0%** | 0.090 | 68.0% | 1.000 / 0.090 |
| EMD-77112 / 13KK | deposited (462) | 0 | null | null | null | null | null | null | null |

13KK is **null, not zero**: no model exists to score (protocol rule).

## Reading these numbers
- **9QOP is the denominator trap from the MICA paper, live:** 100% Cα match on a model
  covering 9% of the structure. Cα quality (0.090) and TM-score (0.031) expose it. Never
  quote Cα match alone.
- **9W66 is near-complete** once scored against the correct reference — the deposited
  dimer assembly. Against the single deposited chain it would look like over-building.
- Map-frame and superposed results agree closely, confirming the frame correction.
- **Pilot coverage depends heavily on the map/reference:** 13KK's deposited model itself
  covers only 32% of its sequence (462 of 1,428 residues).

## Corrections to earlier pilot notes
The "≈ count ratio" in `pilot/PILOT_PLAN.md` (Step 1b) used full FASTA length × copies as
the denominator and the single deposited 9W66 chain. Correct denominators are the observed
reference Cα counts above; 9W66's reference is the two-chain assembly.

---
# EModelX(+AF) pilot scores — 2026-10-07 (same validated evaluator)
Kaggle `mica-pilot-emodelx` v2 (v1 excluded: harness bug, output dir not created).
Frame: EModelX also ignores the MRC ORIGIN; header-derived shift by each canonical map's
origin (0 for 9QOP/9VAX; (77.1, 77.1, 75.6) Å for 9W66 → 98% of Cα within 3 Å after shift).
Final model = `<pdb>_EModelX(+AF)_all_atom_model.pdb` (PULCHRA all-atom).

| Pilot | Method | TM-score | Aligned Cα | Seq. id | Cα match | Cα quality | Seq. match | Map-frame P / R @3 Å | Runtime (T4) |
|---|---|---|---|---|---|---|---|---|---|
| 9W66 | ModelAngelo | 0.982 | 1,142 | 0.987 | 99.7% | 0.982 | 98.8% | 0.997 / 0.982 | 790.6 s |
| 9W66 | EModelX(+AF) | 0.499 | 581 | 0.988 | 96.7% | 0.497 | 98.3% | 0.970 / 0.499 | 402.5 s |
| 9VAX | ModelAngelo | 0.696 | 809 | 0.974 | 98.7% | 0.749 | 97.2% | 0.987 / 0.749 | 836.6 s |
| 9VAX | EModelX(+AF) | 0.872 | 1,042 | 0.907 | 89.4% | 0.924 | 81.0% | 0.904 / 0.935 | 1,033.5 s |
| 9QOP | ModelAngelo | 0.031 | 34 | 0.941 | 100.0% | 0.090 | 68.0% | 1.000 / 0.090 | 976.8 s |
| 9QOP | EModelX(+AF) | 0.618 | 741 | 0.602 | 91.8% | 0.650 | 51.3% | 0.937 / 0.664 | 690.2 s |
| 13KK | both | null | — | — | — | — | — | — | empty final models |

**Pilot-only stoichiometry effect (9W66).** The deposited FASTA lists one copy; the correct
reference is assembly 1 expanded by the deposited operators (2 chains). Methods that follow the
FASTA copy count (EModelX, and MICA via its 1-copy AF3 job) build one chain and are capped near
50% recall; ModelAngelo built both from density. **All 12 benchmark cases have
`oper_expression: "1"`**, so their deposited models already span assembly 1 and FASTA copy
counts match the references — the effect does not arise in the benchmark.

**13KK** is now empty for MICA (zero-AF ablation), ModelAngelo and EModelX — a case property
(only 32% of residues are modelled in the deposited structure), not evidence against any method.

---
# MICA pilot A — 2026-10-07
Kaggle `mica-pilot-mica-a` (MICA @ 726d1bb1, Phenix 1.20.1-4487, AF Server model_0 priors,
README Option A pipeline, T4 + 4 CPU cores). Frame: MICA also ignores the MRC ORIGIN;
header-derived shift → 99% of Cα within 3 Å of the 9W66 reference.

| Pilot | Outcome | TM | Aligned Cα | Seq. id | Cα match | Cα quality | Seq. match | Map-frame P / R @3 Å | Wall time |
|---|---|---|---|---|---|---|---|---|---|
| 13KK | **empty** final model (file present, 0 atoms) | null | — | — | — | — | — | — | 7,945.7 s |
| 9W66 | non-empty, 4,926 atoms (608 Cα, one chain) | 0.263 | 315 | 0.956 | 98.5% | 0.515 | 93.7% | 0.987 / 0.516 | 4,573.6 s |

**Runtime breakdown (9W66):** step 1 (AF3 processing) 0.2 min; **step 2 Phenix domain docking
66.3 min = 87% of 76.2 min** (paper: 85–97%); step 3 run.py 9.7 min (neural prediction 90 s on
GPU; sequence–structure alignment with AF3 462 s). Peak GPU 5.3 GB.

**Why TM 0.263 when map-frame recall ≈ EModelX's.** MICA's single chain is stitched across both
protomers of the dimer: of 603 Cα within 3 Å, 330 lie nearest copy A1 and 273 nearest copy A2,
with one switch along the chain (residue numbering matches the nearest reference residue for
554/603). EModelX's chain stays within one copy (584/584 in A1). The complex-mode TM-score
correctly penalises the cross-protomer chain. Likely contributor: the pilot-only 1-copy FASTA
(see above) — MICA was told one chain and assembled the best-fitting segments from both copies.
Consistent with the paper's own stated limitation (chain registration in complexes, p. 9).
**n = 1 pilot; not a conclusion.** Benchmark cases have correct copy counts.

13KK is now empty for all three methods including nominal MICA.

# MICA pilot B — 2026-10-07
Both non-empty: 9VAX 8,441 atoms (6,750 s), 9QOP 7,801 atoms (6,490 s); peak GPU ~5.3 GB.
Frame: origin 0 for both; 97% / 96% of Cα within 3 Å of the reference.

## All pilot scores (three core methods; pilots are NOT benchmark results)
| Pilot | Method | TM | Aligned Cα | Seq. id | Cα match | Cα quality | Seq. match | Map-frame P / R @3 Å |
|---|---|---|---|---|---|---|---|---|
| 9W66 | ModelAngelo | 0.982 | 1,142 | 0.987 | 99.7% | 0.982 | 98.8% | 0.997 / 0.982 |
| 9W66 | EModelX(+AF) | 0.499 | 581 | 0.988 | 96.7% | 0.497 | 98.3% | 0.970 / 0.499 |
| 9W66 | MICA | 0.263 | 315 | 0.956 | 98.5% | 0.515 | 93.7% | 0.987 / 0.516 |
| 9VAX | ModelAngelo | 0.696 | 809 | 0.974 | 98.7% | 0.749 | 97.2% | 0.987 / 0.749 |
| 9VAX | EModelX(+AF) | 0.872 | 1,042 | 0.907 | 89.4% | 0.924 | 81.0% | 0.904 / 0.935 |
| 9VAX | MICA | 0.861 | 1,008 | 0.938 | 94.4% | 0.885 | 90.5% | 0.946 / 0.888 |
| 9QOP | ModelAngelo | 0.031 | 34 | 0.941 | 100.0% | 0.090 | 68.0% | 1.000 / 0.090 |
| 9QOP | EModelX(+AF) | 0.618 | 741 | 0.602 | 91.8% | 0.650 | 51.3% | 0.937 / 0.664 |
| 9QOP | MICA | 0.797 | 901 | 0.806 | 90.4% | 0.820 | 73.7% | 0.930 / 0.843 |
| 13KK | all three | null | — | — | — | — | — | — |
