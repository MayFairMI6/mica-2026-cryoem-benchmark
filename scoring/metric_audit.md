# MICA metric audit and executable evaluation convention

This document distinguishes the six metrics named in the supplied MICA paper from the additional checks implemented for the 2026 benchmark. All numeric values in `control_results.csv` come from **synthetic evaluator controls** derived from experimental references. They measure evaluator behavior; the method field is empty because no MICA or comparator prediction is represented in those rows.

## The six named metrics

Let `Npred` be the number of retained protein Cα atoms in a prediction, `Nref` the corresponding number in the fixed experimental reference, `M` the number matched by the declared Phenix/cctbx procedure, and `I` the identical amino-acid pairs within that matching. US-align has its own structural alignment and therefore its own aligned count and identical count.

| Metric named by MICA | Definition and tool stated by the paper | Explicit convention implemented here |
|---|---|---|
| TM-score | US-align; normalized by true-structure length; 0–1 | Prediction is Structure_1; full retained experimental reference is Structure_2. Extract the Structure_2-normalized score. Use `-mol prot -mm 1 -ter 1` for the complete protein complex and one global rigid transform. |
| Cα match | Percentage of predicted Cα atoms within a distance threshold of true Cα atoms; Phenix `chain_comparison` | `100 × M/Npred` from official cctbx counts, using explicit 3.0 Å inclusive distance cutoff, all target and query chains, no crystallographic symmetry, and no unique-chain reduction. |
| Cα quality score | Combines Cα match precision and modeled/true residue count | The cited precursor Cryo2Struct gives the explicit formula: `(Cα match/100) × Npred/Nref = M/Nref`, with consistent atom/residue counts. Save a fraction; multiply by 100 for a percentage presentation. |
| Aligned Cα length | Number in the US-align structural alignment | Parse US-align's `Aligned length`. Preserve its raw output; this count is distinct from the 3 Å cctbx count. |
| Sequence identity | US-align identical aligned residues divided by structural alignment length | Parse the separate US-align `Seq_ID=n_identical/n_aligned` fraction. |
| Sequence match | Percentage of identical residues among the Phenix-aligned Cα atoms | Read official cctbx `get_match_percent('close') = 100 × I/M`. Keep the Phenix matching separate from US-align alignment. |

The MICA paper's main comparisons are **ModelAngelo** and **EModelX(+AF)**. Parenthetical `+AF` is part of the comparator specification, rather than a synonymous label for the sequence-only EModelX variant.

## Reproduction ambiguities resolved transparently

The supplied paper's performance-evaluation section is on PDF page 13. It does not give a numeric Cα distance threshold, a complete US-align command, the exact Phenix release, or the target/query uniqueness options. The inspected public MICA repository at commit `726d1bb14be04887a43b7f5aefcff0d843761d01` contains modeling and training scripts, but no evaluation script or US-align command was found. Consequently the convention above implements the published definitions explicitly; it is not a claim to have recovered the authors' exact original command.

**Three angstroms is supported by the predecessor and current implementation.** The cited Cryo2Struct paper describes its Phenix Cα match threshold as 3 Å, and the official `chain_comparison` documentation and source use `max_dist = 3.`. The source excludes a candidate only when its distance exceeds `max_dist`, so the implemented boundary is inclusive. A benchmark should publish this choice, preserve the command and version, and avoid treating it as a threshold explicitly written in the MICA paper.

**The Phenix column label is a material denominator trap.** In the official module, `get_close_to_target_percent()` returns `100 × matched/target_length`. The human-readable `% Found` column therefore has a reference denominator for a reference-first invocation. It cannot be copied into MICA's stated prediction-denominator Cα match field when prediction and reference lengths differ. The module's separate `CA SCORE` is also a different quantity: reference coverage divided by matched-residue RMSD. It is neither MICA's Cα match nor its Cα quality. The adapter extracts the integer counts and calculates the paper's prose definition, while also retaining official `% Found` as a separate diagnostic.

The **half-length 9QOP control** makes the issue observable: 541 of 1083 reference Cα atoms remain at their exact coordinates. Its paper-prose Cα match is 100%, its quality is 541/1083 ≈ 0.499538, and the official `% Found` uses the latter reference denominator. This comparison is evaluator validation, not evidence about any modeling method.

**TM-score is geometric.** The supplied MICA prose describes TM-score as considering sequence match. US-align's default structural alignment is explicitly sequence-independent. A TM-score sums distance-based contributions after alignment; amino-acid identity is a separate output. The control that changes every 9QOP amino-acid identity, while keeping all coordinates unchanged, returns reference-normalized TM-score 1.00000 and US-align sequence identity 0.000. The report should state the metric's operational definition rather than reproduce the ambiguous prose claim.

**US-align aligned length is not a 3 Å true-positive count.** The pinned source uses its own internal alignment and distance-selection logic, including a length-dependent search/final cutoff in the TM-align routines. Its displayed colons describe pairs below 5 Å; dots can still be aligned pairs. The cctbx matched count, US-align aligned count, and added spatial counts are retained separately.

**Homomer handling affects denominators.** The Phenix defaults can evaluate both unique and complete target chains, then select the choice with a higher matching fraction. Explicitly disabling this target reduction preserves the fixed deposited assembly and its copy counts. The official module avoids reusing target atoms through an order-dependent greedy procedure. The supplementary matcher instead maximizes matched cardinality and then minimizes total distance. These algorithms may differ on crowded or ambiguous coordinates; their outputs are not silently interchanged.

## Additional safeguards implemented in `evaluate.py`

The original-frame metrics assess coordinates where the method placed them in the density map. They apply no centering, scaling, periodic wrapping, crystal symmetry expansion, independent chain fit, or superposition. The separately named superposed metrics use the single global US-align transform. This separation prevents a high structural similarity score from concealing a displaced model.

At 1, 2, 3 and 5 Å, a sparse bipartite assignment first maximizes the number of distinct matched prediction/reference Cα pairs and then minimizes their summed distance. Precision is `matched/Npred`, reference recall is `matched/Nref`, and F1 is their harmonic mean. The added sequence-correct reference coverage is `identical matched/Nref`. Unknown amino-acid codes count as non-identical in these supplementary fields. Sequence is deliberately excluded from the geometric assignment objective.

The minimum-cost sparse assignment includes one dummy column per prediction. A dummy cost greater than all possible real-edge cost changes makes cardinality the primary objective. This avoids a dense `Npred × Nref` distance allocation. Duplicate predictions cannot increase the matched count above `Nref`. For the doubled 9QOP control, spatial precision is 0.5, recall is 1.0 and F1 is 2/3, whereas reference-normalized TM-score remains 1.0. Additional precision and excess-atom reporting is consequently necessary even when reference coverage is excellent.

Gemmi parses mmCIF and PDB. The retained subset is the first coordinate model's occupied amino-acid Cα atoms. Author and label chain IDs are jointly preserved, as are insertion codes. The highest-occupancy alternate atom is selected, with blank/A/alphabetical tie-breaking; non-finite or zero-occupancy coordinates are excluded. Modified amino acids recognized by Gemmi map to their parent one-letter codes. Distinct residue identities are checked explicitly. A canonical Cα-only mmCIF is sent to both structural tools so that their denominator counts can be reconciled. No missing coordinates are invented.

## Reference and map-frame findings from the real pilot data

| Pair | Deposited ASU observed protein Cα | ASU sequence residues, including unresolved residues | Evaluation assembly Cα | Correct frame: fraction of observed ASU Cα inside grid | Fraction if integer map starts are ignored |
|---|---:|---:|---:|---:|---:|
| EMD-53266 / 9QOP | 1083 | 1447 | 1083 | 1.000000 | 1.000000 |
| EMD-64908 / 9VAX | 1152 | 1484 | 1152 | 1.000000 | 1.000000 |
| EMD-65678 / 9W66 | 581 | 914 | 1162 | 1.000000 | 0.000000 |
| EMD-77112 / 13KK | 462 | 1428 | 462 | 1.000000 | 0.391775 |

These counts and fractions were computed from the downloaded coordinate files and MRC headers; raw data hashes and derived results are in `pilot_geometry_audit.json`. The map containment test uses the first and last voxel centers and is a geometric consistency check. Containment alone does not establish density support or modeling accuracy.

9W66's mmCIF declares assembly 1 as a dimer and applies operators `1,2`. Operator 2 maps `(x,y,z)` to `(-x+272.24002, -y+272.24002, z)`. EMDB processing metadata declares applied C2 symmetry. Exact deposited-operator expansion creates 1162 Cα atoms in two protein chains, all within the native map grid. The choice comes from deposition metadata, not from selecting whichever assembly gives a larger benchmark score. The generated mmCIF and full operator provenance are in `derived_references/`. The other three pilot entries declare identity-only assembly operations.

13KK contains 14 chains with 33 observed residues per chain, while each full deposited sequence has 102 residues. Its observed-coordinate fraction is 32.35%; 9W66's is 63.57%, 9QOP's 74.84%, and 9VAX's 77.63%. Full sequence length remains the sequence input to a method, with copy counts prescribed by the assembly. Reference-normalized structural metrics use the observed experimental Cα subset. A plausible model completion outside that observed subset cannot be adjudicated as a false physical prediction from coordinates alone. All-output precision must therefore be described as agreement with **deposited observed coordinates**; report unresolved sequence regions and map support separately. An evaluable-residue subset based on sequence/assembly mapping can be added, but should be declared before comparing methods and never selected by nearest-reference proximity to improve a score.

The two cropped maps have nonzero integer starts and zero floating origin. 13KK has starts `(15,15,24)` and approximately 1.65 Å voxels, giving first-voxel Cartesian position `(24.75,24.75,39.6)` Å. 9W66 has starts `(51,51,50)` and approximately 1.51244 Å voxels. The exact mapping must preserve voxel sizes, axis order, starts and physical origin. Nonorthogonal cells and simultaneous nonzero origins/starts are explicitly flagged for additional convention resolution.

A separately executed native MICA normalization stage on 13KK resampled from approximately 1.65 Å to 1 Å while retaining integer starts. Under the standard start-times-voxel coordinate mapping, that changes the physical first-voxel position to `(15,15,24)` Å. This observed preprocessing hazard should be resolved using reference-independent input-header canonicalization and exact transformation provenance before final coordinate metrics are accepted. Its current output is a partial neural stage, so no final model score is inferred from it.

## Actual executed controls and limits

Four references each have identity, large-translation, and half-per-chain truncation controls. 9QOP additionally has chain-order/label permutation, all-residue-type replacement, duplicated predictions, and removal of its largest chain: **16 synthetic controls total**. The controls test frame sensitivity, reference normalization, chain permutation tolerance, amino-acid identity, and overprediction behavior. US-align stdout/stderr, its rotation matrix, canonical input files, official cctbx logs, numeric results and expectation checks are retained per control. The Phenix/cctbx module is run in original map coordinates; aligned spatial metrics remain a separate supplement.

All 16 controls passed all 113 recorded expectations after official cctbx metrics were added to the previously executed US-align and supplementary results. The evaluator has therefore executed all six named metric definitions, with explicit modern options, on these controls. The full Phenix suite was unnecessary for chain comparison because the free official cctbx module provides that operation.

All model-level benchmark performance fields remain pending actual method-generated structures. These controls cannot rank MICA, ModelAngelo, EModelX(+AF), or newer systems. They establish that the measurement pipeline can distinguish several failure modes before expensive inference is committed.

## Sources and pinned implementations

1. Supplied `MICA_cryo-EM.pdf`, especially performance evaluation on p. 13; public publication: https://www.nature.com/articles/s42004-025-01718-5 . The attached PDF is the authoritative task source.
2. Cryo2Struct publication, explicit Cα quality formula and 3 Å Phenix convention: https://www.nature.com/articles/s41467-024-49647-6 .
3. Official Phenix chain-comparison documentation: https://phenix-online.org/documentation/reference/chain_comparison.html .
4. Official cctbx implementation: https://github.com/cctbx/cctbx_project/blob/b0747a4a55f29db3abe04358480d5867e94cb792/mmtbx/validation/chain_comparison.py . The installed `cctbx-base==2025.11` file has exactly the same SHA-256, `9a2dd855f7f5615f274f4e8265180580037de658803ca7bf8a90343fd6a45d55`, as the independently downloaded source. This identifies the executed implementation; it does not identify MICA's historical Phenix release.
5. Official US-align source: https://github.com/pylelab/USalign/tree/1fa25a958fe3095900eba2ef9b48562bdc462251 . Compiled release banner: `Version 20260920`. Flags and compiler details are saved in `runtime_manifest.json` and `usalign_help.txt`.
6. Official MICA repository: https://github.com/jianlin-cheng/MICA/tree/726d1bb14be04887a43b7f5aefcff0d843761d01 . Audited file tree and README are in `source_audit/`.
7. Recent metric context, a Cα-localization benchmark rather than a MICA whole-pipeline comparison: https://academic.oup.com/bioinformatics/article/42/6/btag350/8700507 .
