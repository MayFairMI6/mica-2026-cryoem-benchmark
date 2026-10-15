"""Deterministic 12-case selection from the frozen 24-case candidate set.

Rule (fixed 2026-10-06, before any benchmark method has seen any of these cases):
within each of the 4 resolution strata, take the 3 cases with the lowest
`current_within_stratum_hash_rank` (the project's seeded ordering, seed 20261004).
No size, difficulty, reference-quality or method-output information is used.
Replacement rule: a selected case that fails a pending audit (full-map integrity,
density/reference scope, training/homology) is replaced by the next-ranked case from
the same stratum of the 24; the reason is recorded.
"""
import csv, hashlib, sys
src = sys.argv[1]; out = sys.argv[2]
rows = list(csv.DictReader(open(src)))
sel = []
for stratum in sorted({r["resolution_stratum"] for r in rows}):
    s = sorted((r for r in rows if r["resolution_stratum"] == stratum),
               key=lambda r: int(r["current_within_stratum_hash_rank"]))
    sel += s[:3]
keep = ["resolution_stratum", "current_within_stratum_hash_rank", "case_id", "emdb_id", "pdb_id",
        "resolution_A", "observed_CA_count", "deposited_sequence_residue_count", "citation_doi",
        "QA_flags", "full_map_integrity_status", "detailed_scope_density_review_status",
        "training_homology_audit_status", "coordinate_gzip_sha256", "fasta_sha256"]
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keep); w.writeheader()
    for r in sel: w.writerow({k: r[k] for k in keep})
print("source sha256:", hashlib.sha256(open(src, "rb").read()).hexdigest())
print("output sha256:", hashlib.sha256(open(out, "rb").read()).hexdigest())
