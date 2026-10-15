"""Assign all 24 frozen cases to inclusion rounds (declared 2026-10-06, before any run).

Within each resolution stratum, cases are ordered by `current_within_stratum_hash_rank`
(seed 20261004). Position k (1..6) in that order defines the round:
positions 1-3 -> core rounds 1-3 (the 12 selected cases), position 4 -> round 4,
position 5 -> round 5, position 6 -> round 6. Each round has one case per stratum
(rounds 1-3 are run together as the core).
"""
import csv, hashlib, sys
src, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(src)))
assigned = []
for stratum in sorted({r["resolution_stratum"] for r in rows}):
    s = sorted((r for r in rows if r["resolution_stratum"] == stratum),
               key=lambda r: int(r["current_within_stratum_hash_rank"]))
    for pos, r in enumerate(s, 1):
        assigned.append({"round": pos, "stage": "core" if pos <= 3 else "extension",
                         "resolution_stratum": stratum, "hash_rank": r["current_within_stratum_hash_rank"],
                         "case_id": r["case_id"], "resolution_A": r["resolution_A"],
                         "observed_CA_count": r["observed_CA_count"]})
assigned.sort(key=lambda a: (a["round"], a["resolution_stratum"]))
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(assigned[0])); w.writeheader(); w.writerows(assigned)
print("rounds sha256:", hashlib.sha256(open(out, "rb").read()).hexdigest())
