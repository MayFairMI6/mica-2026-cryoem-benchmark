"""Apply the corrected-build US-align scores to results/*.json (2026-10-13; post hoc evaluator correction).
For every scored model, usalign_fixed_build_audit.json holds two evaluations made with US-align Version 20261012
(commit 746ad8a): the complex-mode search and the map-frame evaluation (-se). The primary TM-score is the higher of
the two (TM-score is a maximum over superpositions, and the map frame is a valid superposition for a model built
into the map). Aligned Calpha length and sequence identity are taken from the same evaluation. The value recorded at
first scoring with Version 20260920 is kept as tm_score_recorded. Idempotent."""
import json, os
A = {(r["set"], r["case"]): r for r in json.load(open("usalign_fixed_build_audit.json"))}
n = 0; miss = []
for s in sorted(set(k[0] for k in A)):
    p = f"results/{s}.json"
    if not os.path.exists(p): continue
    R = json.load(open(p))
    for c, x in R.items():
        m = x.get("metrics") or {}
        if m.get("tm_score") is None: continue
        a = A.get((s, c))
        if not a or "tm" not in a["search"] or "tm" not in a["mapframe"]: miss.append((s, c)); continue
        if "tm_score_recorded" not in m: m["tm_score_recorded"], m["aligned_ca_recorded"], m["seq_identity_recorded"] = m["tm_score"], m.get("aligned_ca"), m.get("seq_identity")
        for k in ("tm_score_plain", "tm_source"): m.pop(k, None)                      # fields of the superseded three-evaluation step
        m["tm_score_search"], m["tm_score_mapframe"] = a["search"]["tm"], a["mapframe"]["tm"]
        best = "search" if a["search"]["tm"] >= a["mapframe"]["tm"] else "mapframe"
        m["tm_score"], m["aligned_ca"], m["seq_identity"], m["tm_source"] = a[best]["tm"], a[best]["aligned"], a[best]["seq_id"], best
        m["usalign_build"] = "20261012 (746ad8a)"; n += 1
    json.dump(R, open(p, "w"), indent=1)
print("models updated:", n, "| missing from the audit:", miss)
