"""Apply the US-align robustness audit to results/*.json (2026-10-13).
TM-score is by definition the maximum over superpositions; US-align's complex-mode search does not always find it for
fragmented models (see usalign_robustness_audit.py). For every scored model three US-align evaluations of the same CA
files are available: the value recorded at scoring time, a repeat without the '-m' option, and the map-frame value
('-se', no superposition search). The primary TM-score is the highest of the three; aligned Calpha length and sequence
identity are taken from that same evaluation. All three values are kept in the result file. Idempotent."""
import json, os
A = {}
for f in ("usalign_robustness_audit.json", "usalign_robustness_audit_extra.json", "usalign_robustness_audit_9P9D.json"):   # later files override
    if os.path.exists(f):
        for r in json.load(open(f)): A[(r["set"], r["case"])] = r
n = ch = 0
for s in sorted(set(k[0] for k in A)):
    p = f"results/{s}.json"; R = json.load(open(p))
    for c, x in R.items():
        m = x.get("metrics") or {}; a = A.get((s, c))
        if not a or m.get("tm_score") is None: continue
        rec_tm = m.get("tm_score_recorded", m["tm_score"])
        if "tm_score_recorded" not in m: m["tm_score_recorded"], m["aligned_ca_recorded"], m["seq_identity_recorded"] = m["tm_score"], m.get("aligned_ca"), m.get("seq_identity")
        assert abs(a["recorded"]["tm"] - m["tm_score_recorded"]) < 1e-9, (s, c, a["recorded"]["tm"], m["tm_score_recorded"])
        cand = [("recorded", m["tm_score_recorded"], m["aligned_ca_recorded"], m["seq_identity_recorded"])]
        for k in ("plain", "mapframe"):
            if "tm" in a[k]: cand.append((k, a[k]["tm"], a[k]["aligned"], a[k]["seq_id"])); m[f"tm_score_{k}"] = a[k]["tm"]
        src, tm, al, sid = max(cand, key=lambda t: (t[1], t[0] == "recorded"))
        m["tm_score"], m["aligned_ca"], m["seq_identity"], m["tm_source"] = tm, al, sid, src
        n += 1; ch += tm - m["tm_score_recorded"] > 0.01
    json.dump(R, open(p, "w"), indent=1)
print("models updated:", n, "| primary TM higher than recorded by > 0.01 in:", ch)
