"""Audit 3 (training/homology proxy) — RCSB sequence search (MMseqs2), identity >= 25%, e <= 0.1, per benchmark entity,
restricted to PDB entries released before 2023-05-01 (MICA's documented training cutoff) and before 2026-01-01.
Best hit = max(identity x query coverage). Usage: audit3_homology.py OUT.json FASTA [FASTA ...]"""
import json, sys, time, urllib.request
def search(seq, before):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "sequence", "parameters": {"evalue_cutoff": 0.1, "identity_cutoff": 0.25, "sequence_type": "protein", "value": seq}},
            {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_accession_info.initial_release_date", "operator": "less", "value": before}}]},
         "return_type": "polymer_entity", "request_options": {"results_verbosity": "verbose", "paginate": {"start": 0, "rows": 200}, "scoring_strategy": "sequence"}}
    req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query", data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                if r.status == 204: return {"best": None, "n": 0}
                d = json.loads(r.read() or b"{}")
            break
        except urllib.error.HTTPError as e:
            if e.code == 204: return {"best": None, "n": 0}
            time.sleep(5 * (attempt + 1))
        except Exception: time.sleep(5 * (attempt + 1))
    else: return {"error": "request failed"}
    best = None
    for h in d.get("result_set", []):
        for svc in h.get("services", []):
            for nd in svc.get("nodes", []):
                for mc in nd.get("match_context", []) or []:
                    ident = mc.get("sequence_identity"); cov = (mc["query_end"] - mc["query_beg"] + 1) / len(seq) if "query_end" in mc else mc.get("query_coverage")
                    if ident is None or cov is None: continue
                    if best is None or ident * cov > best["identity"] * best["query_coverage"]:
                        best = {"entry": h["identifier"], "identity": round(ident, 3), "query_coverage": round(cov, 3)}
    return {"best": best, "n": d.get("total_count", 0)}
out = {}
for fa in sys.argv[2:]:
    for b in open(fa).read().strip().split(">")[1:]:
        h, *s = b.split("\n"); eid = h.split("|")[0]; seq = "".join(s).strip()
        if len(seq) < 10: out[eid] = {"length": len(seq), "pre_2023-05": {"best": None, "n": 0}, "pre_2026": {"best": None, "n": 0}, "note": "too short for sequence search"}; continue
        out[eid] = {"length": len(seq), "pre_2023-05": search(seq, "2023-05-01"), "pre_2026": search(seq, "2026-01-01")}
        print(eid, len(seq), out[eid]["pre_2023-05"], flush=True); time.sleep(0.5)
json.dump(out, open(sys.argv[1], "w"), indent=1)
