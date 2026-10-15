"""Re-score every model with the corrected US-align build (Version 20261012, commit 746ad8a; see
USalign_20261012_746ad8a/BUILD_INFO.txt). For each scored model the CA files the evaluator wrote are re-used unchanged:
  search   : USalign pred ref -mol prot -mm 1 -ter 1
  mapframe : the same with -se (no superposition search; model scored where it sits in the map frame)
Writes usalign_fixed_build_audit.json (resumable). Nothing in results/ is modified."""
import json, glob, os, subprocess, sys, concurrent.futures as cf
EXE = os.path.abspath("USalign_20261012_746ad8a/USalign"); OUT = "usalign_fixed_build_audit.json"
SETS = sys.argv[1:] or ["MICA", "EModelX(+AF)", "ModelAngelo", "CryoAtom", "EMProt", "CryoZeta", "CryoZeta_partial", "MICA_unrefined", "EMProt_raw_single_chain", "EMProt_all_fragments_split"]
def one(d):
    out = {"set": d.split("/")[1], "case": d.split("/")[2]}
    for name, extra in (("search", []), ("mapframe", ["-se"])):
        try:
            o = subprocess.run([EXE, "prediction_CA.cif", "reference_CA.cif", "-mol", "prot", "-mm", "1", "-ter", "1"] + extra, cwd=d + "/usalign", capture_output=True, text=True, timeout=14400).stdout
            al = [l for l in o.split("\n") if l.startswith("Aligned length")][0]; tm = [l for l in o.split("\n") if l.startswith("TM-score") and "Structure_2" in l][0]
            out[name] = {"tm": float(tm.split()[1]), "aligned": int(al.split("=")[1].split(",")[0]), "seq_id": float(al.split("=")[-1])}
        except Exception as e: out[name] = {"error": repr(e)[:80]}
    return out
res = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(r["set"], r["case"]) for r in res if "tm" in r.get("search", {}) and "tm" in r.get("mapframe", {})}; res = [r for r in res if (r["set"], r["case"]) in done]
dirs = [d for s in SETS for d in sorted(glob.glob(f"eval_bench/{s}/EMD-*")) if os.path.exists(d + "/usalign/prediction_CA.cif") and (s, d.split("/")[2]) not in done]
dirs.sort(key=lambda d: -os.path.getsize(d + "/usalign/reference_CA.cif"))       # largest first
print("to do:", len(dirs), "already done:", len(done), flush=True)
with cf.ThreadPoolExecutor(7) as ex:
    for f in cf.as_completed([ex.submit(one, d) for d in dirs]):
        r = f.result(); res.append(r); json.dump(res, open(OUT, "w"), indent=1); print(r["set"], r["case"], r["search"].get("tm"), r["mapframe"].get("tm"), flush=True)
print("FIXED_AUDIT_DONE", flush=True)
