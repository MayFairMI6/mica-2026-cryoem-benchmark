"""Robustness audit of the complex TM-score (2026-10-13). For every scored model, the CA files that the evaluator gave
to US-align are re-used unchanged and US-align is run again in three ways:
  recorded : value stored at scoring time (command included '-m <matrix file>')
  plain    : same command without the '-m' option
  mapframe : same command with '-se' (no superposition search; the model is scored where it sits in the map frame)
Writes usalign_robustness_audit.json. Nothing in results/ is modified."""
import json, glob, os, subprocess, sys, concurrent.futures as cf
EXE = os.path.abspath("USalign/USalign"); SETS = sys.argv[1:] or ["MICA", "EModelX(+AF)", "ModelAngelo", "CryoAtom", "EMProt", "CryoZeta"]
def one(d):
    m, c = d.split("/")[1], d.split("/")[2]; u = d + "/usalign"; out = {"set": m, "case": c}
    try:
        rec = json.load(open(d + "/metrics.json"))["paper_named_usalign_metrics"]
        out["recorded"] = {"tm": rec.get("tm_score_reference_normalized"), "aligned": rec.get("aligned_ca_length"), "seq_id": rec.get("sequence_identity")}
    except Exception as e: out["recorded"] = {"error": repr(e)[:80]}
    for name, extra in (("plain", []), ("mapframe", ["-se"])):
        try:
            o = subprocess.run([EXE, "prediction_CA.cif", "reference_CA.cif", "-mol", "prot", "-mm", "1", "-ter", "1"] + extra, cwd=u, capture_output=True, text=True, timeout=14400).stdout
            al = [l for l in o.split("\n") if l.startswith("Aligned length")][0]; tm = [l for l in o.split("\n") if l.startswith("TM-score") and "Structure_2" in l][0]
            out[name] = {"tm": float(tm.split()[1]), "aligned": int(al.split("=")[1].split(",")[0]), "seq_id": float(al.split("=")[-1])}
        except Exception as e: out[name] = {"error": repr(e)[:80]}
    return out
dirs = [d for s in SETS for d in sorted(glob.glob(f"eval_bench/{s}/EMD-*")) if os.path.exists(d + "/usalign/prediction_CA.cif")]
if os.environ.get("AUDIT_CASE"): dirs = [d for d in dirs if d.endswith(os.environ["AUDIT_CASE"])]
OUT = os.environ.get("AUDIT_OUT", "usalign_robustness_audit.json")
res = json.load(open(OUT)) if os.path.exists(OUT) and os.environ.get("AUDIT_RESUME") else []
done = {(r["set"], r["case"]) for r in res}; dirs = [d for d in dirs if (d.split("/")[1], d.split("/")[2]) not in done]
with cf.ThreadPoolExecutor(6) as ex:
    for r in ex.map(one, dirs):
        res.append(r); json.dump(res, open(os.environ.get("AUDIT_OUT", "usalign_robustness_audit.json"), "w"), indent=1); print(r["set"], r["case"], r["recorded"].get("tm"), r["plain"].get("tm"), r["mapframe"].get("tm"), flush=True)
