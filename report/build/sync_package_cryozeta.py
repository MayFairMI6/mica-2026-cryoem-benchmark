"""Bring the CryoZeta part of the reproducibility package in line with scoring/results (2026-10-13). Idempotent.
For every scored model of the result sets CryoZeta (package name CryoZeta_interim) and CryoZeta_partial
(CryoZeta_partial_models): copy the model file, check its SHA-256 against the one stored with the score, copy the run
records (also for cases without a model) and rewrite those rows of models/MODEL_MANIFEST.json. Called by rebuild_report.sh."""
import glob, hashlib, json, os, shutil, sys
R = "/Users/anony/Downloads/MICA_GAIVI_run"; O = f"{R}/benchmark/gaivi_bench/out"; P = f"{R}/reproducibility_package"
SETS = (("CryoZeta", "CryoZeta_interim", "cryozeta_interim"), ("CryoZeta_partial", "CryoZeta_partial_models", "cryozeta_partial_models"))
man = [r for r in json.load(open(f"{P}/models/MODEL_MANIFEST.json")) if r["result_set"] not in [s[1] for s in SETS]]; bad = []
for src, name, recdir in SETS:
    res = json.load(open(f"{R}/scoring/results/{src}.json")); os.makedirs(f"{P}/models/{name}", exist_ok=True)
    for d in sorted(glob.glob(f"{O}/{src}/*/")):                                   # run records, model or not
        c = os.path.basename(d.rstrip("/")); os.makedirs(f"{P}/run_records/{recdir}/{c}", exist_ok=True)
        for f in ("run_record.json", "stage_npz_paths.txt"):
            if os.path.exists(d + f): shutil.copy2(d + f, f"{P}/run_records/{recdir}/{c}/{f}")
    for c, v in res.items():
        m = v.get("metrics") or {}
        if m.get("tm_score") is None: continue
        f = f"{O}/{src}/{c}/{c}.cif"; sha = hashlib.sha256(open(f, "rb").read()).hexdigest(); ok = sha == v.get("model_sha256")
        if not ok: bad.append(f"{src} {c}: file on disk is not the scored model")
        shutil.copy2(f, f"{P}/models/{name}/{c}.cif")
        man.append({"result_set": name, "case": c, "file": f"models/{name}/{c}.cif", "sha256": sha, "sha256_matches_scored_model": ok,
                    "tm_score": m["tm_score"], "tm_score_search": m.get("tm_score_search"), "tm_score_mapframe": m.get("tm_score_mapframe")})
    extra = set(os.path.basename(x)[:-4] for x in glob.glob(f"{P}/models/{name}/*.cif")) - set(r["case"] for r in man if r["result_set"] == name)
    if extra: bad.append(f"{name}: model files in the package without a score: {sorted(extra)}")
json.dump(man, open(f"{P}/models/MODEL_MANIFEST.json", "w"), indent=1)
CF = f"{R}/benchmark/cryozeta/cluster_files"   # job scripts as run on the cluster and the five-line recording change to cycle_predict.py
for src_, dst in (("cz_run.sbatch", "cz_run.sbatch"), ("cz_partial.sbatch", "cz_partial.sbatch"), ("envs/cz/patch/cycle_predict.diff", "cycle_predict_recording_change.diff")):
    if os.path.exists(f"{CF}/{src_}"): shutil.copy2(f"{CF}/{src_}", f"{P}/job_scripts/cryozeta/{dst}")
print("package CryoZeta models:", {n: sum(r["result_set"] == n for r in man) for _, n, _ in SETS}, "| manifest rows:", len(man))
if bad: print("PACKAGE CHECK FAILED:", bad); sys.exit(1)
