"""Score one method's benchmark predictions with the validated evaluator — identical procedure for every method.

usage: python score_benchmark.py METHOD MAP_KIND MODEL_GLOB
  METHOD     label (e.g. ModelAngelo)
  MAP_KIND   'canonical' or 'canonical_cube' — the map the method was given (sets the origin shift)
  MODEL_GLOB pattern with {case} and/or {pdb}/{emd} placeholders locating the FINAL model per case

Steps per case: (1) find the final model; missing/empty -> null (not zero); (2) shift by the input map's
recorded ORIGIN (header-derived; no fitting) and record the fraction of Cα within 3 Å of the reference
as a frame sanity check; (3) run evaluate.py against the deposited coordinates (assembly 1 = deposited
chains for every benchmark case); (4) append one row to results/<METHOD>.json.
"""
import glob, gzip, hashlib, json, os, shutil, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "vendor"))
import gemmi
BENCH = os.path.join(HERE, "..", "benchmark")
method, map_kind, pattern = sys.argv[1], sys.argv[2], sys.argv[3]
integ = json.load(open(f"{BENCH}/kaggle_bench12_prep/output/bench12/integrity_report.json"))["cases"]
for _f in sorted(glob.glob(f"{BENCH}/gaivi_bench/records/integrity_report_*.json")):          # extension rounds 4-6 (GAIVI)
    for _c, _v in json.load(open(_f))["cases"].items():
        if _c not in integ and _v.get("full_map_integrity") == "pass": integ[_c] = _v
os.makedirs(f"{HERE}/results", exist_ok=True); os.makedirs(f"{HERE}/bench_refs", exist_ok=True)
out_path = f"{HERE}/results/{method}.json"
res = json.load(open(out_path)) if os.path.exists(out_path) else {}
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
for case, rep in integ.items():
    emd, pdb = case.split("_")[0].replace("EMD-", ""), case.split("_")[1]
    hits = sorted(glob.glob(pattern.format(case=case, pdb=pdb, pdbl=pdb.lower(), emd=emd)))
    if not hits:
        continue                                   # not run (yet) for this case: no row
    src = hits[0]; row = {"case": case, "model": src}
    if case in res and res[case].get("model_sha256") == sha(src) and (res[case].get("metrics") or {}).get("tm_score") is not None:
        continue                                   # same model file already scored
    st = gemmi.read_structure(src)
    nat = sum(1 for m in st for c in m for r in c for a in r)
    if nat == 0:
        row.update(outcome="empty_model", metrics=None); res[case] = row; continue
    ref = f"{HERE}/bench_refs/{pdb}.cif"
    if not os.path.exists(ref):
        gz = (glob.glob(f"{BENCH}/kaggle_bench12_inputs/{case}/*.cif.gz") + glob.glob(f"{BENCH}/ext_inputs/{case}/*.cif.gz"))[0]
        with gzip.open(gz) as fi, open(ref, "wb") as fo: shutil.copyfileobj(fi, fo)
    o = rep[map_kind]["output_origin_xyz_A_serialized"]
    st.setup_entities()
    for m in st:
        for c in m:
            for r in c:
                for a in r: a.pos = gemmi.Position(a.pos.x + o[0], a.pos.y + o[1], a.pos.z + o[2])
    dst = f"{HERE}/predictions_map_frame/bench_{case}_{method}.cif"; st.make_mmcif_document().write_file(dst)
    P = np.array([a.pos.tolist() for c in st[0] for r in c for a in r if a.name == "CA"])
    R = np.array([a.pos.tolist() for c in gemmi.read_structure(ref)[0] for r in c for a in r if a.name == "CA"])
    within = float(np.mean([np.min(np.linalg.norm(R - p, axis=1)) <= 3 for p in P])) if len(P) else 0.0
    ev = f"{HERE}/eval_bench/{method}/{case}"; shutil.rmtree(ev, ignore_errors=True)
    p = subprocess.run([sys.executable, f"{HERE}/evaluate.py", "--prediction", dst, "--reference", ref, "--output", ev], capture_output=True, text=True,
                       env=dict(os.environ, MICA_EVAL_TIMEOUT_S=os.environ.get("MICA_EVAL_TIMEOUT_S", "14400")))
    row.update(outcome="non_empty_model", origin_shift_A=o, frame_check_frac_ca_within_3A=round(within, 3),
               model_sha256=sha(src), evaluator_exit=p.returncode)
    if p.returncode == 0:
        mt = json.load(open(f"{ev}/metrics.json")); u = mt["paper_named_usalign_metrics"]; ph = mt["phenix_paper_metrics"]
        s3 = [x for x in mt["supplementary_map_frame"] if x["threshold_angstrom_inclusive"] == 3][0]
        row["usalign_status"] = u.get("status")          # a US-align/cctbx time-out leaves that metric null, never zero
        row["metrics"] = {"tm_score": u.get("tm_score_reference_normalized"), "aligned_ca": u.get("aligned_ca_length"),
                          "seq_identity": u.get("sequence_identity"), "ca_match_pct": ph.get("ca_match_percent_paper_prose_denominator"),
                          "ca_quality": ph.get("ca_quality_fraction_paper_formula"), "seq_match_pct": ph.get("sequence_match_percent"),
                          "mapframe_precision_3A": s3["precision"], "mapframe_recall_3A": s3["recall"],
                          "pred_ca": mt["prediction"]["protein_ca_count"], "ref_ca": mt["reference"]["protein_ca_count"], "cctbx": ph.get("status")}
    res[case] = row; json.dump(res, open(out_path, "w"), indent=2)      # save after every case
    print(case, row["outcome"], row.get("frame_check_frac_ca_within_3A"), (row.get("metrics") or {}).get("tm_score"), flush=True)
json.dump(res, open(out_path, "w"), indent=2)
