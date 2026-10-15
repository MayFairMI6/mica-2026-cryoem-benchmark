"""MICA benchmark, steps 1-2 of MICA_pipeline.sh on a CPU session (no GPU quota).

Step 1: utils/process_AF3_results.py -f FASTA -a AF3_results
Step 2: utils/dock_in_map.py -m MAP -c CONTOUR -r RES -f FASTA -a AF3_results --phenix_act PHENIX
Exactly the commands MICA_pipeline.sh runs (commit 726d1bb1); step 3 (run.py) runs later on a GPU.
Map = the project's canonical map (lossless canonicalization, frozen and integrity-checked).
Outputs: the whole input/<EMDB>/ directory after steps 1-2, tarred, for the GPU step.
"""
import glob, json, os, shutil, subprocess, time
CASES = __CASES__
COMMIT = "726d1bb14be04887a43b7f5aefcff0d843761d01"
WORK = "/tmp/mica"; VENV = f"{WORK}/venv"; REPO = f"{WORK}/MICA"; PHX = "/tmp/phenix"; SAVE = "/kaggle/working"
INP = os.path.dirname(next(iter(glob.glob("/kaggle/input/**/mica-bench-mica-inputs/manifest.json", recursive=True)) or
                           glob.glob("/kaggle/input/**/manifest.json", recursive=True)))
MAPS = os.path.dirname(next(iter(glob.glob("/kaggle/input/**/bench12/integrity_report.json", recursive=True))))
PHX_TAR = next(iter(glob.glob("/kaggle/input/**/phenix-installer-1.20.1-4487*.tgz.bin", recursive=True)))
rec = {"step": "MICA steps 1-2 (CPU)", "commit": COMMIT, "setup": {}, "cases": {}}
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def save(): json.dump(rec, open(f"{SAVE}/run_record.json", "w"), indent=2)
def run(name, cmd, cwd=None, env=None, timeout=None):
    t = time.time()
    try:
        p = subprocess.run(cmd, shell=True, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout); code, out, err = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e: code, out, err = "timeout", str(e.stdout or ""), str(e.stderr or "")
    dt = round(time.time() - t, 1); open(f"{SAVE}/{name}.log", "w").write(out + "\n--- stderr ---\n" + err)
    log(f"{name}: exit {code} in {dt}s"); 
    if code != 0: print(err[-2000:])
    return code, dt
pins = ("numpy==1.19.1 scipy==1.5.2 pandas==1.1.3 biopython==1.78 mrcfile==1.3.0 open3d==0.13.0 superpose3d==1.1.1 "
        "tqdm==4.61.2 scikit-image==0.18.2 networkx==2.6.3 psutil einops rotary-embedding-torch natsort wandb matplotlib scikit-learn requests")
os.makedirs(WORK, exist_ok=True)
for name, cmd in [("install_uv", "pip install -q uv"), ("make_venv", f"uv venv -q -p 3.8 {VENV}"),
                  ("install_env", f"uv pip install -q -p {VENV}/bin/python torch==1.8.1+cu111 {pins} --find-links https://download.pytorch.org/whl/torch_stable.html"),
                  ("clone", f"git clone -q https://github.com/jianlin-cheng/MICA.git {REPO} && git -C {REPO} checkout -q {COMMIT}"),
                  ("unpack_phenix", f"mkdir -p {WORK}/phx_src && tar -xzf {PHX_TAR} -C {WORK}/phx_src"),
                  ("install_phenix", f"cd {WORK}/phx_src/phenix-installer-1.20.1-4487-* && ./install --prefix={PHX} --nproc=4")]:
    code, dt = run(name, cmd); rec["setup"][name] = {"exit": code, "seconds": dt}; save()
    if code != 0: rec["outcome"] = f"failed at {name}"; save(); raise SystemExit(1)
phenix_env = glob.glob(f"{PHX}/phenix-1.20.1-4487/phenix_env.sh")[0]
man = json.load(open(f"{INP}/manifest.json"))
env = dict(os.environ, PATH=f"{VENV}/bin:" + os.environ["PATH"], OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", WANDB_MODE="disabled")
for emd in CASES:
    m = man[emd]; d = f"{REPO}/input/{emd}"; shutil.copytree(f"{INP}/{emd}", d)
    shutil.copy(f"{MAPS}/{m['case']}/canonical.mrc", f"{d}/emd_{emd}.map")
    fa, af = f"input/{emd}/{m['fasta']}", f"input/{emd}/AF3_results"
    c1, t1 = run(f"step1_{emd}", f"python utils/process_AF3_results.py -f {fa} -a {af}", cwd=REPO, env=env, timeout=2 * 3600)
    c2, t2 = (run(f"step2_{emd}", f"python utils/dock_in_map.py -m input/{emd}/emd_{emd}.map -c {m['contour_level']} -r {m['resolution_A']} "
                  f"-f {fa} -a {af} --phenix_act {phenix_env}", cwd=REPO, env=env, timeout=10 * 3600) if c1 == 0 else ("skipped", 0))
    shutil.make_archive(f"{SAVE}/mica_input_{emd}", "gztar", f"{REPO}/input", emd)
    rec["cases"][emd] = {"case": m["case"], "step1": {"exit": c1, "seconds": t1}, "step2": {"exit": c2, "seconds": t2},
                         "contents": sorted(os.listdir(d))}
    save()
rec["outcome"] = "done"; save()
