"""MICA benchmark, step 3 of MICA_pipeline.sh (run.py) on a Kaggle T4 GPU.

Inputs: input/<EMDB>/ directories produced by steps 1-2 on CPU sessions (mica-bench-dock-*),
i.e. AF3 processing + Phenix domain docking already done. Command = MICA_pipeline.sh step 3 verbatim:
run.py -m MAP -f FASTA -a AF3_results --run_pulchra --pulchra_path --run_phenix --phenix_act --resolution --device cuda
Recorded deviation: batch_threshold raised (README-invited) -> per-sample path (batch 1); avoids T4 OOM;
result-neutral (eval / no_grad). Same environment as the pilots (uv venv py3.8, torch 1.8.1+cu111).
"""
import glob, json, os, shutil, subprocess, threading, time
CASES = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "cases.json"))) \
    if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), "cases.json")) else None
COMMIT = "726d1bb14be04887a43b7f5aefcff0d843761d01"; ZEN = "https://zenodo.org/records/15756654/files"
WORK = "/tmp/mica"; VENV = f"{WORK}/venv"; REPO = f"{WORK}/MICA"; PHX = "/tmp/phenix"; SAVE = "/kaggle/working"
INP = os.path.dirname(next(iter(glob.glob("/kaggle/input/**/mica-bench-mica-inputs/manifest.json", recursive=True)) or
                           glob.glob("/kaggle/input/**/manifest.json", recursive=True)))
PHX_TAR = next(iter(glob.glob("/kaggle/input/**/phenix-installer-1.20.1-4487*.tgz.bin", recursive=True)))
rec = {"method": "MICA", "step": "3 (run.py, GPU)", "commit": COMMIT, "phenix": "1.20.1-4487", "setup": {}, "cases": {}}
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def save(): json.dump(rec, open(f"{SAVE}/run_record.json", "w"), indent=2)
def run(name, cmd, env=None, cwd=None, timeout=None):
    t = time.time()
    try:
        p = subprocess.run(cmd, shell=True, env=env, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        code, out, err = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e: code, out, err = "timeout", str(e.stdout or ""), str(e.stderr or "")
    dt = round(time.time() - t, 1); open(f"{SAVE}/{name}.log", "w").write(out + "\n--- stderr ---\n" + err)
    log(f"{name}: exit {code} in {dt}s")
    if code != 0: print(err[-2500:])
    return code, dt
pins = ("numpy==1.19.1 scipy==1.5.2 pandas==1.1.3 biopython==1.78 mrcfile==1.3.0 open3d==0.13.0 superpose3d==1.1.1 "
        "tqdm==4.61.2 scikit-image==0.18.2 networkx==2.6.3 psutil einops rotary-embedding-torch natsort wandb matplotlib scikit-learn requests")
os.makedirs(WORK, exist_ok=True)
for name, cmd, cwd in [
    ("install_uv", "pip install -q uv", None), ("make_venv", f"uv venv -q -p 3.8 {VENV}", None),
    ("install_env", f"uv pip install -q -p {VENV}/bin/python torch==1.8.1+cu111 {pins} --find-links https://download.pytorch.org/whl/torch_stable.html", None),
    ("clone", f"git clone -q https://github.com/jianlin-cheng/MICA.git {REPO} && git -C {REPO} checkout -q {COMMIT}", None),
    ("patch_batch_threshold", "sed -i 's/self.batch_threshold = 200/self.batch_threshold = 10**9  # per-sample path (batch 1); avoids CUDA OOM on 15 GB T4; result-neutral in eval\\/no_grad/' utils/predict.py && grep -q 'batch_threshold = 10\\*\\*9' utils/predict.py", REPO),
    ("download_models", f"curl -sfL '{ZEN}/trained_models.tar.gz?download=1' -o tm.tgz && tar -xzf tm.tgz && rm tm.tgz", REPO),
    ("unpack_phenix", f"mkdir -p {WORK}/phx_src && tar -xzf {PHX_TAR} -C {WORK}/phx_src", None),
    ("install_phenix", f"cd {WORK}/phx_src/phenix-installer-1.20.1-4487-* && ./install --prefix={PHX} --nproc=4", None),
    ("cuda_check", f"{VENV}/bin/python -c \"import torch;assert torch.cuda.is_available();print(torch.__version__, torch.cuda.get_device_name(0))\"", None)]:
    code, dt = run(name, cmd, cwd=cwd); rec["setup"][name] = {"exit": code, "seconds": dt}; save()
    if code != 0: rec["outcome"] = f"failed at {name}"; save(); raise SystemExit(1)
os.chmod(f"{REPO}/modules/pulchra304/src/pulchra", 0o755)
phenix_env = glob.glob(f"{PHX}/phenix-1.20.1-4487/phenix_env.sh")[0]
subprocess.run(f"uv pip freeze -p {VENV}/bin/python > {SAVE}/env_freeze.txt", shell=True)
man = json.load(open(f"{INP}/manifest.json"))
env = dict(os.environ, PATH=f"{VENV}/bin:" + os.environ["PATH"], OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", WANDB_MODE="disabled")
for emd in CASES:
    m = man[emd]; tars = sorted(glob.glob(f"/kaggle/input/**/mica_input_{emd}.tar.gz", recursive=True))
    if not tars: rec["cases"][emd] = {"case": m["case"], "outcome": "missing_docked_input"}; save(); continue
    os.makedirs(f"{REPO}/input", exist_ok=True); subprocess.run(f"tar -xzf {tars[0]} -C {REPO}/input", shell=True, check=True)
    d = f"input/{emd}"; docked = glob.glob(f"{REPO}/{d}/{emd}_af3_docked.pdb")
    peak = {"gpu_mib": 0}; stop = threading.Event()
    def sample():
        while not stop.is_set():
            r = subprocess.run("nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits", shell=True, capture_output=True, text=True).stdout.split()
            if r: peak["gpu_mib"] = max(peak["gpu_mib"], max(int(x) for x in r))
            stop.wait(5)
    th = threading.Thread(target=sample, daemon=True); th.start()
    out = f"{WORK}/out_{emd}"
    code, dt = run(f"step3_{emd}", f"python run.py -m {d}/emd_{emd}.map -f {d}/{m['fasta']} -a {d}/AF3_results --run_pulchra "
                   f"--pulchra_path=modules/pulchra304/src/pulchra --run_phenix --phenix_act={phenix_env} --resolution={m['resolution_A']} "
                   f"--device=cuda -o {out}", env=env, cwd=REPO, timeout=4 * 3600)
    stop.set(); th.join()
    outdir = f"{SAVE}/outputs/{emd}"; os.makedirs(outdir, exist_ok=True); files = {}
    for f in glob.glob(f"{out}/*"):
        if os.path.isfile(f):
            shutil.copy(f, outdir)
            files[os.path.basename(f)] = sum(1 for l in open(f, errors="ignore") if l.startswith(("ATOM", "HETATM"))) if f.endswith((".pdb", ".cif")) else None
    final = [k for k in files if k.endswith("_MICA_all_atom_model.pdb")]; atoms = files[final[0]] if final else 0
    rec["cases"][emd] = {"case": m["case"], "af3_docked_present": bool(docked), "exit": code, "seconds": dt, "peak_gpu_mib": peak["gpu_mib"],
                         "outputs": files, "final_model": final[0] if final else None, "final_atoms": atoms,
                         "outcome": "non_empty_model" if code == 0 and atoms > 0 else "empty_model" if code == 0 else "failed"}
    save(); log(f"{emd}: {rec['cases'][emd]['outcome']} atoms={atoms} {dt}s")
    shutil.rmtree(f"{REPO}/{d}", ignore_errors=True)
rec["outcome"] = "done"; save()
