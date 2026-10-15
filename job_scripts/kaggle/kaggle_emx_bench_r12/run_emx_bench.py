"""MICA BENCHMARK rounds 1-2: EModelX(+AF) on 8 cases (frozen settings from the pilot).

Frozen settings from MICA_2026_benchmark_package baseline_runtime/handoff_commands.sh:
commit 53713cf8, --protocol=temp_flex, canonical.mrc (non-cube), deposited FASTA,
exact-sequence AF2 template at <template_dir>/<PDBID_n>/ranked_0.pdb (no AFDB retrieval),
--seed 2022 --mul_proc_num 4 --run_pulchra. Environment mirrors EModelX.yml
(Python 3.8, torch 1.8.1+cu111, pinned scientific stack) in a uv venv.
Deviation: atom3d==0.2.4 is pinned in EModelX.yml but imported by no repository module;
omitted. Feasibility record only: no scoring here.
"""
import glob, hashlib, json, os, shutil, subprocess, threading, time

COMMIT = "53713cf820d10d78b1ca6669ef197850cf907843"
WEIGHTS_GDRIVE_ID = "13BKzEBfL0uubYgcJJTGSZ-PkO9hXAB9X"
IN_DIR = os.path.dirname(next(p for p in glob.glob("/kaggle/input/**/selected12.csv", recursive=True)))      # FASTA
MAP_DIR = os.path.dirname(next(p for p in glob.glob("/kaggle/input/**/bench12/integrity_report.json", recursive=True)))  # maps
PRIORS = sorted(glob.glob("/kaggle/input/**/priors/af2", recursive=True))
WORK = "/tmp/emx"; VENV = f"{WORK}/venv"; REPO = f"{WORK}/EModelX"; SAVE = "/kaggle/working"
CASES = {"EMD-66010_9WIX": "9WIX.fasta", "EMD-54783_9SDB": "9SDB.fasta", "EMD-68274_22HP": "22HP.fasta", "EMD-71359_9P80": "9P80.fasta", "EMD-69851_24VT": "24VT.fasta", "EMD-71407_9P9D": "9P9D.fasta", "EMD-73885_9Z81": "9Z81.fasta", "EMD-63895_9U62": "9U62.fasta"}  # rounds 1-2, small first
rec = {"method": "EModelX(+AF)", "commit": COMMIT, "setup": {}, "cases": {},
       "deviations": ["uv venv instead of conda (same Python minor / torch / CUDA build)",
                      "atom3d omitted: pinned in EModelX.yml, imported by no repository module",
                      "AF2 templates from ColabFold 1.6.3 (see mica-bench-af2-r12)"]}

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def save(): json.dump(rec, open(f"{SAVE}/run_record.json", "w"), indent=2)
def run(name, cmd, env=None, cwd=None, timeout=None):
    t = time.time()
    try:
        p = subprocess.run(cmd, shell=True, env=env, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        code, out, err = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        code, out, err = "timeout", str(e.stdout or ""), str(e.stderr or "")
    dt = round(time.time() - t, 1)
    open(f"{SAVE}/{name}.log", "w").write(out + "\n--- stderr ---\n" + err)
    log(f"{name}: exit {code} in {dt}s")
    if code != 0: print(err[-2500:])
    return code, dt

assert PRIORS, "AF2 priors not mounted"
PRIORS = PRIORS[0]
os.makedirs(WORK, exist_ok=True)
pins = ("numpy==1.19.1 scipy==1.5.2 pandas==1.1.3 biopython==1.78 mrcfile==1.3.0 open3d==0.13.0 "
        "superpose3d==1.1.1 tqdm==4.61.2 scikit-image==0.18.2 networkx==2.6.3 argparse==1.4.0 requests gdown")
for name, cmd in [
    ("install_uv", "pip install -q uv"),
    ("make_venv", f"uv venv -q -p 3.8 {VENV}"),
    ("install_torch", f"uv pip install -q -p {VENV}/bin/python torch==1.8.1+cu111 "
                      f"--find-links https://download.pytorch.org/whl/torch_stable.html"),
    ("install_pins", f"uv pip install -q -p {VENV}/bin/python {pins}"),
    ("clone", f"git clone -q https://github.com/biomed-AI/EModelX.git {REPO} && git -C {REPO} checkout -q {COMMIT}"),
    ("weights", f"{VENV}/bin/gdown -q {WEIGHTS_GDRIVE_ID} -O {WORK}/models.zip && mkdir -p {REPO}/models "
                f"&& unzip -o -q {WORK}/models.zip -d {WORK}/models_unzip "
                f"&& find {WORK}/models_unzip -name 'best_*_model.ckpt' -exec cp {{}} {REPO}/models/ \;"),
    ("cuda_check", f"{VENV}/bin/python -c \"import torch;print(torch.__version__, torch.version.cuda, "
                   f"torch.cuda.is_available(), torch.cuda.get_device_name(0))\""),
]:
    code, dt = run(name, cmd)
    rec["setup"][name] = {"exit": code, "seconds": dt}; save()
    if code != 0:
        rec["outcome"] = f"failed at {name}"; save(); raise SystemExit(1)
rec["setup"]["checkpoints"] = sorted(os.listdir(f"{REPO}/models"))
rec["setup"]["head"] = subprocess.run(f"git -C {REPO} rev-parse HEAD", shell=True, capture_output=True, text=True).stdout.strip()
os.chmod(f"{REPO}/modules/pulchra304/src/pulchra", 0o755)
subprocess.run(f"uv pip freeze -p {VENV}/bin/python > {SAVE}/env_freeze.txt", shell=True)
save()

def atoms(path):
    try: return sum(1 for l in open(path, errors="ignore") if l.startswith(("ATOM", "HETATM")))
    except Exception: return 0

env = dict(os.environ, OMP_NUM_THREADS="4", MKL_NUM_THREADS="4")
for case, fa in CASES.items():
    tdir = f"{PRIORS}/{case}"
    out = f"{WORK}/pred/{case}"
    # Every entity in the FASTA needs its own sequence-validated template; never run a
    # case with a partial template set.
    ents = [b.split("|")[0].strip() for b in open(f"{IN_DIR}/{case}/{fa}").read().split(">")[1:]]
    missing = [e for e in ents if not os.path.exists(f"{tdir}/{e}/ranked_0.pdb")]
    if missing:
        rec["cases"][case] = {"outcome": "skipped_incomplete_templates", "missing": missing}; save(); continue
    os.makedirs(out, exist_ok=True)  # EModelX does not create its output directory
    peak = {"gpu_mib": 0}; stop = threading.Event()
    def sample():
        while not stop.is_set():
            r = subprocess.run("nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits",
                               shell=True, capture_output=True, text=True).stdout.split()
            if r: peak["gpu_mib"] = max(peak["gpu_mib"], max(int(x) for x in r))
            stop.wait(2)
    th = threading.Thread(target=sample, daemon=True); th.start()
    code, dt = run(f"emodelx_{case}",
        f"{VENV}/bin/python run.py --protocol=temp_flex --EM_map {MAP_DIR}/{case}/canonical.mrc "
        f"--fasta {IN_DIR}/{case}/{fa} --template_dir {tdir} --output_dir {out} "
        f"--seed 2022 --mul_proc_num 4 --run_pulchra --pulchra_path modules/pulchra304/src/pulchra",
        env=env, cwd=REPO, timeout=3 * 3600)
    stop.set(); th.join()
    files = {}
    for f in glob.glob(f"{out}/**/*", recursive=True):
        if f.endswith((".pdb", ".cif")):
            files[os.path.relpath(f, out)] = atoms(f)
            d = f"{SAVE}/outputs/{case}"; os.makedirs(d, exist_ok=True); shutil.copy(f, d)
    rec["cases"][case] = {"exit": code, "seconds": dt, "peak_gpu_mib": peak["gpu_mib"],
                          "structure_files_atoms": files}
    save()
rec["outcome"] = "see cases; final-model classification done after inspecting output names"
save()
