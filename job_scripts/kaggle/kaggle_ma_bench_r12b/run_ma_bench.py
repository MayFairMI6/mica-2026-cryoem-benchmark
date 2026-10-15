"""MICA BENCHMARK rounds 1-2: ModelAngelo 1.0.18 on 8 cases (frozen settings from the pilot).

Same frozen settings as step 1 (commit 994945bd, `nucleotides` bundle, canonical_cube.mrc,
deposited FASTA, `build ... --device 0`; uv venv with Python 3.11 / torch 2.9.1).
One weight download, amortized across cases. Outcome is judged on the FINAL pruned
model only. Feasibility record: no scoring, no ranking.
"""
import glob, hashlib, json, os, resource, shutil, subprocess, threading, time

CASES = {"EMD-71359_9P80": "9P80.fasta", "EMD-63895_9U62": "9U62.fasta", "EMD-66010_9WIX": "9WIX.fasta", "EMD-69851_24VT": "24VT.fasta", "EMD-73885_9Z81": "9Z81.fasta", "EMD-54783_9SDB": "9SDB.fasta", "EMD-68274_22HP": "22HP.fasta", "EMD-71407_9P9D": "9P9D.fasta"}
COMMIT = "994945bdfa6e5368e0d62349a47792f4864eebc3"
IN_DIR = os.path.dirname(next(p for p in glob.glob("/kaggle/input/**/selected12.csv", recursive=True)))      # FASTA
MAP_DIR = os.path.dirname(next(p for p in glob.glob("/kaggle/input/**/bench12/integrity_report.json", recursive=True)))  # maps
WORK = "/tmp/mica"; VENV = f"{WORK}/ma-venv"; TORCH_HOME = f"{WORK}/weights/modelangelo"
SAVE = "/kaggle/working"
rec = {"method": "ModelAngelo", "commit": COMMIT, "setup": {}, "cases": {}}

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def save():
    json.dump(rec, open(f"{SAVE}/run_record.json", "w"), indent=2)

def run(name, cmd, env, timeout=None):
    log(f"START {name}")
    t = time.time()
    try:
        p = subprocess.run(cmd, shell=True, env=env, capture_output=True, text=True, timeout=timeout)
        code, out, err = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        code, out, err = "timeout", str(e.stdout or ""), str(e.stderr or "")
    dt = round(time.time() - t, 1)
    open(f"{SAVE}/{name}.log", "w").write(out + "\n--- stderr ---\n" + err)
    log(f"END {name}: exit {code} in {dt}s")
    if code != 0:
        print(err[-3000:])
    return code, dt

os.makedirs(WORK, exist_ok=True)
integ = json.load(open(f"{MAP_DIR}/integrity_report.json"))["cases"]
for case in CASES:
    assert integ[case]["full_map_integrity"] == "pass", case
    got = hashlib.sha256(open(f"{MAP_DIR}/{case}/canonical_cube.mrc", "rb").read()).hexdigest()
    assert got == integ[case]["canonical_cube"]["destination_sha256"], f"cube hash mismatch {case}"
rec["inputs_verified"] = True
rec["gpu"] = subprocess.run("nvidia-smi --query-gpu=name,memory.total --format=csv,noheader",
                            shell=True, capture_output=True, text=True).stdout.strip()

env = dict(os.environ, TORCH_HOME=TORCH_HOME, OMP_NUM_THREADS="4", MKL_NUM_THREADS="4")
for name, cmd in [
    ("install_uv", "pip install -q uv"),
    ("make_venv", f"uv venv -q -p 3.11 {VENV}"),
    ("install_torch", f"uv pip install -q -p {VENV}/bin/python torch==2.9.1 torchvision"),
    ("install_modelangelo", f"uv pip install -q -p {VENV}/bin/python "
                            f"git+https://github.com/3dem/model-angelo@{COMMIT}"),
    ("setup_weights", f"{VENV}/bin/model_angelo setup_weights --bundle-name nucleotides"),
]:
    code, dt = run(name, cmd, env)
    rec["setup"][name] = {"exit": code, "seconds": dt}
    if code != 0:
        rec["outcome"] = f"failed at {name}"; save(); raise SystemExit(1)
subprocess.run(f"uv pip freeze -p {VENV}/bin/python > {SAVE}/env_freeze.txt", shell=True)
rec["ma_version"] = subprocess.run(f"{VENV}/bin/model_angelo --version", shell=True,
                                   capture_output=True, text=True, env=env).stdout.strip()
save()

def atoms_in(path):
    return sum(1 for l in open(path) if l.startswith(("ATOM", "HETATM"))) if os.path.exists(path) else 0

for case, fa in CASES.items():
    out = f"{WORK}/pred/{case}"
    peak = {"gpu_mib": 0}; stop = threading.Event()
    def sample():
        while not stop.is_set():
            r = subprocess.run("nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits",
                               shell=True, capture_output=True, text=True).stdout.split()
            if r: peak["gpu_mib"] = max(peak["gpu_mib"], max(int(x) for x in r))
            stop.wait(2)
    th = threading.Thread(target=sample, daemon=True); th.start()
    code, dt = run(f"build_{case}", f"{VENV}/bin/model_angelo build -v {MAP_DIR}/{case}/canonical_cube.mrc "
                                    f"-pf {IN_DIR}/{case}/{fa} -o {out} --device 0", env, timeout=3 * 3600)
    stop.set(); th.join()
    final, raw = atoms_in(f"{out}/{case}.cif"), atoms_in(f"{out}/{case}_raw.cif")
    for f in (f"{out}/{case}.cif", f"{out}/{case}_raw.cif"):
        if os.path.exists(f): shutil.copy(f, SAVE)
    rec["cases"][case] = {
        "exit": code, "build_seconds": dt, "peak_gpu_mib": peak["gpu_mib"],
        "final_atoms": final, "raw_atoms": raw,
        "outcome": "non_empty_model" if code == 0 and final > 0
                   else "empty_model" if code == 0 else "build_failed",
    }
    save()
    log(f"{case}: {rec['cases'][case]['outcome']} | final atoms {final} | raw {raw} | {dt}s")

rec["peak_child_rss_gib_session"] = round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1048576, 2)
save()
