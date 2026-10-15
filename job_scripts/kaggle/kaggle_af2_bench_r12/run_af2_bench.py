"""EModelX(+AF) prerequisite: exact-sequence AlphaFold2 templates, BENCHMARK rounds 1-2 (8 cases).

The frozen protocol stages one AF2 single-chain model per entity at
priors/af2/<case>/<PDBID_n>/ranked_0.pdb and sequence-validates it, replacing EModelX's
default 60%-identity AFDB retrieval. Generated with ColabFold 1.6.3 (AlphaFold2 weights,
MMseqs2 MSA server), 5 models, top-ranked kept.
v2: driver-only LD_LIBRARY_PATH (v1 silently ran on CPU) + hard GPU check. Two T4s run jobs in parallel.
"""
import glob, json, os, queue, shutil, subprocess, threading, time

CF_VERSION = "1.6.3"
INPUTS = next(p for p in glob.glob("/kaggle/input/**/selected12.csv", recursive=True))
IN_DIR = os.path.dirname(INPUTS)
WORK = "/tmp/af2"; VENV = f"{WORK}/cf-venv"; SAVE = "/kaggle/working"
CASES = {"EMD-71359_9P80": "9P80.fasta", "EMD-63895_9U62": "9U62.fasta", "EMD-66010_9WIX": "9WIX.fasta", "EMD-69851_24VT": "24VT.fasta", "EMD-73885_9Z81": "9Z81.fasta", "EMD-54783_9SDB": "9SDB.fasta", "EMD-68274_22HP": "22HP.fasta", "EMD-71407_9P9D": "9P9D.fasta"}  # rounds 1-2
rec = {"tool": f"ColabFold {CF_VERSION}", "settings": {"num_models": 5, "num_recycle": 3,
       "random_seed": 0, "msa": "MMseqs2 server (api.colabfold.com)"}, "setup": {}, "jobs": {}}

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def save(): json.dump(rec, open(f"{SAVE}/run_record.json", "w"), indent=2)

def run(cmd, env=None, timeout=None, logname=None):
    t = time.time()
    try:
        p = subprocess.run(cmd, shell=True, env=env, capture_output=True, text=True, timeout=timeout)
        code, out, err = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        code, out, err = "timeout", str(e.stdout or ""), str(e.stderr or "")
    if logname:
        open(f"{SAVE}/{logname}.log", "w").write(out + "\n--- stderr ---\n" + err)
    return code, round(time.time() - t, 1), err

os.makedirs(WORK, exist_ok=True)
for name, cmd in [
    ("install_uv", "pip install -q uv"),
    ("make_venv", f"uv venv -q -p 3.11 {VENV}"),
    ("install_colabfold", f'uv pip install -q -p {VENV}/bin/python "colabfold[alphafold]=={CF_VERSION}" '
                          f'"jax[cuda12]>=0.6.2,<0.12"'),
]:
    code, dt, err = run(cmd, logname=name)
    rec["setup"][name] = {"exit": code, "seconds": dt}; log(f"{name}: exit {code} ({dt}s)")
    if code != 0:
        print(err[-3000:]); rec["outcome"] = f"failed at {name}"; save(); raise SystemExit(1)
run(f"uv pip freeze -p {VENV}/bin/python > {SAVE}/env_freeze.txt")
gpu_env = dict(os.environ, LD_LIBRARY_PATH="/usr/local/nvidia/lib64")
chk = subprocess.run(f"{VENV}/bin/python -c 'import jax; print(jax.default_backend(), jax.devices())'",
                     shell=True, capture_output=True, text=True, env=gpu_env)
rec["jax_backend_check"] = chk.stdout.strip() + " | " + chk.stderr.strip()[-300:]; save()
log("JAX check: " + rec["jax_backend_check"][:200])
if not chk.stdout.startswith("gpu"):
    rec["outcome"] = "aborted: JAX sees no GPU"; save(); raise SystemExit(1)

# One single-sequence FASTA per entity, named by the header ID before '|'
jobs = []
for case, fa in CASES.items():
    text = open(f"{IN_DIR}/{case}/{fa}").read().strip().split(">")[1:]
    for block in text:
        head, *seq = block.strip().split("\n"); seq = "".join(seq).strip()
        eid = head.split("|")[0]
        d = f"{WORK}/in/{case}"; os.makedirs(d, exist_ok=True)
        open(f"{d}/{eid}.fasta", "w").write(f">{eid}\n{seq}\n")
        jobs.append((case, eid, seq))
jobs.sort(key=lambda j: -len(j[2]))  # longest first so the two GPUs finish together
rec["n_jobs"] = len(jobs); save()

def one(case, eid, seq, gpu):
    env = dict(os.environ, LD_LIBRARY_PATH="/usr/local/nvidia/lib64", CUDA_VISIBLE_DEVICES=str(gpu), TF_FORCE_UNIFIED_MEMORY="1",
               XLA_PYTHON_CLIENT_MEM_FRACTION="4.0", XDG_CACHE_HOME=f"{WORK}/cache")
    out = f"{WORK}/out/{case}/{eid}"
    code, dt, err = run(f"{VENV}/bin/colabfold_batch --num-models 5 --num-recycle 3 --random-seed 0 "
                        f"{WORK}/in/{case}/{eid}.fasta {out}", env=env, timeout=5 * 3600,
                        logname=f"colabfold_{eid}")
    pdbs = sorted(glob.glob(f"{out}/{eid}*_rank_001_*.pdb"))
    r = {"case": case, "length": len(seq), "gpu": gpu, "exit": code, "seconds": dt}
    if code == 0 and pdbs:
        dst = f"{SAVE}/priors/af2/{case}/{eid}"; os.makedirs(dst, exist_ok=True)
        shutil.copy(pdbs[0], f"{dst}/ranked_0.pdb")
        # sequence validation against the frozen FASTA
        three = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
                 "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
        got, last, plddt = [], None, []
        for l in open(pdbs[0]):
            if l.startswith("ATOM") and l[12:16].strip() == "CA":
                got.append(three.get(l[17:20], "X")); plddt.append(float(l[60:66]))
        r.update(source=os.path.basename(pdbs[0]), sequence_matches_fasta="".join(got) == seq,
                 mean_plddt=round(sum(plddt) / len(plddt), 1) if plddt else None)
        for s in glob.glob(f"{out}/{eid}*scores_rank_001*.json"):
            shutil.copy(s, dst)
    else:
        r["error_tail"] = err[-1500:]
    rec["jobs"][eid] = r; save()
    log(f"{eid} ({len(seq)} aa, GPU{gpu}): exit {code} in {dt}s, seq match {r.get('sequence_matches_fasta')}, pLDDT {r.get('mean_plddt')}")

q = queue.Queue()
for j in jobs: q.put(j)
def worker(gpu):
    while True:
        try: case, eid, seq = q.get_nowait()
        except queue.Empty: return
        one(case, eid, seq, gpu)
threads = [threading.Thread(target=worker, args=(g,)) for g in (0, 1)]
for t in threads: t.start()
for t in threads: t.join()

ok = [e for e, r in rec["jobs"].items() if r.get("sequence_matches_fasta")]
rec["outcome"] = f"{len(ok)}/{len(jobs)} templates produced and sequence-validated"
save(); log(rec["outcome"])
shutil.make_archive(f"{SAVE}/priors_af2", "zip", f"{SAVE}/priors")
