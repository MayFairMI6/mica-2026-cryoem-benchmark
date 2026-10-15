"""CryoZeta inputs for the 24 benchmark cases (post hoc addition, 2026-10-12). Runs locally on macOS.
  1. per case: entities (sequence, copies) from the same FASTA every other method received; total residues;
     EMDB recommended contour level (EMDB API); resolution from case_rounds.csv
  2. per unique protein sequence: MSA from the public ColabFold MMseqs2 server (mode env), saved in the layout of
     CryoZeta's own example: <sha256>/uniref100_hits.a3m (= uniref.a3m) and mmseqs_other_hits.a3m (= environmental db)
  3. one input JSON per case (format of assets/examples/example.json)
  4. upload to GAIVI, hard-link the canonical maps as <PDB>.map, submit cz_run.sbatch in 3 lanes, smallest case first
Usage: cz_prepare.py plan | all        (plan = steps 1 only, prints the table)"""
import csv, hashlib, io, json, os, re, subprocess, sys, tarfile, time, urllib.parse, urllib.request
R = "/Users/anony/Downloads/MICA_GAIVI_run"; B = f"{R}/benchmark"; D = f"{B}/cryozeta"; os.makedirs(f"{D}/msa", exist_ok=True); os.makedirs(f"{D}/json", exist_ok=True)
UA = {"User-Agent": "colabfold/1.5.5 (academic benchmark; low volume)"}
def log(*a): print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, flush=True)
def get(url, data=None, timeout=120):
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode() if data else None, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r: return r.read()
rounds = {r["case_id"]: r for r in csv.DictReader(open(f"{B}/case_rounds.csv"))}
cases = []
for c in sorted(rounds, key=lambda c: (int(rounds[c]["round"]), float(rounds[c]["resolution_A"]))):
    pdb = c.split("_")[1]; emd = c.split("_")[0].replace("EMD-", "")
    fa = next(p for p in (f"{B}/kaggle_bench12_inputs/{c}/{pdb}.fasta", f"{B}/ext_inputs/{c}/{pdb}.fasta") if os.path.exists(p))
    ents = []
    for b in open(fa).read().strip().split(">")[1:]:
        h, *s = b.split("\n"); seq = "".join(s).strip().upper()
        ents.append({"id": h.split("|")[0], "copies": len(re.sub(r"Chains?\s*", "", h.split("|")[1]).split(",")), "seq": seq,
                     "nucleic": set(seq) <= set("ACGUTN"), "sha": hashlib.sha256(seq.encode()).hexdigest()})
    cases.append({"case": c, "pdb": pdb, "emd": emd, "res": float(rounds[c]["resolution_A"]), "ents": ents, "tokens": sum(len(e["seq"]) * e["copies"] for e in ents)})
cache = f"{D}/contours.json"; contours = json.load(open(cache)) if os.path.exists(cache) else {}
for k in cases:
    if k["emd"] not in contours:
        j = json.loads(get(f"https://www.ebi.ac.uk/emdb/api/entry/EMD-{k['emd']}"))
        cl = j["map"]["contour_list"]["contour"]; pri = [x for x in cl if x.get("primary")] or cl
        contours[k["emd"]] = float(pri[0]["level"]); json.dump(contours, open(cache, "w"), indent=1)
    k["contour"] = contours[k["emd"]]; k["mode"] = "normal" if k["tokens"] <= 2800 else "large"
print("case   EMD     res  contour  entities  residues(seq x copies)  mode   nucleic")
for k in sorted(cases, key=lambda k: k["tokens"]):
    print(f"{k['pdb']:6s} {k['emd']:6s} {k['res']:.2f} {k['contour']:8.4g} {len(k['ents']):6d} {k['tokens']:12d}            {k['mode']:6s} {[e['id'] for e in k['ents'] if e['nucleic']]}")
uniq = {e["sha"]: e["seq"] for k in cases for e in k["ents"] if not e["nucleic"]}
print("unique protein sequences:", len(uniq), "| cases normal/large:", sum(k["mode"] == "normal" for k in cases), sum(k["mode"] == "large" for k in cases))
json.dump([{x: k[x] for x in ("case", "pdb", "emd", "res", "contour", "tokens", "mode")} for k in cases], open(f"{D}/plan.json", "w"), indent=1)
if len(sys.argv) < 2 or sys.argv[1] != "all": sys.exit(0)
if any(e["nucleic"] for k in cases for e in k["ents"]): log("WARNING: nucleic-acid entities present; they are written as dnaSequence/rnaSequence without MSA")

CZ = "/home/a/aravindav/bench/envs/cz"   # absolute paths: inside the container relative paths resolve against /app/CryoZeta, not the run directory
API = "https://api.colabfold.com"
def msa(sha, seq):
    d = f"{D}/msa/{sha}"
    if all(os.path.getsize(f"{d}/{f}") > 0 for f in ("uniref100_hits.a3m", "mmseqs_other_hits.a3m") if os.path.exists(f"{d}/{f}")) and os.path.exists(f"{d}/mmseqs_other_hits.a3m") and os.path.exists(f"{d}/uniref100_hits.a3m"): return "cached"
    for attempt in range(8):
        try:
            t = json.loads(get(f"{API}/ticket/msa", {"q": f">101\n{seq}\n", "mode": "env"}))
            while t.get("status") in ("PENDING", "RUNNING", "UNKNOWN", "RATELIMIT"):
                time.sleep(60 if t["status"] == "RATELIMIT" else 8)
                t = json.loads(get(f"{API}/ticket/msa", {"q": f">101\n{seq}\n", "mode": "env"})) if t["status"] == "RATELIMIT" else json.loads(get(f"{API}/ticket/{t['id']}"))
            if t.get("status") != "COMPLETE": raise RuntimeError(f"status {t.get('status')}")
            tf = tarfile.open(fileobj=io.BytesIO(get(f"{API}/result/download/{t['id']}", timeout=600)))
            os.makedirs(d, exist_ok=True); names = tf.getnames()
            for src, dst in (("uniref.a3m", "uniref100_hits.a3m"), ("bfd.mgnify30.metaeuk30.smag30.a3m", "mmseqs_other_hits.a3m")):
                txt = tf.extractfile(src).read().decode().replace("\x00", "")
                lines = txt.split("\n"); assert lines[0].startswith(">") and lines[1].replace("-", "") == seq, "query mismatch"
                lines[0] = ">" + sha; open(f"{d}/{dst}", "w").write("\n".join(lines).rstrip("\n") + "\n")
            return f"fetched ({sum(1 for l in open(f'{d}/uniref100_hits.a3m') if l[0] == '>')} + {sum(1 for l in open(f'{d}/mmseqs_other_hits.a3m') if l[0] == '>')} seqs)"
        except Exception as e:
            log("  retry", attempt + 1, sha[:10], repr(e)[:150]); time.sleep(30 * (attempt + 1))
    return "FAILED"
bad = []
for i, (sha, seq) in enumerate(uniq.items()):
    r = msa(sha, seq); log(f"MSA {i + 1}/{len(uniq)} {sha[:10]} len {len(seq)}: {r}")
    if r == "FAILED": bad.append(sha)
    elif r != "cached": time.sleep(5)
ready = []
for k in cases:
    if any(e["sha"] in bad for e in k["ents"]): log("no JSON for", k["pdb"], "(MSA missing)"); continue
    seqs = []
    for e in k["ents"]:
        if e["nucleic"]: seqs.append({("rnaSequence" if "U" in e["seq"] else "dnaSequence"): {"sequence": e["seq"], "count": e["copies"]}})
        else: seqs.append({"proteinChain": {"sequence": e["seq"], "count": e["copies"], "msa": {"precomputed_msa_dir": f"{CZ}/inputs/msa/{e['sha']}", "pairing_db": "uniref100"}}})
    json.dump([{"name": k["pdb"], "modelSeeds": [], "map_path": f"{CZ}/inputs/maps/{k['pdb']}.map", "resolution": k["res"], "contour_level": k["contour"], "sequences": seqs}],
              open(f"{D}/json/{k['pdb']}.json", "w"), indent=1); ready.append(k)
log("JSON written for", len(ready), "cases; uploading")
SSH = ["ssh", "-o", f"ControlPath={os.path.expanduser('~/.ssh/cm-gaivi')}", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", "aravindav@gaivi.cse.usf.edu"]
tar = subprocess.Popen(["tar", "-C", D, "-cf", "-", "msa", "json", "plan.json"], stdout=subprocess.PIPE)
up = subprocess.run(SSH + ["mkdir -p ~/bench/envs/cz/inputs && tar -C ~/bench/envs/cz/inputs -xf - && echo UPLOADED"], stdin=tar.stdout, capture_output=True, text=True)
if "UPLOADED" not in up.stdout: log("UPLOAD FAILED (log in to GAIVI, then rerun: python3 automation/cz_prepare.py all)", up.stderr[-300:]); sys.exit(3)
order = sorted(ready, key=lambda k: k["tokens"]); cmds = ["cd ~/bench", "mkdir -p envs/cz/inputs/maps out/CryoZeta"]
for i, k in enumerate(order):
    cmds.append(f"ln -f maps/{k['case']}/canonical.mrc envs/cz/inputs/maps/{k['pdb']}.map")
    cmds.append(f"[ -f out/CryoZeta/{k['case']}/run_record.json ] || squeue -u $USER -h -o '%k' | grep -q 'cz:{k['case']}$' || sbatch -J cz{i % 3 + 1} --comment=cz:{k['case']} --dependency=singleton cz_run.sbatch {k['case']} {k['mode']} | tee -a chain_ids.txt")
cmds.append("echo \"CryoZeta chain submitted $(date -Is)\" >> chain_ids.txt; squeue -u $USER -h | wc -l")
sub = subprocess.run(SSH + ["; ".join(cmds)], capture_output=True, text=True); log("SUBMIT:", sub.stdout[-900:].replace("\n", " | "), sub.stderr[-200:])
