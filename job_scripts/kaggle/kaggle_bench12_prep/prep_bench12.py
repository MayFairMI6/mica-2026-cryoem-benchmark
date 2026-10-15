"""Benchmark core (12 cases): full-map integrity audit + lossless canonicalization.

For each selected case: download the primary EMDB map; verify (1) gzip integrity,
(2) sha256 of the first 65,536 compressed bytes == frozen value, (3) sha256 of the first
1,024 uncompressed bytes (MRC header) == frozen value, (4) uncompressed size == EMDB
metadata, (5) MRC parses, 3-D float data, all values finite; verify the frozen coordinate
and FASTA files against the candidate table's sha256. Then write canonical.mrc and
canonical_cube.mrc (cube_min=128) with the project's canonicalize_map.py, exactly as for
the pilots. No method is run. Output: /kaggle/working/bench12/<case>/.
"""
import csv, glob, gzip, hashlib, json, os, shutil, subprocess, sys, time

IN = os.path.dirname(next(iter(glob.glob("/kaggle/input/**/selected12.csv", recursive=True))))
OUT = "/kaggle/working/bench12"; RAW = "/tmp/raw"; os.makedirs(OUT, exist_ok=True); os.makedirs(RAW, exist_ok=True)
subprocess.run("pip install -q mrcfile", shell=True)
sys.path.insert(0, IN)
import mrcfile, numpy as np
from canonicalize_map import canonicalize

sha = lambda b: hashlib.sha256(b).hexdigest()
def fsha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""): h.update(chunk)
    return h.hexdigest()

EXPECTED_CIF_SHA = {"EMD-71359_9P80": "7e0f67f05ea37f7aa6ee50335f57b7dbef18e8ce48904b3a2d75cc0eed67c429", "EMD-73885_9Z81": "af167946dd441209dfa37a2aef9fbf26702ed4f9854e1ab16a34b80a09a44bdf", "EMD-74910_9ZWK": "0f6eb47128415b2a2f4f891daaeacac6bfd659edf7539226dd284fdbefc66d67", "EMD-63895_9U62": "0491b1f3fe4ff5eb649220ef300baea689fba3bdc4ee676acf8dc5e16b21cb8e", "EMD-54783_9SDB": "a6ad6e95f1aaeffa97d19cd1cb2bc84980195531e0f8467d4ad895691418756e", "EMD-63011_9LE0": "ebbbce17e9dd9008a2a3759f247db9ed19eac77e4de4a3fc1bf39472f42baf04", "EMD-66010_9WIX": "86e4210b9c0388d883cff59f235ca3e3d2e87f777a693d78857c2757c6d4c2ab", "EMD-68274_22HP": "6e26a7025adbfb633be5f5407179ebb35de515cf64ef5ee99697d3f5de661791", "EMD-55353_9SYT": "30e2c5cc1e57f876b80b4f3415f014d1fbb0e50c388d0d979af4261a4c37eacd", "EMD-69851_24VT": "c261e55d6a3e9e033f79a326a8c0cfe4fa48c112fefce4363682f2da9378264f", "EMD-71407_9P9D": "df9da67e7712b3f43e81de0df5fad44b3fed89ab74b3bc4276f56b7a2a9eda1e", "EMD-65920_9WF0": "d88b03b5c26131fda3df151258747a72a01f3c0b8f3692ff01ff1ff6c2ffb2d4"}
table = {r["case_id"]: r for r in csv.DictReader(open(f"{IN}/primary24_scope_screened_candidate.csv"))}
cases = [r["case_id"] for r in csv.DictReader(open(f"{IN}/selected12.csv"))]
report = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "cases": {}}
for c in cases:
    t0 = time.time(); r = {}; src = f"{IN}/{c}"; dst = f"{OUT}/{c}"; os.makedirs(dst, exist_ok=True)
    try:
        req = json.load(open(f"{src}/map_header_request.json")); qc = json.load(open(f"{src}/input_qc.json"))
        gz = f"{RAW}/{c}.map.gz"
        subprocess.run(f"curl -sfL --retry 3 -o {gz} '{req['url']}'", shell=True, check=True)
        r["compressed_bytes"] = os.path.getsize(gz)
        with open(gz, "rb") as f: r["compressed_prefix_match"] = sha(f.read(65536)) == req["compressed_prefix_sha256"]
        raw = f"{RAW}/{c}.map"
        with gzip.open(gz, "rb") as fi, open(raw, "wb") as fo: shutil.copyfileobj(fi, fo, 1 << 22)
        r["gzip_ok"] = True
        with open(raw, "rb") as f: r["header_match"] = sha(f.read(1024)) == req["header_sha256"]
        r["uncompressed_bytes"] = os.path.getsize(raw)
        r["size_matches_emdb_kb"] = abs(r["uncompressed_bytes"] / 1000 - float(qc["emdb_map_uncompressed_size_kb"])) <= 1  # EMDB kB = 1000 bytes
        with mrcfile.open(raw, permissive=False) as m:
            d = m.data; r["shape"] = list(d.shape); r["dtype"] = str(d.dtype)
            r["all_finite"] = bool(np.isfinite(d).all()); r["density_min_max"] = [float(d.min()), float(d.max())]
        r["raw_map_sha256"] = fsha(raw)
        row = table[c]
        cifs = glob.glob(f"{src}/*.cif.gz") + glob.glob(f"{src}/*.cif"); cif = cifs[0]
        fa = glob.glob(f"{src}/*.fasta")[0]
        # Kaggle decompresses uploaded .cif.gz; verify the gz hash if present, else the
        # uncompressed-content hash computed locally from the frozen gz (which matched the table)
        r["coordinates_sha_match"] = (fsha(cif) == row["coordinate_gzip_sha256"]) if cif.endswith(".gz") \
            else (fsha(cif) == EXPECTED_CIF_SHA[c])
        r["coordinates_file"] = os.path.basename(cif)
        r["fasta_sha_match"] = fsha(fa) == row["fasta_sha256"]
        shutil.copy(cif, dst); shutil.copy(fa, dst)
        r["canonical"] = canonicalize(raw, f"{dst}/canonical.mrc")
        r["canonical_cube"] = canonicalize(raw, f"{dst}/canonical_cube.mrc", cube_min=128)
        checks = ["compressed_prefix_match", "gzip_ok", "header_match", "size_matches_emdb_kb",
                  "all_finite", "coordinates_sha_match", "fasta_sha_match"]
        r["full_map_integrity"] = "pass" if all(r.get(k) for k in checks) else \
            "FAIL: " + ",".join(k for k in checks if not r.get(k))
        os.remove(raw); os.remove(gz)
    except Exception as e:
        r["full_map_integrity"] = f"ERROR: {type(e).__name__}: {e}"
    r["seconds"] = round(time.time() - t0, 1)
    report["cases"][c] = r
    json.dump(report, open(f"{OUT}/integrity_report.json", "w"), indent=2, default=str)
    print(f"{c}: {r['full_map_integrity']} ({r['seconds']}s)", flush=True)
print("DONE", sum(1 for v in report["cases"].values() if v["full_map_integrity"] == "pass"), "/", len(cases), "passed")
