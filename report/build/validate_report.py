"""Validation of the built report against the project's quality-control checklist (checks A-I and rendering; 2026-10-13).
Run after rebuild_report.sh has built the Word file and the package (rebuild_report.sh calls it last).
  validate_report.py                      run all checks, write report/qa/QA_RECORD.md, exit 1 on any FAIL
  validate_report.py --record-visual "…"  record that the rendered pages of the CURRENT Word file were inspected (note = what was looked at)
Statuses: PASS / FAIL / NOT CHECKED. Nothing is marked PASS without the check having been executed in this run.
The numbers themselves are recomputed by sync_report_tables.py (run here as a subprocess); this script checks the
dataset, metric fields, alignment selection, the exported DOCX, the figures, the package and the regression cases."""
import csv, datetime, glob, hashlib, html, json, os, re, statistics as st, subprocess, sys, zipfile
R = "/Users/anony/Downloads/MICA_GAIVI_run"; B = f"{R}/benchmark"; QA = f"{R}/report/qa"; os.makedirs(QA, exist_ok=True)
DOCX = f"{R}/report/MICA_benchmark_report.docx"; ZIP = f"{R}/MICA_benchmark_reproducibility_package.zip"; P = f"{R}/reproducibility_package"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
now = datetime.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
if "--record-visual" in sys.argv:
    note = sys.argv[sys.argv.index("--record-visual") + 1]
    json.dump({"docx_sha256": sha(DOCX), "inspected": now, "note": note}, open(f"{QA}/visual_inspection.json", "w"), indent=1)
rounds = {r["case_id"]: r for r in csv.DictReader(open(f"{B}/case_rounds.csv"))}; cases = sorted(rounds); pid = lambda c: c.split("_")[1]
M = ["MICA", "EModelX(+AF)", "ModelAngelo", "CryoAtom", "EMProt"]
SETS = [os.path.basename(p)[:-5] for p in sorted(glob.glob(f"{R}/scoring/results/*.json"))]; res = {s: json.load(open(f"{R}/scoring/results/{s}.json")) for s in SETS}
met = lambda s, c: (res[s].get(c, {}).get("metrics") or {})
out = []   # (section, check, status, evidence)
def chk(sec, name, ok, ev): out.append((sec, name, "PASS" if ok else "FAIL", ev))
def skip(sec, name, ev): out.append((sec, name, "NOT CHECKED", ev))

# A. dataset
sel = {r["case_id"]: r for r in csv.DictReader(open(f"{B}/primary24_scope_screened_candidate.csv"))}
bands = {}
for c in cases: bands.setdefault(rounds[c]["resolution_stratum"], []).append(float(rounds[c]["resolution_A"]))
lim = {"1.5_to_lt2.5": (1.5, 2.5), "2.5_to_lt3.0": (2.5, 3.0), "3.0_to_lt3.5": (3.0, 3.5), "3.5_to_4.0": (3.5, 4.0001)}
chk("A", "24 unique cases; EMDB and PDB identifiers unique; 4 cases in each of 6 rounds", len(cases) == 24 and len({c.split("_")[0] for c in cases}) == 24 and len({pid(c) for c in cases}) == 24
    and sorted(sum(1 for c in cases if rounds[c]["round"] == str(r)) for r in range(1, 7)) == [4] * 6, "benchmark/case_rounds.csv")
chk("A", "6 cases per resolution band, every resolution inside its band", all(len(v) == 6 and all(lim[b][0] <= x < lim[b][1] for x in v) for b, v in bands.items()) and len(bands) == 4, "case_rounds.csv: " + str({b: sorted(v) for b, v in bands.items()}))
dates = [(pid(c), sel[c]["pdb_initial_release_date"][:10], sel[c]["emdb_initial_map_release_date"][:10]) for c in cases if c in sel]
chk("A", "PDB and EMDB initial release dates both in 2026 for all 24 cases (as recorded at screening)", len(dates) == 24 and all(a.startswith("2026") and b.startswith("2026") for _, a, b in dates),
    "benchmark/primary24_scope_screened_candidate.csv (recorded 2026-10-04 from RCSB PDB and EMDB; not re-queried in this run)")
refs = {c: {met(s, c).get("ref_ca") for s in SETS if met(s, c).get("ref_ca") is not None} for c in cases}
chk("A", "polymer reference Cα count identical in every result set, per case", all(len(v) == 1 for v in refs.values()), "scoring/results/*.json ref_ca; sets: " + ", ".join(SETS))

# B. metrics
six = ["tm_score", "aligned_ca", "seq_identity", "ca_match_pct", "ca_quality", "seq_match_pct"]
bad = [(m, pid(c), k) for m in M for c in cases for k in six if not isinstance(met(m, c).get(k), (int, float))]
chk("B", "all six declared metrics present for 120 five-method runs", not bad, "missing: " + str(bad[:5]) if bad else "tm_score, aligned_ca, seq_identity, ca_match_pct, ca_quality, seq_match_pct")
rng = [(m, pid(c), k) for m in M for c in cases for k, hi in (("tm_score", 1), ("seq_identity", 1), ("ca_quality", 1), ("ca_match_pct", 100), ("seq_match_pct", 100)) if not 0 <= met(m, c)[k] <= hi * 1.0001]
rng += [(m, pid(c), "aligned_ca") for m in M for c in cases if not 0 < met(m, c)["aligned_ca"] <= min(met(m, c)["ref_ca"], met(m, c)["pred_ca"])]
chk("B", "metric values inside their defined ranges (fractions 0-1, percentages 0-100, aligned length <= smaller model)", not rng, str(rng[:5]) if rng else "120 runs x 6 metrics")
rt = open(f"{R}/report/RESULTS_TABLES.md").read()
chk("B", "each of the six metrics has its own appendix table", all(t in rt for t in ("Table 1. Complex TM-score", "Table 13. Aligned Cα length", "Table 14. Sequence identity", "Table 15. Cα match", "Table 2. Cα quality", "Table 3. Sequence match")), "report/RESULTS_TABLES.md Tables 1, 13, 14, 15, 2, 3")

# C-F. aggregates, intervals, populations, completion accounting: recomputed by sync_report_tables.py
r = subprocess.run([sys.executable, f"{R}/report/build/sync_report_tables.py", "--registry"], capture_output=True, text=True); last = (r.stdout.strip().split("\n") or [""])[-1]
chk("C", "means, medians, win counts, subgroup means and sensitivity values in the text recomputed from full-precision results", r.returncode == 0, "sync_report_tables.py: " + last[:160])
ties = []
for c in cases:
    v = sorted(((met(m, c)["tm_score"], m) for m in M), reverse=True)
    if f"{v[0][0]:.3f}" == f"{v[1][0]:.3f}": ties.append(f"{pid(c)} {v[0][1]} {v[0][0]:.5f} / {v[1][1]} {v[1][0]:.5f}")
draft = open(f"{R}/report/REPORT_DRAFT.md").read()
chk("C", "tie handling documented (best value decided on unrounded scores)", "unrounded" in draft and "unrounded" in rt, f"displayed ties for best TM among the five methods: {ties or 'none'}; rule stated in section 6.1 and in the appendix header")
norm = lambda s: s.replace("−", "-")
iv = re.findall(r"[+-]\d\.\d{3} \([+-]\d\.\d{3} to [+-]\d\.\d{3}\)", norm(draft[draft.index("### 6.2"):draft.index("### 6.4")])); t10 = norm(rt[rt.index("### Table 10."):rt.index("### Table 13.")])
def in_app(x):   # the appendix shows the mean and the interval in separate cells of one row
    mean_, rng_ = re.fullmatch(r"([+-]\d\.\d{3}) \((.*)\)", x).groups(); return any(mean_ in l and rng_ in l for l in t10.split("\n"))
chk("D", "confidence intervals of sections 6.2-6.3 identical to appendix Tables 10 and 12 (two code paths, one seeded procedure)", len(iv) >= 11 and all(in_app(x) for x in iv), f"{len(iv)} intervals; cases resampled, 10,000 resamples, percentile interval, random.Random(20261004) per comparison")
chk("E", "subset comparisons use one population for every method (CryoZeta 14-case and 20-case figures)", r.returncode == 0, "computed in sync_report_tables.py from the same case list for all methods")
acc = [l for l in r.stdout.split("\n") if l.startswith("CryoZeta outcomes")]
n = [int(x) for x in re.findall(r"(\d+) (?:complete|without|not attempted|no record)", acc[0])] if acc else []
chk("F", "CryoZeta statuses are mutually exclusive and sum to 24", len(n) == 4 and sum(n) == 24, acc[0] if acc else "no accounting line")
part = set(res.get("CryoZeta_partial", {})) & {c for c in cases if met("CryoZeta", c).get("tm_score") is not None}
chk("F", "no case counted both as complete and as partial", not part, str(sorted(part)) or "disjoint")

# G. the exported DOCX
z = zipfile.ZipFile(DOCX); x = z.read("word/document.xml").decode()
def docx_tables():
    T = []
    for tb in re.findall(r"<w:tbl>.*?</w:tbl>", x, flags=re.S):
        T.append([[html.unescape("".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", tc))).strip() for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, flags=re.S)] for tr in re.findall(r"<w:tr[ >].*?</w:tr>", tb, flags=re.S)])
    return T
def md_tables(text):
    T = []; cur = None
    for l in text.split("\n"):
        if l.startswith("|"):
            row = [re.sub(r"\s+", " ", re.sub(r"\*\*|`", "", c)).strip() for c in l.strip().strip("|").split("|")]
            if cur is None: cur = []
            if not all(re.fullmatch(r":?-+:?", c) for c in row): cur.append(row)
        elif cur is not None: T.append(cur); cur = None
    if cur is not None: T.append(cur)
    return T
src = md_tables(draft) + md_tables(rt); dt = [[[re.sub(r"\s+", " ", c) for c in row] for row in t] for t in docx_tables()]
diff = []
for i, (a, b) in enumerate(zip(src, dt)):
    if a != b:
        for ra, rb in zip(a, b):
            if ra != rb: diff.append(f"table {i + 1} row {ra[0]}: source {ra} / docx {rb}"); break
        else: diff.append(f"table {i + 1}: row count {len(a)} / {len(b)}")
chk("G", "every table of the Word file equals its generated source, cell for cell", len(src) == len(dt) and not diff, f"{len(dt)} tables, {sum(len(t) for t in dt)} rows" if not diff and len(src) == len(dt) else f"{len(src)} source / {len(dt)} docx tables; " + "; ".join(diff[:3]))
media = sorted(n for n in z.namelist() if n.startswith("word/media/") and not n.endswith("/")); figs = sorted(glob.glob(f"{R}/report/figures/image[0-9].png")); fh = {sha(f) for f in figs}
chk("G", "figures in the Word file are the files regenerated from the current results in this build", len(media) == len(figs) and {hashlib.sha256(z.read(n)).hexdigest() for n in media} == fh
    and all(os.path.getmtime(f) >= max(os.path.getmtime(f"{R}/scoring/results/{m}.json") for m in M) for f in figs), f"{len(media)} embedded images = report/figures/image1.png, image2.png (make_figures.py reads scoring/results)")
chk("G", "no tracked changes or comments in the Word file", not re.search(r"<w:(ins|del) ", x) and "word/comments.xml" not in z.namelist(), "word/document.xml")
figrefs = set(re.findall(r"Figure (\d)", draft)); tabrefs = {int(t) for t in re.findall(r"Tables? (\d+)", draft)} | {int(t) for a, b in re.findall(r"Tables (\d+)[–-](\d+)", draft) for t in (a, b)}
have = {int(t) for t in re.findall(r"### Table (\d+)", rt)}
chk("G", "every table and figure cited in the text exists", tabrefs <= have and len(figrefs) <= len(figs), f"tables cited {sorted(tabrefs)}; figures cited {sorted(figrefs)}")

# H. alignment selection
A = {(a["set"], a["case"]): a for a in json.load(open(f"{R}/scoring/usalign_fixed_build_audit.json"))}; hb = []; nH = 0
for s in SETS:
    for c, v in res[s].items():
        m = v.get("metrics") or {}
        if m.get("tm_score") is None: continue
        a = A.get((s, c)); nH += 1
        if not a: hb.append((s, pid(c), "not in audit")); continue
        src_ = m.get("tm_source"); e = a.get(src_, {})
        if not (m["tm_score"] == e.get("tm") == max(a["search"]["tm"], a["mapframe"]["tm"]) and m["aligned_ca"] == e.get("aligned") and m["seq_identity"] == e.get("seq_id") and "746ad8a" in m.get("usalign_build", "")): hb.append((s, pid(c), src_))
chk("H", "TM-score, aligned length and sequence identity come from the same evaluation; same candidate rule and build for every model", not hb, f"{nH} models: primary = higher of fixed-build search and map-frame, US-align 20261012 (746ad8a)" if not hb else str(hb[:5]))

# I. package
man = json.load(open(f"{P}/models/MODEL_MANIFEST.json")); mb = [m["file"] for m in man if not os.path.exists(f"{P}/{m['file']}") or sha(f"{P}/{m['file']}") != m["sha256"] or not m["sha256_matches_scored_model"]]
chk("I", "every model file in the package matches its manifest hash and the hash stored with its score", not mb, f"{len(man)} manifest rows" if not mb else str(mb[:5]))
zf = zipfile.ZipFile(ZIP); zn = {n for n in zf.namelist() if not n.endswith("/")}; dn = {os.path.relpath(os.path.join(d, f), R) for d, _, fs in os.walk(P) for f in fs if f != ".DS_Store"}
chk("I", "zip archive contains exactly the package folder", zn == dn, f"{len(zn)} files; zip sha256 {sha(ZIP)[:16]}…" if zn == dn else f"only in zip {sorted(zn - dn)[:3]}; only in folder {sorted(dn - zn)[:3]}")
chk("I", "report in the package is this Word file; its results are the current result files", hashlib.sha256(zf.read("reproducibility_package/report/MICA_benchmark_report.docx")).hexdigest() == sha(DOCX)
    and all(open(f"{R}/scoring/results/{s}.json", "rb").read() == zf.read("reproducibility_package/scoring/results/" + {"CryoZeta": "CryoZeta_interim", "CryoZeta_partial": "CryoZeta_partial_models"}.get(s, s) + ".json") for s in SETS), "sha256 comparison inside the zip")
chk("I", "restricted and external inputs are listed as excluded (maps, Phenix, AlphaFold Server outputs)", all(k in open(f"{P}/README.md").read() for k in ("Cryo-EM maps", "Phenix", "AlphaFold Server")) and not [n for n in zn if re.search(r"phenix-installer|\.mrc$|\.map(\.gz)?$", n)], "reproducibility_package/README.md, section 'Not included, and why'")

# regression cases named in the protocol
chk("REG", "9P9D: polymer-only reference count in every result set (no ligand residues counted)", len(refs[[c for c in cases if pid(c) == "9P9D"][0]]) == 1, f"ref_ca = {refs[[c for c in cases if pid(c) == '9P9D'][0]]}")
c22 = [c for c in cases if pid(c) == "22HP"][0]; a22 = A[("CryoAtom", c22)]
chk("REG", "22HP / CryoAtom: fixed-build search and map-frame evaluations agree (the old build returned 0.277-0.733)", abs(a22["search"]["tm"] - a22["mapframe"]["tm"]) < 0.01 and met("CryoAtom", c22)["tm_score"] == max(a22["search"]["tm"], a22["mapframe"]["tm"]),
    f"search {a22['search']['tm']}, map-frame {a22['mapframe']['tm']}, sequence identity {met('CryoAtom', c22)['seq_identity']}")
cell = lambda t, case, col: next(row[col] for row in t if row[0] == case)
t61 = next(t for t in dt if t[0][:5] == ["Case", "Round", "Res. (Å)", "Ref. Cα", "MICA"]); tA = [t for t in dt if t[0][:5] == ["Case", "Round", "Res. (Å)", "Ref. Cα", "MICA"]][1]
want = {("9U62", 7): f"{met('CryoAtom', [c for c in cases if pid(c) == '9U62'][0])['tm_score']:.3f}", ("9EC7", 8): f"{met('EMProt', [c for c in cases if pid(c) == '9EC7'][0])['tm_score']:.3f}"}
chk("REG", "9U62/CryoAtom and 9EC7/EMProt: section 6.1 cell = appendix Table 1 cell = result file", all(cell(t61, k, j) == cell(tA, k, j) == v for (k, j), v in want.items()), str({k: (cell(t61, k, j), cell(tA, k, j), v) for (k, j), v in want.items()}))

# section 9: structure and rendering
sizes = [int(s) for tb in re.findall(r"<w:tbl>.*?</w:tbl>", x, flags=re.S) for s in re.findall(r'<w:sz w:val="(\d+)"', tb)]; _b = re.findall(r'<w:sz w:val="(\d+)"', re.sub(r"<w:tbl>.*?</w:tbl>", "", x, flags=re.S)); body = re.match(r"(\d+)", max(set(_b), key=_b.count))   # most frequent run size outside tables
ntb = len(re.findall(r"<w:tbl>", x)); nhead = len(re.findall(r"<w:tblHeader/>", x)); nbody = len(re.findall(r'<w:tblHeader w:val="off"/>', x))
chk("9", "table header row (and only it) repeats on every page; rows do not split across pages", nhead == ntb and nbody == sum(len(t) for t in dt) - ntb and len(re.findall(r"<w:cantSplit", x)) >= sum(len(t) for t in dt), f"{ntb} tables, {nhead} repeating header rows, {nbody} body rows")
chk("9", "font sizes: body about 11 pt, tables 8.5 pt or larger", body and int(body.group(1)) >= 21 and min(sizes) >= 17, f"body {int(body.group(1)) / 2} pt; smallest table text {min(sizes) / 2} pt")
pdf = f"{QA}/render/MICA_benchmark_report.pdf"; pages = None
if os.path.exists(pdf) and os.path.getmtime(pdf) >= os.path.getmtime(DOCX):
    pages = int(re.search(r"Pages:\s+(\d+)", subprocess.run(["pdfinfo", pdf], capture_output=True, text=True).stdout).group(1))
    txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True).stdout
    chk("9", "Word file renders to PDF; no replacement glyphs; no leftover markdown or draft markers", pages > 0 and "�" not in txt and not re.search("|".join(l.strip() for l in open(f"{QA}/leftover_markers.txt") if l.strip()), txt), f"{pages} pages (LibreOffice), report/qa/render/")
else: skip("9", "Word file renders to PDF", "no render newer than the Word file (run rebuild_report.sh)")
vi = json.load(open(f"{QA}/visual_inspection.json")) if os.path.exists(f"{QA}/visual_inspection.json") else {}
if vi.get("docx_sha256") == sha(DOCX): out.append(("9", "page images inspected for clipping, overlap, page breaks, table continuation, near-empty pages", "PASS", f"{vi['inspected']}: {vi['note']}"))
else: skip("9", "page images inspected for clipping, overlap, page breaks, table continuation, near-empty pages", "not yet inspected for this exact file")
skip("6", "references re-verified against publishers / repositories", "last verified 2026-10-13 during review fixes; not re-queried in this run")
skip("1", "baseline document 'MICA_benchmark_report(4).docx'", "no file of that name on this computer; the report is built from report/REPORT_DRAFT.md + scoring/results")

L = ["# QA record — MICA benchmark report", "", f"Generated {now} by `report/build/validate_report.py` (checklist: `REPORT_QC_PROTOCOL.md`).", "",
     f"- Report: `report/MICA_benchmark_report.docx`, sha256 `{sha(DOCX)}`" + (f", {pages} pages" if pages else ""),
     f"- Package: `MICA_benchmark_reproducibility_package.zip`, sha256 `{sha(ZIP)}`, {len(zn)} files",
     "- Data snapshot: `scoring/results/*.json` (" + ", ".join(f"{s} {sum(1 for v in res[s].values() if (v.get('metrics') or {}).get('tm_score') is not None)}" for s in SETS) + " scored models)",
     "- Evaluator: `scoring/evaluate.py` with US-align Version 20261012 (commit 746ad8a); cctbx chain_comparison", "",
     "| Checklist section | Check | Status | Evidence |", "|---|---|---|---|"] + [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in out]
cnt = {s: sum(1 for o in out if o[2] == s) for s in ("PASS", "FAIL", "NOT CHECKED")}
L += ["", f"Totals: {cnt['PASS']} PASS, {cnt['FAIL']} FAIL, {cnt['NOT CHECKED']} NOT CHECKED.", "", "Open and accepted items: `report/qa/ISSUE_LEDGER.md`. Claim-value registry: `report/qa/claim_registry.json`.", ""]
open(f"{QA}/QA_RECORD.md", "w").write("\n".join(L))
for a, b, c, d in out:
    if c != "PASS": print(f"{c}: [{a}] {b} — {d[:300]}")
print(f"QA: {cnt['PASS']} PASS, {cnt['FAIL']} FAIL, {cnt['NOT CHECKED']} NOT CHECKED -> report/qa/QA_RECORD.md")
sys.exit(1 if cnt["FAIL"] else 0)
