"""Generate the number-bearing tables of report section 6 (and a list of facts for the prose) from scoring/results.
Same primary-set rule as scoring/make_tables.py (complete rounds; final failures from scoring/final_failures.json).
Writes report/section6_generated.md; prints facts. No number in the report is typed from memory."""
import csv, json, os, random, statistics as st
R = "/Users/anony/Downloads/MICA_GAIVI_run"; B = f"{R}/benchmark"
rounds = {r["case_id"]: r for r in csv.DictReader(open(f"{B}/case_rounds.csv"))}
M = ["MICA", "EModelX(+AF)", "ModelAngelo", "CryoAtom", "EMProt"]
res = {m: json.load(open(f"{R}/scoring/results/{m}.json")) for m in M + ["MICA_unrefined"]}
FF = json.load(open(f"{R}/scoring/final_failures.json")) if os.path.exists(f"{R}/scoring/final_failures.json") else {}
low = json.load(open(f"{B}/audit3_low_homology_by_case.json")); low.update(json.load(open(f"{B}/audit3_low_homology_by_case_ext.json")))
g = lambda m, c, k: (res[m].get(c, {}).get("metrics") or {}).get(k)
failed = lambda m, c: c in FF.get(m, {}) and g(m, c, "tm_score") is None
done = lambda m, c: g(m, c, "tm_score") is not None or failed(m, c)
by_round = {}
for c in sorted(rounds, key=lambda c: (int(rounds[c]["round"]), float(rounds[c]["resolution_A"]))): by_round.setdefault(int(rounds[c]["round"]), []).append(c)
complete = [r for r, cs in by_round.items() if all(done(m, c) for m in M for c in cs)]
cases = [c for r in complete for c in by_round[r]]; N = len(cases); pid = lambda c: c.split("_")[1]
V = lambda m, k, cs=None: [g(m, c, k) for c in (cases if cs is None else cs) if g(m, c, k) is not None]
def wins(cs): return {m: sum(1 for c in cs if g(m, c, "tm_score") is not None and g(m, c, "tm_score") == max(g(x, c, "tm_score") for x in M if g(x, c, "tm_score") is not None)) for m in M}
def boot(d, n=10000):
    rng = random.Random(20261004); ms = sorted(st.mean(rng.choices(d, k=len(d))) for _ in range(n)); return ms[int(0.025 * n)], ms[int(0.975 * n) - 1]
def paired(m, k, cs=None): return [g("MICA", c, k) - g(m, c, k) for c in (cases if cs is None else cs) if g("MICA", c, k) is not None and g(m, c, k) is not None]
w = wins(cases); o = []
o += ["| Method | Cases scored | Mean TM | Median TM | Mean Cα quality | Mean sequence match (%) | Mean map-frame precision / recall @3 Å | Best TM in |", "|---|---|---|---|---|---|---|---|"]
for m in M:
    o.append(f"| {m} | {len(V(m, 'tm_score'))}/{N} | {st.mean(V(m, 'tm_score')):.3f} | {st.median(V(m, 'tm_score')):.3f} | {st.mean(V(m, 'ca_quality')):.3f} | {st.mean(V(m, 'seq_match_pct')):.1f} | "
             f"{st.mean(V(m, 'mapframe_precision_3A')):.3f} / {st.mean(V(m, 'mapframe_recall_3A')):.3f} | {w[m]} |")
o += ["", "PER-CASE", "", "| Case | Round | Res. (Å) | Ref. Cα | " + " | ".join(M) + " |", "|---|---|---|---|" + "---|" * len(M)]
for c in cases:
    v = [g(m, c, "tm_score") for m in M]; b = max(x for x in v if x is not None); ref = next(g(m, c, "ref_ca") for m in M if g(m, c, "ref_ca"))
    o.append(f"| {pid(c)} | {rounds[c]['round']} | {float(rounds[c]['resolution_A']):.2f} | {ref} | " + " | ".join(("fail" if failed(m, c) else "–") if x is None else (f"**{x:.3f}**" if x == b else f"{x:.3f}") for m, x in zip(M, v)) + " |")
o += ["", "PAIRED", "", "| MICA − … | TM-score: mean difference (95% interval) | Cα quality: mean difference (95% interval) |", "|---|---|---|"]
for m in M[1:]:
    cell = []
    for k in ("tm_score", "ca_quality"):
        d = paired(m, k); lo, hi = boot(d); cell.append(f"{st.mean(d):+.3f} ({lo:+.3f} to {hi:+.3f}); MICA higher in {sum(x > 0 for x in d)}/{len(d)}")
    o.append(f"| {m} | {cell[0]} | {cell[1]} |")
open(f"{R}/report/section6_generated.md", "w").write("\n".join(o) + "\n")
print("\n".join(o)); print("\nFACTS"); print("complete rounds", complete, "N", N, "runs with model", sum(len(V(m, "tm_score")) for m in M), "of", N * len(M))
common = [c for c in cases if all(g(m, c, "tm_score") is not None for m in M)]
print("mean TM on cases scored for all methods (n=%d):" % len(common), {m: round(st.mean(g(m, c, "tm_score") for c in common), 3) for m in M})
d = sorted(((g("MICA", c, "tm_score") - g("EModelX(+AF)", c, "tm_score"), pid(c)) for c in cases if g("EModelX(+AF)", c, "tm_score") is not None))
print("MICA-EMX extremes:", [(p, round(x, 3)) for x, p in d[:3]], [(p, round(x, 3)) for x, p in d[-3:]])
for name, st_ in (("core", "core"), ("extension", "extension")):
    cs = [c for c in cases if rounds[c]["stage"] == st_]
    if not cs: continue
    print(name, len(cs), "mean TM", {m: round(st.mean(V(m, "tm_score", cs)), 3) for m in M}, "wins", wins(cs))
    for m in M[1:]:
        dd = paired(m, "tm_score", cs); lo, hi = boot(dd); print("   MICA -", m, f"{st.mean(dd):+.3f} ({lo:+.3f} to {hi:+.3f}) {sum(x > 0 for x in dd)}/{len(dd)}")
bands = {}
for c in cases: bands.setdefault(rounds[c]["resolution_stratum"], []).append(c)
for b in sorted(bands): print("band", b, len(bands[b]), {m: (round(st.mean(V(m, "tm_score", bands[b])), 3), len(V(m, "tm_score", bands[b]))) for m in M})
sub = [c for c in cases if low[c]["frac_residues_low_pre2023_05"] >= 0.5]; rest = [c for c in cases if c not in sub]
print("low-homology", [pid(c) for c in sub], {m: (round(st.mean(V(m, "tm_score", sub)), 3), len(V(m, "tm_score", sub))) for m in M})
print("other", len(rest), {m: round(st.mean(V(m, "tm_score", rest)), 3) for m in M})
for m in M[1:]:
    dd = paired(m, "tm_score", sub)
    if len(dd) > 1: lo, hi = boot(dd); print("   low-hom MICA -", m, f"{st.mean(dd):+.3f} ({lo:+.3f} to {hi:+.3f}) {sum(x > 0 for x in dd)}/{len(dd)}")
for m in M:
    r = [(g(m, c, "pred_ca") / g(m, c, "ref_ca"), pid(c)) for c in cases if g(m, c, "pred_ca")]
    print("size ratio", m, "min %.2f %s max %.2f %s" % (*min(r), *max(r)), "| >1.2x:", [(p, round(x, 2)) for x, p in r if x > 1.2])
dr = [(abs(g("MICA", c, "tm_score") - g("MICA_unrefined", c, "tm_score")), pid(c)) for c in cases if g("MICA_unrefined", c, "tm_score") is not None]
print("refinement max |dTM| %.4f (%s), n=%d" % (*max(dr), len(dr)))
for m in M:
    fc = [res[m][c].get("frame_check_frac_ca_within_3A") for c in cases if c in res[m] and res[m][c].get("frame_check_frac_ca_within_3A") is not None]
    print("frame check", m, min(fc), max(fc), "| <0.85:", [(pid(c), res[m][c]["frame_check_frac_ca_within_3A"]) for c in cases if c in res[m] and (res[m][c].get("frame_check_frac_ca_within_3A") or 1) < 0.85])
print("lowest precision cases:", {m: sorted(((round(g(m, c, "mapframe_precision_3A"), 3), pid(c)) for c in cases if g(m, c, "mapframe_precision_3A") is not None))[:2] for m in M})
print("ca quality best counts:", {m: sum(1 for c in cases if g(m, c, "ca_quality") is not None and g(m, c, "ca_quality") == max(g(x, c, "ca_quality") for x in M if g(x, c, "ca_quality") is not None)) for m in M})
