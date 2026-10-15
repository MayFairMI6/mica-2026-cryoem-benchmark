"""Keep the number tables in report/REPORT_DRAFT.md identical to scoring/results (2026-10-13).
Every table of section 6 (6.1 summary and per-case, 6.2 paired, 6.3 core/extension, 6.8 CryoZeta complete / means / partial)
is recomputed here with the same rules as scoring/make_tables.py and compared cell by cell with the text.
  sync_report_tables.py          check only; prints every difference; exit 1 if any (the build stops)
  sync_report_tables.py --fix    rewrite the differing rows in REPORT_DRAFT.md, then check again
Also checked: the per-case TM table of the text against Table 1 of report/RESULTS_TABLES.md (the appendix), and the
sentences that repeat table numbers. No number is typed by hand."""
import csv, glob, json, os, random, re, statistics as st, sys
R = "/Users/anony/Downloads/MICA_GAIVI_run"; B = f"{R}/benchmark"; DRAFT = f"{R}/report/REPORT_DRAFT.md"
rounds = {r["case_id"]: r for r in csv.DictReader(open(f"{B}/case_rounds.csv"))}
M = ["MICA", "EModelX(+AF)", "ModelAngelo", "CryoAtom", "EMProt"]
load = lambda n: json.load(open(f"{R}/scoring/results/{n}.json")) if os.path.exists(f"{R}/scoring/results/{n}.json") else {}
res = {m: load(m) for m in M + ["CryoZeta", "CryoZeta_partial"]}
g = lambda m, c, k: (res[m].get(c, {}).get("metrics") or {}).get(k)
cases = sorted(rounds, key=lambda c: (int(rounds[c]["round"]), float(rounds[c]["resolution_A"]))); N = len(cases); pid = lambda c: c.split("_")[1]
assert all(g(m, c, "tm_score") is not None for m in M for c in cases), "a five-method result is missing"
V = lambda m, k, cs=cases: [g(m, c, k) for c in cs if g(m, c, k) is not None]
f3 = lambda x: f"{x:.3f}"; sg = lambda x: f"{x:+.3f}".replace("-", "−")
def best(c, ms): return max(g(m, c, "tm_score") for m in ms if g(m, c, "tm_score") is not None)
def wins(cs): return {m: sum(1 for c in cs if g(m, c, "tm_score") == best(c, M)) for m in M}
def boot(d, n=10000):
    rng = random.Random(20261004); ms = sorted(st.mean(rng.choices(d, k=len(d))) for _ in range(n)); return ms[int(0.025 * n)], ms[int(0.975 * n) - 1]
def paired(m, k, cs=cases): return [g("MICA", c, k) - g(m, c, k) for c in cs]
def diffcell(d, long): lo, hi = boot(d); return f"{sg(st.mean(d))} ({sg(lo)} to {sg(hi)}); {'MICA higher in ' if long else ''}{sum(x > 0 for x in d)}/{len(d)}"
def bold(x, b): return f"**{f3(x)}**" if x == b else f3(x)
def rec(kind, c):
    p = f"{B}/gaivi_bench/out/{kind}/{c}/run_record.json"; return json.load(open(p)) if os.path.exists(p) else {}

T = {}   # name -> (text that identifies the header row, expected data rows as lists of cells)
w = wins(cases)
T["6.1 summary"] = ("| Method | Mean TM |", [[m, f3(st.mean(V(m, "tm_score"))), f3(st.median(V(m, "tm_score"))), f3(st.mean(V(m, "ca_quality"))), f"{st.mean(V(m, 'seq_match_pct')):.1f}",
    f"{f3(st.mean(V(m, 'mapframe_precision_3A')))} / {f3(st.mean(V(m, 'mapframe_recall_3A')))}", "0" if w[m] == 0 else f"{w[m]} case{'s' if w[m] > 1 else ''}"] for m in M])
T["6.1 per case"] = ("| Case | Round | Res. (Å) | Ref. Cα | MICA |", [[pid(c), rounds[c]["round"], f"{float(rounds[c]['resolution_A']):.2f}", str(g("MICA", c, "ref_ca"))] +
    [bold(g(m, c, "tm_score"), best(c, M)) for m in M] for c in cases])
T["6.2 paired"] = ("| MICA − … |", [[m, diffcell(paired(m, "tm_score"), m == M[1]), diffcell(paired(m, "ca_quality"), False)] for m in M[1:]])
rows = []
for name, cs in (("Core (rounds 1–3)", [c for c in cases if rounds[c]["stage"] == "core"]), ("Extension (rounds 4–6)", [c for c in cases if rounds[c]["stage"] == "extension"]), ("All", cases)):
    d = paired("EModelX(+AF)", "tm_score", cs); lo, hi = boot(d); ww = wins(cs)
    rows.append([name, str(len(cs))] + [f3(st.mean(V(m, "tm_score", cs))) for m in M] + [f"{sg(st.mean(d))} ({sg(lo)} to {sg(hi)})", f"{ww['MICA']} / {ww['EModelX(+AF)']} / {sum(ww[m] for m in M[2:])}"])
T["6.3 core/extension"] = ("| Set | n | MICA |", rows)
cz = [c for c in cases if g("CryoZeta", c, "tm_score") is not None]; S = ["CryoZeta", "MICA", "EModelX(+AF)"]; rows = []
for c in cz:
    r = rec("CryoZeta", c); mode = r.get("mode_used", "?") + (" (after normal mode ran out of GPU memory)" if r.get("normal_mode_oom") else "")
    rows.append([pid(c), f"{float(rounds[c]['resolution_A']):.2f}", str(g("CryoZeta", c, "ref_ca")), mode, f"{g('CryoZeta', c, 'pred_ca') / g('CryoZeta', c, 'ref_ca'):.2f}×"] + [bold(g(m, c, "tm_score"), best(c, S)) for m in S])
T["6.8 CryoZeta complete"] = ("| Case | Res. (Å) | Ref. Cα | CryoZeta mode |", rows)
A = ["CryoZeta"] + M
T["6.8 CryoZeta means"] = ("| Mean over the ", [[lab] + [f3(st.mean(V(m, k, cz))) for m in A] for lab, k in (("TM-score", "tm_score"), ("Cα quality", "ca_quality"),
    ("Map-frame precision @3 Å", "mapframe_precision_3A"), ("Map-frame recall @3 Å", "mapframe_recall_3A"))])
part = [c for c in cases if c in res["CryoZeta_partial"] and g("CryoZeta_partial", c, "tm_score") is not None and c not in cz]; rows = []
for c in part:
    r = rec("CryoZeta_partial", c); P = "CryoZeta_partial"
    rows.append([pid(c), f"{float(rounds[c]['resolution_A']):.2f}", str(g(P, c, "ref_ca")), f"{r.get('stages_completed', '?')} of {r.get('stages_total', '?')}", f"{g(P, c, 'pred_ca') / g(P, c, 'ref_ca'):.2f}×",
                 f3(g(P, c, "tm_score")), f3(g(P, c, "ca_quality")), f3(g("MICA", c, "tm_score")), f3(g("EModelX(+AF)", c, "tm_score"))])
T["6.8 CryoZeta partial"] = ("| Case | Res. (Å) | Ref. Cα | Stages completed |", rows)

cells = lambda line: [x.strip() for x in line.strip().strip("|").split("|")]
def run(fix):
    L = open(DRAFT).read().split("\n"); problems = []
    for name, (head, exp) in T.items():
        hits = [i for i, l in enumerate(L) if l.startswith(head)]
        if len(hits) != 1: problems.append(f"{name}: header found {len(hits)} times"); continue
        i = hits[0]
        if name == "6.8 CryoZeta means":
            want = f"| Mean over the {len(cz)} cases | " + " | ".join(A) + " |"
            if L[i] != want:
                problems.append(f"{name}: header says '{cells(L[i])[0]}', results have {len(cz)} cases")
                if fix: L[i] = want
        j = i + 2; k = j
        while k < len(L) and L[k].startswith("|"): k += 1
        have = [cells(l) for l in L[j:k]]
        if have != exp:
            for a in range(max(len(have), len(exp))):
                h = have[a] if a < len(have) else None; e = exp[a] if a < len(exp) else None
                if h != e:
                    if h is None or e is None or len(h) != len(e) or h[0] != e[0]: problems.append(f"{name}: row {a + 1} text {h and h[0]} / results {e and e[0]} (row added, removed or reordered)")
                    else: problems += [f"{name}: {e[0]}, column {b + 1}: text {x} / results {y}" for b, (x, y) in enumerate(zip(h, e)) if x != y]
            if fix: L[j:k] = ["| " + " | ".join(r) + " |" for r in exp]
    if fix: open(DRAFT, "w").write("\n".join(L))
    return problems
p = run("--fix" in sys.argv)
if "--fix" in sys.argv:
    print(f"fixed {len(p)} difference(s):" if p else "tables already identical"); print("\n".join("  " + x for x in p)); p = run(False)

# appendix Table 1 (report/RESULTS_TABLES.md) against the same results
rt = open(f"{R}/report/RESULTS_TABLES.md").read().split("\n"); i = next(i for i, l in enumerate(rt) if l.startswith("| Case |") and "MICA" in l); t1 = {}
for l in rt[i + 2:]:
    if not l.startswith("|"): break
    c = cells(l); t1[c[0]] = [x.replace("*", "") for x in c[-5:]]
for c in cases:
    e = [f3(g(m, c, "tm_score")) for m in M]
    if t1.get(pid(c)) != e: p.append(f"appendix Table 1: {pid(c)} has {t1.get(pid(c))}, results {e} (run scoring/make_tables.py)")

# sentences that repeat table numbers
import datetime as _dt
_sel = {x["case_id"]: x for x in csv.DictReader(open(f"{B}/primary24_scope_screened_candidate.csv"))}
_rel = sorted(min(_sel[c]["pdb_initial_release_date"][:10], _sel[c]["emdb_initial_map_release_date"][:10]) for c in cases); _hi = max(max(_sel[c]["pdb_initial_release_date"][:10], _sel[c]["emdb_initial_map_release_date"][:10]) for c in cases)
_d = lambda x: f"{_dt.date.fromisoformat(x).day} {_dt.date.fromisoformat(x).strftime('%B %Y')}"
assert _rel[0] > "2025-06-30", "a benchmark entry was released inside CryoZeta's training period"
names = lambda cs: ", ".join(pid(c) for c in cs)
txt = open(DRAFT).read(); mean = lambda m, cs=cases: f3(st.mean(V(m, "tm_score", cs)))
facts = [f"first released between {_d(_rel[0])} and {_d(_hi)}", f"highest mean complex TM-score ({mean('MICA')})", f"{mean('MICA')} for MICA and {mean('EModelX(+AF)')} for EModelX(+AF)"]
if cz:
    hi = sum(g("CryoZeta", c, "tm_score") > g("MICA", c, "tm_score") for c in cz); top = sum(g("CryoZeta", c, "tm_score") == best(c, A) for c in cz)
    facts += [f"({mean('CryoZeta', cz)} against {mean('MICA', cz)})", f"higher than MICA in {hi} of {len(cz)} cases", f"highest of the six methods in {top}", f"{len(cz)} cases had produced a complete model", f"Of the {len(cz)} complete models"]
# derived statements: leave-one-out on rounded Table 1 values (section 6.9), refinement effect, gap before the evaluator correction
rd = {c: round(g("MICA", c, "tm_score"), 3) - round(g("EModelX(+AF)", c, "tm_score"), 3) for c in cases}; loo = lambda drop: f"+{st.mean(v for c, v in rd.items() if pid(c) not in drop):.4f}"
facts += [f"mean difference is {loo([])}", f"changes it to {loo(['9EFH'])} (n = 23)", f"changes it to {loo(['24VT'])} (n = 23)", f"Omitting both yields {loo(['9EFH', '24VT'])} (n = 22)"]
unref = load("MICA_unrefined"); dr = max(abs(g("MICA", c, "tm_score") - unref[c]["metrics"]["tm_score"]) for c in cases if (unref.get(c, {}).get("metrics") or {}).get("tm_score") is not None)
facts.append(f"by no more than {dr:.4f} on any case")
old = lambda m: st.mean(res[m][c]["metrics"]["tm_score_recorded"] for c in cases)
facts.append(f"gap is {st.mean(V('MICA', 'tm_score')) - st.mean(V('CryoAtom', 'tm_score')):.3f} instead of {old('MICA') - old('CryoAtom'):.3f}")
# CryoZeta outcome counts of section 6.8, from the run records: complete / attempted without a final model / not attempted / still running
if cz:
    W = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen"]
    rr = {c: rec("CryoZeta", c) for c in cases}; nomodel = [c for c in cases if rr[c] and rr[c].get("exit") != -1 and c not in cz]; skipped = [c for c in cases if rr[c].get("exit") == -1]
    running = [c for c in cases if not rr[c]]
    facts += [f"{W[len(nomodel)]} attempted cases had ended without a final model", f"Of the {W[len(nomodel)]} cases without a final model"]
    if running: facts += [f"{W[len(running)]} {'was' if len(running) == 1 else 'were'} still running", " and ".join([", ".join(pid(c) for c in running[:-1]), pid(running[-1])] if len(running) > 1 else [pid(running[0])]) + (" were" if len(running) > 1 else " was") + " still running"]
    facts.append(f"The {len(cz)} complete CryoZeta models and the {len(part)} partial models")
    print(f"CryoZeta outcomes: {len(cz)} complete, {len(nomodel)} without a final model ({names(nomodel)}), {len(skipped)} not attempted ({names(skipped)}), {len(running)} no record yet ({names(running)})")
# analysis sentences of section 6.8 (added 2026-10-13): paired differences, results by run mode, precision and size, partial models
if cz:
    Z = "CryoZeta"; tm = lambda m, c: g(m, c, "tm_score"); mode = {c: rec(Z, c).get("mode_used") for c in cz}; nm = [c for c in cz if mode[c] == "normal"]; lg = [c for c in cz if mode[c] != "normal"]
    def pd(m): d = [tm(Z, c) - tm(m, c) for c in cz]; lo, hi = boot(d); return f"{sg(st.mean(d))} ({sg(lo)} to {sg(hi)})"
    d3 = [sg(st.mean(tm(Z, c) - tm(m, c) for c in cz)) for m in M[2:]]
    plan = {x["case"]: x["tokens"] for x in json.load(open(f"{B}/cryozeta/plan.json"))}; th = lambda n: f"{n:,}"
    facts += [f"difference from MICA is {pd('MICA')}", f"from EModelX(+AF) it is {pd('EModelX(+AF)')}", f"({d3[0]}, {d3[1]} and {d3[2]})",
              f"The {W[len(nm)]} cases built in normal mode ({th(min(plan[c] for c in nm))}–{th(max(plan[c] for c in nm))} input residues) all score between {f3(min(tm(Z, c) for c in nm))} and {f3(max(tm(Z, c) for c in nm))}",
              f"(mean {mean(Z, nm)}, against {mean('MICA', nm)} for MICA and {mean('EModelX(+AF)', nm)} for EModelX(+AF) on the same cases)",
              f"above MICA in all {W[len(nm)]}" if all(tm(Z, c) > tm("MICA", c) for c in nm) else f"above MICA in {sum(tm(Z, c) > tm('MICA', c) for c in nm)} of {len(nm)}",
              f"their mean is {mean(Z, lg)} against {mean('MICA', lg)} for MICA, and CryoZeta is above MICA in {W[sum(tm(Z, c) > tm('MICA', c) for c in lg)]}",
              f"the lowest of the six" if all(st.mean(V(Z, "mapframe_precision_3A", cz)) < st.mean(V(m, "mapframe_precision_3A", cz)) for m in M) else "PRECISION RANK CHANGED",
              f"lower mean precision ({f3(st.mean(V(Z, 'mapframe_precision_3A', cz)))}"]
    for c in cz:   # every per-case value quoted in the analysis text
        k = pid(c); q = lambda key, m=Z: f3(g(m, c, key))
        if k == "9WF0": facts.append(f"9WF0 ({q('tm_score')} against {q('tm_score', 'MICA')} for MICA and {q('tm_score', 'EModelX(+AF)')} for EModelX(+AF))")
        if k == "9P0O": facts += [f"({q('tm_score')} against {q('tm_score', 'MICA')};", f"Cα quality there is {q('ca_quality')}"]
        if k == "9W9E": facts += [f"{q('tm_score')} against {q('tm_score', 'MICA')})", f"and {q('ca_quality')}, so on these two cases"]
        if k == "9P80": facts.append(f"9P80: recall {q('mapframe_recall_3A')}, precision {q('mapframe_precision_3A')} with {g(Z, c, 'pred_ca') / g(Z, c, 'ref_ca'):.2f}× the reference")
        if k == "22HP": facts.append(f"precision is {q('mapframe_precision_3A')} and")
        if k == "11OC": facts.append(f"and {q('mapframe_precision_3A')}; where")
    big1 = max(cz, key=lambda c: g(Z, c, "ref_ca")); facts.append(f"{pid(big1)} ({th(g(Z, big1, 'ref_ca'))} reference Cα), was completed with a TM-score of {f3(tm(Z, big1))} against {f3(tm('MICA', big1))} for MICA")
    facts.append(f"The {W[len(lg)]} built in large-complex mode")
    lo3 = sorted(tm(Z, c) for c in lg); facts.append(f"between {f3(lo3[0])} and {f3(max(x for x in lo3 if x < 0.9))}")
    pk = [rec("CryoZeta_partial", c).get("peak_gpu_mib") for c in part]; facts.append(f"peak use {th(min(pk))}–{th(max(pk))} MiB")
    Pm = "CryoZeta_partial"; hiP = [c for c in part if tm(Pm, c) > tm("MICA", c)]; loP = [c for c in part if c not in hiP]
    facts += [f"higher TM-score than MICA's complete model in {W[len(hiP)]} cases ({names(hiP)})", f"a lower one in {W[len(loP)]} (" + "; ".join(f"{pid(c)}: {f3(tm(Pm, c))} against {f3(tm('MICA', c))}" for c in loP) + ")",
              "(" + "; ".join(f"{f3(g(Pm, c, 'ca_quality'))} against {f3(g('MICA', c, 'ca_quality'))}" for c in loP) + ")",
              f"within {max(abs(g(Pm, c, 'ca_quality') - g('MICA', c, 'ca_quality')) for c in hiP) + 0.005:.2f} of MICA's"]
    last = [c for c in part if rec(Pm, c)["stages_total"] - rec(Pm, c)["stages_completed"] == 1]; facts.append(f"in {W[len(last)]} of the {W[len(part)]} cases at the last stage ({names(last)})")
    big = [c for c in part if g(Pm, c, "pred_ca") / g(Pm, c, "ref_ca") > 1.5]
    facts.append("(" + ", ".join(f"{pid(c)} {g(Pm, c, 'pred_ca') / g(Pm, c, 'ref_ca'):.2f}×" for c in big) + ") with map-frame precision of " + " and ".join(f3(g(Pm, c, "mapframe_precision_3A")) for c in big))
    facts.append("(" + " and ".join(th(plan[c]) for c in big) + " against " + " and ".join(th(g(Pm, c, "ref_ca")) for c in big) + ")")
    both = cz + part; mz = st.mean([tm(Z, c) for c in cz] + [tm(Pm, c) for c in part])
    facts.append(f"gives a mean TM-score of {f3(mz)} over {len(both)} cases, against {mean('MICA', both)} for MICA on the same {len(both)}")
    none = [c for c in cases if c not in both]; facts.append(f"({names(none)}; {th(min(g('MICA', c, 'ref_ca') for c in none))}–{th(max(g('MICA', c, 'ref_ca') for c in none))} reference Cα)")
    facts.append(f"its mean TM-score ({mean(Z, cz)}) is above MICA's ({mean('MICA', cz)})"); facts.append(f"all {W[len(nm)]} cases that ran in its normal mode scored {f3(min(tm(Z, c) for c in nm))} or higher")
    if "--facts" in sys.argv: print("\n".join("FACT " + f for f in facts))
p += [f"sentence not found in the text (numbers changed?): '{f}'" for f in facts if f not in txt]
print(f"CryoZeta: {len(cz)} complete ({names(cz)}); {len(part)} partial ({names(part)})")
if "--registry" in sys.argv:   # claim-value registry: every generated table row and every checked sentence, with its status
    os.makedirs(f"{R}/report/qa", exist_ok=True)
    json.dump({"source": "scoring/results/*.json", "tables": {k: e for k, (_, e) in T.items()}, "sentences": [{"text": f, "found_in_report_text": f in txt} for f in facts]},
              open(f"{R}/report/qa/claim_registry.json", "w"), indent=1, ensure_ascii=False)
if p: print(f"\n{len(p)} INCONSISTENCY(IES):"); print("\n".join("  " + x for x in p)); sys.exit(1)
print(f"CONSISTENT: {sum(len(e) for _, e in T.values())} table rows in the text, {N} rows of appendix Table 1 and {len(facts)} sentences match scoring/results")
