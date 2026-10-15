"""Case-level low-homology summary for the extension cases (rounds 4-6), same rule as the core file
audit3_low_homology_by_case.json: an entity is 'low' when its best pre-cutoff PDB hit (max identity x coverage, from
audit3_homology_rcsb_ext.json) does not reach identity >= 0.5 AND query coverage >= 0.5; the case-level value is the
fraction of residues (entity length x copies, copies from the input FASTA header) in low entities.
Writes audit3_low_homology_by_case_ext.json. No result file is read."""
import json, os, re
B = os.path.dirname(os.path.abspath(__file__))
aud = json.load(open(f"{B}/audit3_homology_rcsb_ext.json"))
low = lambda best: best is None or not (best["identity"] >= 0.5 and best["query_coverage"] >= 0.5)
out = {}
for case in [l.strip() for l in open(f"{B}/gaivi_bench/cases_ext.txt") if l.strip()]:
    pdb = case.split("_")[1]; ents = []
    for b in open(f"{B}/ext_inputs/{case}/{pdb}.fasta").read().strip().split(">")[1:]:
        h, *s = b.split("\n"); eid = h.split("|")[0]; seq = "".join(s).strip()
        copies = len(re.sub(r"Chains?\s*", "", h.split("|")[1]).split(","))
        a = aud[eid]; assert "error" not in a["pre_2023-05"] and "error" not in a["pre_2026"], eid
        assert a["length"] == len(seq), eid
        ents.append({"id": eid, "len": len(seq), "copies": copies,
                     "low_pre2023_05": low(a["pre_2023-05"]["best"]), "best_pre2023_05": a["pre_2023-05"]["best"],
                     "low_pre2026": low(a["pre_2026"]["best"]), "best_pre2026": a["pre_2026"]["best"]})
    tot = sum(e["len"] * e["copies"] for e in ents)
    out[case] = {"frac_residues_low_pre2023_05": round(sum(e["len"] * e["copies"] for e in ents if e["low_pre2023_05"]) / tot, 3),
                 "frac_residues_low_pre2026": round(sum(e["len"] * e["copies"] for e in ents if e["low_pre2026"]) / tot, 3), "entities": ents}
    print(case, out[case]["frac_residues_low_pre2023_05"], out[case]["frac_residues_low_pre2026"])
json.dump(out, open(f"{B}/audit3_low_homology_by_case_ext.json", "w"), indent=1)
print("low-homology (>= 0.5, pre-2023-05):", [c.split("_")[1] for c, v in out.items() if v["frac_residues_low_pre2023_05"] >= 0.5])
