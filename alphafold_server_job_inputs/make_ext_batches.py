"""AlphaFold Server job lists for extension rounds 4-6, built with MICA's own fasta_to_AF3_json.py (@726d1bb1) exactly
as for rounds 1-3: one job per entity, one proteinChain entry per chain copy, server-chosen seed.
Recorded deviation (same rule as 9WF0_1 in round 3): an entity whose copies exceed the server's 5,000-token limit is
submitted as a single copy."""
import csv, json, os, re, sys, logging
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
logging.disable(logging.CRITICAL)
from fasta_to_AF3_json import FastaToJsonConverter
B = "/Users/anony/Downloads/MICA_GAIVI_run/benchmark"; conv = FastaToJsonConverter()
rounds = {r["case_id"]: r["round"] for r in csv.DictReader(open(f"{B}/case_rounds.csv"))}
out = {"4": [], "5": [], "6": []}; notes = []
for case, rd in rounds.items():
    if rd not in out: continue
    pdb = case.split("_")[1]; blocks = open(f"{B}/ext_inputs/{case}/{pdb}.fasta").read().strip().split(">")[1:]
    blocks.sort(key=lambda b: int(b.split("|")[0].rsplit("_", 1)[1]))           # README: entities in entity-number order
    ordered = "".join(">" + b.strip() + "\n" for b in blocks); open(f"{B}/af3_jobs/fasta_ordered/{pdb.lower()}.fasta", "w").write(ordered)
    for p in conv.process_files(ordered, pdb):
        chains, seq = p["chains"], p["sequence"]
        if len(seq) * len(chains) > 5000:
            notes.append(f"{p['protein_id']}: {len(chains)} x {len(seq)} = {len(seq) * len(chains)} tokens > 5000 -> single copy"); chains = chains[:1]
        jobs = conv.generate_json_multichain(p["protein_id"], chains, seq)
        if jobs is None: notes.append(f"{p['protein_id']}: converter produced no job (skipped by its own filters)"); continue
        out[rd] += jobs
for rd, jobs in out.items():
    json.dump(jobs, open(f"{B}/af3_jobs/batch_ext_round{rd}.json", "w"))
    print("round", rd, len(jobs), "jobs:", [(j["name"], len(j["sequences"]), len(j["sequences"][0]["proteinChain"]["sequence"])) for j in jobs])
print("\n".join(notes)); json.dump(notes, open(f"{B}/af3_jobs/batch_ext_notes.json", "w"), indent=1)
