"""Evaluation-side format normalisation for EMProt v1.2 output (no coordinate is changed; reference never used).

EMProt's pipeline/reorderx.py writes every built fragment that matches the same INPUT SEQUENCE into one chain and
separates fragments only by a residue-number jump (delta = 10). For a sequence present in several copies, all copies
therefore share one chain ID, which a chain-level complex TM-score cannot interpret.
Rule (uses only the input FASTA that EMProt itself was given):
  * chain matched to a sequence with ONE copy  -> kept as EMProt wrote it (it is one real chain);
  * chain matched to a sequence with >1 copies -> each continuous fragment (residue-number jump != 1) gets its own chain ID;
  * mode "all" splits every chain; mode "none" is EMProt's raw file. All three are scored and reported.
Usage: emprot_split_chains.py IN.cif OUT.cif FASTA [copyaware|all]"""
import re, sys
from difflib import SequenceMatcher
sys.path.insert(0, __file__.rsplit("/", 1)[0] + "/vendor")
import gemmi
inp, outp, fasta = sys.argv[1:4]; mode = sys.argv[4] if len(sys.argv) > 4 else "copyaware"
ents = []
for b in open(fasta).read().strip().split(">")[1:]:
    h, *s = b.split("\n"); ents.append((h.split("|")[0], len(re.sub(r"Chains?\s*", "", h.split("|")[1]).split(",")), "".join(s).strip()))
st = gemmi.read_structure(inp); model = st[0]; new = gemmi.Model("1"); n = 0; rep = []
for ch in model:
    frags, prev = [], None
    for res in ch:
        if prev is None or res.seqid.num != prev + 1: frags.append([])
        frags[-1].append(res); prev = res.seqid.num
    seq = "".join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for f in frags for r in f)
    # entity = the input sequence that explains most residues of this chain (longest-matching-block coverage per fragment)
    def cover(es):
        return sum(sum(b.size for b in SequenceMatcher(None, "".join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for r in f), es, autojunk=False).get_matching_blocks()) for f in frags)
    eid, copies, _ = max(ents, key=lambda e: cover(e[2]))
    split = mode == "all" or (mode == "copyaware" and copies > 1)
    groups = frags if split else [[r for f in frags for r in f]]
    for grp in groups:
        n += 1; new.add_chain(gemmi.Chain(f"F{n}")); c = new[len(new) - 1]
        for r in grp: c.add_residue(r)
    rep.append(f"{ch.name}:{eid}x{copies}:{len(frags)}frag->{'split' if split else 'kept'}")
out = gemmi.Structure(); out.cell = st.cell; out.spacegroup_hm = st.spacegroup_hm; out.name = st.name
out.add_model(new); out.setup_entities(); out.assign_label_seq_id()
a0 = sum(1 for _ in st[0].all()); a1 = sum(1 for _ in out[0].all()); assert a0 == a1, (a0, a1)
out.make_mmcif_document().write_file(outp)
print(f"{inp.split('/')[-1]} [{mode}] atoms {a1} unchanged, chains {len(st[0])} -> {len(out[0])}: " + " ".join(rep))
