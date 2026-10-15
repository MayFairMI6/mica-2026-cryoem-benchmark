#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 23:00:00
#SBATCH -J af2
#SBATCH -o /home/a/aravindav/bench/logs/af2_%A_%a.out
# EModelX(+AF) prerequisite, same protocol as rounds 1-2: ColabFold 1.6.3, 5 models, 3 recycles, seed 0, MMseqs2 server;
# one single-chain model per entity, rank_001 kept as template_dir/<PDBID_n>/ranked_0.pdb and sequence-validated.
METHOD=af2_templates; source $HOME/bench/common.sh
export XDG_CACHE_HOME=$B/weights/cf_cache TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
mkdir -p $OUT/in
$E/cf/bin/python - $FASTA $OUT/in <<'PY'
import sys
for b in open(sys.argv[1]).read().strip().split(">")[1:]:
    h, *s = b.strip().split("\n"); e = h.split("|")[0]; open(f"{sys.argv[2]}/{e}.fasta", "w").write(f">{e}\n{''.join(s).strip()}\n")
PY
T0=$(date +%s); RC=0
for f in $OUT/in/*.fasta; do e=$(basename $f .fasta)
  [ -f $OUT/$e/ranked_0.pdb ] && continue
  $E/cf/bin/colabfold_batch --num-models 5 --num-recycle 3 --random-seed 0 $f $OUT/cf_$e || RC=1
  P=$(ls $OUT/cf_$e/${e}*_rank_001_*.pdb 2>/dev/null | head -1); [ -n "$P" ] && mkdir -p $OUT/$e && cp $P $OUT/$e/ranked_0.pdb && cp $OUT/cf_$e/${e}*scores_rank_001*.json $OUT/$e/ 2>/dev/null
done
$E/cf/bin/python - $OUT <<'PY'
import glob, json, os, sys
three = dict(ALA="A",ARG="R",ASN="N",ASP="D",CYS="C",GLN="Q",GLU="E",GLY="G",HIS="H",ILE="I",LEU="L",LYS="K",MET="M",PHE="F",PRO="P",SER="S",THR="T",TRP="W",TYR="Y",VAL="V")
o = sys.argv[1]; rec = {}
for f in glob.glob(f"{o}/in/*.fasta"):
    e = os.path.basename(f)[:-6]; seq = open(f).read().split("\n", 1)[1].replace("\n", ""); p = f"{o}/{e}/ranked_0.pdb"
    if not os.path.exists(p): rec[e] = {"template": None}; continue
    got = [three.get(l[17:20], "X") for l in open(p) if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    pl = [float(l[60:66]) for l in open(p) if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    rec[e] = {"length": len(seq), "sequence_matches_fasta": "".join(got) == seq, "mean_plddt": round(sum(pl) / len(pl), 1)}
json.dump(rec, open(f"{o}/templates_record.json", "w"), indent=1); print(rec)
PY
echo "af2 exit $RC seconds $(( $(date +%s)-T0 ))"
