#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 6
#SBATCH --mem=48G
#SBATCH -t 23:30:00
#SBATCH -J cz
#SBATCH -o /home/a/aravindav/bench/logs/cz_%j_%x.out
# CryoZeta (post hoc addition, 2026-10-12): authors container image + wrapper scripts, default settings
# (seed 101, 5 samples, 20 steps, 10 cycles). Usage: sbatch -J czN --dependency=singleton cz_run.sbatch <CASE> <normal|large>
# normal = run_cryozeta.sh (README: up to ~2,800 residues); large = run_cryozeta_large.sh. If normal mode runs out of
# GPU memory on the 24 GB card, large mode is tried once in the same job and the record says which mode gave the model.
CASE=$1; MODE=$2; B=$HOME/bench; C=$B/envs/cz; cd $C
export PATH=$C/bin:$PATH TMPDIR=$B/tmp SINGULARITY_TMPDIR=$B/tmp SINGULARITY_CACHEDIR=$C/sing_cache
PDB=${CASE##*_}; OUT=$B/out/CryoZeta/$CASE; mkdir -p $OUT $C/work
[ -f $OUT/run_record.json ] && { echo "recorded already"; exit 0; }
J=inputs/json/$PDB.json; [ -f $J ] || { echo "no input json for $CASE"; exit 2; }
# Known-limit rule (2026-10-12): a case is not attempted if a case of the same or smaller size has already hit a
# measured hardware limit on this GPU (23.5 h time limit, or GPU out of memory in large-complex mode). Predicted limits
# are not used; other kinds of failure do not trigger it.
SKIP=$(python3 - "$CASE" <<PY
import glob, json, os, sys
B = os.path.expanduser("~/bench"); plan = {p["case"]: p for p in json.load(open(B + "/envs/cz/inputs/plan.json"))}; me = plan[sys.argv[1]]
for f in sorted(glob.glob(B + "/out/CryoZeta/*/run_record.json")):
    r = json.load(open(f)); o = plan.get(r.get("case"))
    if not o or r.get("final_atoms", 0) > 0 or r.get("skipped"): continue
    # 2026-10-13: large-complex mode predicts one chain per stage, so its memory limit depends on the longest chain,
    # not on the total size (9U62 ran out of memory on its 1,663-residue chain). The memory branch compares longest chains.
    mx = lambda q: max(len(list(e.values())[0]["sequence"]) for e in json.load(open(B + "/envs/cz/inputs/json/%s.json" % q["pdb"]))[0]["sequences"])
    if r.get("timed_out") == 1 and o["tokens"] <= me["tokens"]:
        print("not attempted: larger than %s (%d residues), which exceeded the 23.5 h limit on the 24 GB GPU" % (o["pdb"], o["tokens"])); break
    if r.get("mode_used") == "large" and r.get("gpu_oom_lines", 0) > 0 and mx(o) <= mx(me):
        print("not attempted: its longest chain (%d residues) is at least as long as that of %s (%d), which ran out of GPU memory in large-complex mode on the 24 GB GPU" % (mx(me), o["pdb"], mx(o))); break
PY
)
if [ -n "$SKIP" ]; then echo "{\"method\":\"CryoZeta\",\"case\":\"$CASE\",\"exit\":-1,\"seconds\":0,\"skipped\":\"$SKIP\",\"final_atoms\":0}" | tee $OUT/run_record.json; exit 0; fi
GPU=${CUDA_VISIBLE_DEVICES%%,*}; SIF=$C/assets/apptainer/CryoZeta.sif; G=$B/tmp/cz_gpu_$SLURM_JOB_ID.txt
echo "case=$CASE mode=$MODE host=$(hostname) gpu=$GPU start=$(date -Is)"; nvidia-smi --query-gpu=name,memory.total --format=csv,noheader -i $GPU
( while :; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i $GPU >> $G 2>/dev/null; sleep 10; done ) & MON=$!
run() { W=$C/work/${PDB}_$1; rm -rf $W
  if [ "$1" = normal ]; then timeout 82000 bash assets/apptainer/run_cryozeta.sh --sif $SIF --gpu $GPU $J $W $C/assets > $OUT/cryozeta_$1.log 2>&1; RC=$?
    F=$(ls $W/$PDB/CryoZeta-Final/${PDB}_sample_0.cif 2>/dev/null)
  else timeout 82000 bash assets/apptainer/run_cryozeta_large.sh --sif $SIF --gpu $GPU $J $W $C/assets > $OUT/cryozeta_$1.log 2>&1; RC=$?
    F=$(find $W -maxdepth 2 -name "combined*.cif" 2>/dev/null | head -1); fi
  OOM=$(grep -c -i "out of memory" $OUT/cryozeta_$1.log); }
T0=$(date +%s); run $MODE; USED=$MODE
if [ -z "$F" ] && [ "$MODE" = normal ] && [ "$OOM" -gt 0 ]; then echo "normal mode: GPU out of memory after $(( $(date +%s)-T0 )) s; trying large mode"; FIRST_OOM=1; run large; USED=large; fi
T1=$(date +%s); kill $MON 2>/dev/null
A=0; if [ -n "$F" ]; then cp "$F" $OUT/$CASE.cif; A=$(grep -c "^ATOM" $OUT/$CASE.cif); mkdir -p $OUT/samples; cp $(dirname "$F")/*.cif $OUT/samples/ 2>/dev/null; fi
[ "$RC" = 124 ] && TO=1 || TO=0
echo "{\"method\":\"CryoZeta\",\"case\":\"$CASE\",\"exit\":$RC,\"seconds\":$((T1-T0)),\"mode_requested\":\"$MODE\",\"mode_used\":\"$USED\",\"normal_mode_oom\":${FIRST_OOM:-0},\"gpu_oom_lines\":$OOM,\"timed_out\":$TO,\"peak_gpu_mib\":$(sort -n $G | tail -1),\"final_model\":\"$(basename "$F")\",\"final_atoms\":$A,\"host\":\"$(hostname)\",\"gpu\":\"NVIDIA TITAN RTX 24GB\"}" | tee $OUT/run_record.json
rm -f $G; tail -n 5 $OUT/cryozeta_$USED.log | cut -c1-300; rm -rf $C/work/${PDB}_normal $C/work/${PDB}_large
