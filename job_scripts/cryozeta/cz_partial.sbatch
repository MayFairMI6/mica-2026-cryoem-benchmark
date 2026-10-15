#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 6
#SBATCH --mem=48G
#SBATCH -t 23:30:00
#SBATCH -J czp
#SBATCH -o /home/a/aravindav/bench/logs/czp_%j_%x.out
# CryoZeta large-complex mode with a partial-model fallback (2026-10-13, user decision).
# Usage: sbatch -J czN --dependency=singleton cz_partial.sbatch <CASE>
#
# Why: `cryozeta-cycle-predict` builds one chain copy per stage and exits with an error if one stage cannot be
# registered to the remaining density; the list of completed stages is written only after the last stage, so the
# authors' combine step cannot be run on what was built. Here the SAME pipeline and settings are run with one change
# to cycle_predict.py: the list of completed stages is written at the start of every stage (no computation changes;
# on a complete run the original code overwrites the same files at the end). If the pipeline stops early,
# `cryozeta-combine-stages` is run on the completed stages and the result is recorded as a PARTIAL model.
#
# Records: a complete run -> ~/bench/out/CryoZeta/<CASE>/ (normal result, if none exists yet).
#          a partial run  -> ~/bench/out/CryoZeta_partial/<CASE>/ only (the main record stays "no model").
CASE=$1; B=$HOME/bench; C=$B/envs/cz; cd $C
export PATH=$C/bin:$PATH TMPDIR=$B/tmp SINGULARITY_TMPDIR=$B/tmp SINGULARITY_CACHEDIR=$C/sing_cache
PDB=${CASE##*_}; MAIN=$B/out/CryoZeta/$CASE; PART=$B/out/CryoZeta_partial/$CASE; mkdir -p $MAIN $PART $C/work $C/patch
[ -f $PART/run_record.json ] && { echo "partial-run record exists already"; exit 0; }
J=inputs/json/$PDB.json; [ -f $J ] || { echo "no input json for $CASE"; exit 2; }
GPU=${CUDA_VISIBLE_DEVICES%%,*}; SIF=$C/assets/apptainer/CryoZeta.sif; G=$B/tmp/czp_gpu_$SLURM_JOB_ID.txt; W=$C/work/${PDB}_partial; rm -rf $W
echo "case=$CASE host=$(hostname) gpu=$GPU start=$(date -Is)"

# one-time: copy cycle_predict.py out of the container and add the three recording lines
ACT="cd /app/CryoZeta && source /opt/conda/etc/profile.d/conda.sh && conda activate cryozeta"   # same activation as the authors' wrapper
CP=$(apptainer exec $SIF bash -c "$ACT && python -c 'import cryozeta.runner.cycle_predict as m; print(m.__file__)'" 2>/dev/null | tail -1)
[ -n "$CP" ] || { echo "cannot locate cycle_predict.py in the container"; exit 3; }
( flock 9
  if [ ! -f $C/patch/cycle_predict.py ]; then
    apptainer exec $SIF cat "$CP" > $C/patch/cycle_predict.orig.py
    python3 - $C/patch/cycle_predict.orig.py $C/patch/cycle_predict.py <<'PY'
import sys
s = open(sys.argv[1]).read()
a = "    for stage_idx, (seq_idx, copy_idx, seq_entry) in enumerate(chains, start=1):\n"
assert s.count(a) == 1, "anchor line not found exactly once"
add = ('        # benchmark patch (2026-10-13): record the stages completed so far; the original writes these lists only after the last stage\n'
       '        with open(dump_root / "stage_npz_paths.txt", "w") as _f:\n'
       '            _f.write("\\n".join(str(p) for p in stage_npz_paths if p is not None))\n'
       '        with open(dump_root / "stage_names.txt", "w") as _f:\n'
       '            _f.write("\\n".join(stage_names))\n')
open(sys.argv[2], "w").write(s.replace(a, a + add))
PY
    diff $C/patch/cycle_predict.orig.py $C/patch/cycle_predict.py > $C/patch/cycle_predict.diff; echo "patch written ($(wc -l < $C/patch/cycle_predict.diff) diff lines)"
  fi ) 9>$C/patch/.lock
[ -s $C/patch/cycle_predict.py ] || { echo "patched file missing"; exit 3; }
export SINGULARITY_BIND="$C/patch/cycle_predict.py:$CP" APPTAINER_BIND="$C/patch/cycle_predict.py:$CP"
apptainer exec $SIF grep -c "benchmark patch" "$CP" | sed 's/^/patch visible in container (count): /'

( while :; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i $GPU >> $G 2>/dev/null; sleep 10; done ) & MON=$!
T0=$(date +%s)
timeout 82000 bash assets/apptainer/run_cryozeta_large.sh --sif $SIF --gpu $GPU $J $W $C/assets > $PART/cryozeta_large.log 2>&1; RC=$?
T1=$(date +%s); kill $MON 2>/dev/null
OOM=$(grep -c -i "out of memory" $PART/cryozeta_large.log); [ "$RC" = 124 ] && TO=1 || TO=0
TOTAL=$(grep -o "Processing stage [0-9]*/[0-9]*" $PART/cryozeta_large.log | tail -1 | grep -o "[0-9]*$"); DONE=0; [ -s $W/stage_npz_paths.txt ] && DONE=$(grep -c . $W/stage_npz_paths.txt)
F=$(find $W -maxdepth 1 -name "combined.cif" 2>/dev/null | head -1); PARTIAL=0
if [ -z "$F" ] && [ "$DONE" -gt 0 ]; then
  echo "pipeline stopped (exit $RC) after $DONE of ${TOTAL:-?} stages; combining the completed stages"
  apptainer exec --bind $C:$C --bind $C/assets:/app/CryoZeta/assets $SIF bash -c "$ACT && cryozeta-combine-stages --dump_root $W --output $W/combined_partial.cif" > $PART/combine_partial.log 2>&1
  F=$(ls $W/combined_partial.cif 2>/dev/null); PARTIAL=1
fi
A=0; [ -n "$F" ] && { cp "$F" $PART/$CASE.cif; A=$(grep -c "^ATOM" $PART/$CASE.cif); }
cp $W/stage_npz_paths.txt $PART/ 2>/dev/null
REC="{\"method\":\"CryoZeta\",\"case\":\"$CASE\",\"exit\":$RC,\"seconds\":$((T1-T0)),\"mode_used\":\"large\",\"partial_model\":$PARTIAL,\"stages_completed\":$DONE,\"stages_total\":${TOTAL:-0},\"gpu_oom_lines\":$OOM,\"timed_out\":$TO,\"peak_gpu_mib\":$(sort -n $G | tail -1),\"final_model\":\"$(basename "$F")\",\"final_atoms\":$A,\"normal_mode_oom\":0,\"mode_requested\":\"large\",\"host\":\"$(hostname)\",\"gpu\":\"NVIDIA TITAN RTX 24GB\",\"cycle_predict_patch\":\"completed-stage list written per stage\"}"
echo "$REC" | tee $PART/run_record.json
if [ ! -f $MAIN/run_record.json ]; then          # case had no result under the unmodified pipeline yet
  if [ "$PARTIAL" = 0 ] && [ "$A" -gt 0 ]; then cp $PART/$CASE.cif $MAIN/$CASE.cif; echo "$REC" > $MAIN/run_record.json
  else echo "$REC" | sed "s/\"final_atoms\":$A/\"final_atoms\":0/; s/\"final_model\":\"[^\"]*\"/\"final_model\":\"\"/" > $MAIN/run_record.json; fi
  cp $PART/cryozeta_large.log $MAIN/cryozeta_large.log
fi
rm -f $G; tail -n 3 $PART/cryozeta_large.log | cut -c1-250; rm -rf $W
