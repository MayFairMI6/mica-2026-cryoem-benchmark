#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -J emprot
#SBATCH -o /home/a/aravindav/bench/logs/emprot_%A_%a.out
# EMProt v1.2 @9489d52e, map+sequence mode (README example 1): emprot build --map MAP --seq FASTA --device 0 --output OUT.
# Map = canonical.mrc; deposited FASTA; no predicted models (--chain not used). Final model = OUT/output.cif.
METHOD=EMProt; source $HOME/bench/common.sh
source /apps/anaconda3/etc/profile.d/conda.sh; conda activate $E/emprot_env
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
peak_gpu; T0=$(date +%s)
emprot build --map $MAP --seq $FASTA --device 0 --output $OUT/run > $OUT/emprot.log 2>&1; RC=$?; cat $OUT/emprot.log
BIN=$E/EMProt/emprot/bin
if [ $RC -ne 0 ] && grep -qiE "segmentation fault|illegal instruction|cannot convert|core dumped" $OUT/emprot.log && [ ! -d $BIN/static_original ]; then
  # README "Notes": static binaries may crash on non-CentOS7/Intel hosts; use emprot/bin/recompiled instead (done once, for all cases)
  echo "static binary crashed -> switching to recompiled binaries (README guidance) and retrying once"
  mkdir -p $BIN/static_original && for f in $BIN/recompiled/*; do b=$(basename $f); [ -f $BIN/$b ] && cp $BIN/$b $BIN/static_original/ && cp $f $BIN/$b && chmod +x $BIN/$b; done
  rm -rf /home/a/aravindav/bench/out/EMProt/$CASE/run
  emprot build --map $MAP --seq $FASTA --device 0 --output $OUT/run > $OUT/emprot_retry.log 2>&1; RC=$?; cat $OUT/emprot_retry.log
fi
cp $OUT/run/output.cif $OUT/$CASE.cif 2>/dev/null; cp $OUT/run/output_denovo.cif $OUT/${CASE}_denovo.cif 2>/dev/null
record $RC $OUT/$CASE.cif $(( $(date +%s)-T0 ))
