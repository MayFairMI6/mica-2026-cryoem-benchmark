#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -J cryoatom
#SBATCH -o /home/a/aravindav/bench/logs/cryoatom_%A_%a.out
# CryoAtom v1.0.0 @e2a31f2b, README usage: cryoatom build -s FASTA -v MAP -o OUT -d cuda:0 (defaults otherwise).
# Map = canonical.mrc (CryoAtom cubes/resamples internally). ESM weights cached under TORCH_HOME on first use.
METHOD=CryoAtom; source $HOME/bench/common.sh
export TORCH_HOME=$B/weights/torch_hub OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
peak_gpu; T0=$(date +%s)
cd $OUT && $E/ca/bin/cryoatom build -s $FASTA -v $MAP -o $OUT/$CASE -d cuda:0; RC=$?
cp $OUT/$CASE/$CASE.cif $OUT/$CASE/${CASE}_raw.cif $OUT/ 2>/dev/null
record $RC $OUT/$CASE.cif $(( $(date +%s)-T0 ))
