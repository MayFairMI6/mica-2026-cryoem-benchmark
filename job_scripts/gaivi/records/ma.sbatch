#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -J ma
#SBATCH -o /home/a/aravindav/bench/logs/ma_%A_%a.out
# ModelAngelo 1.0.18 @994945bd, frozen settings: canonical_cube.mrc, deposited FASTA, build ... --device 0.
METHOD=ModelAngelo; source $HOME/bench/common.sh
export TORCH_HOME=$B/weights/modelangelo OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
peak_gpu; T0=$(date +%s)
$E/ma/bin/model_angelo build -v $CUBE -pf $FASTA -o $OUT/$CASE --device 0; RC=$?
cp $OUT/$CASE/$CASE.cif $OUT/$CASE/${CASE}_raw.cif $OUT/ 2>/dev/null
record $RC $OUT/$CASE.cif $(( $(date +%s)-T0 ))
