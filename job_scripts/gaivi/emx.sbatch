#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -J emx
#SBATCH -o /home/a/aravindav/bench/logs/emx_%A_%a.out
# EModelX(+AF) @53713cf8, frozen settings: --protocol=temp_flex, canonical.mrc, deposited FASTA, exact-sequence AF2
# templates (template_dir/<PDBID_n>/ranked_0.pdb), --seed 2022 --mul_proc_num 4 --run_pulchra. Never run with a partial template set.
METHOD="EModelX(+AF)"; source $HOME/bench/common.sh
T=$B/out/af2_templates/$CASE
for e in $(grep ">" $FASTA | sed 's/>//; s/|.*//'); do [ -f $T/$e/ranked_0.pdb ] || { echo "missing template $e"; record 2 none 0; exit 0; }; done
mkdir -p $OUT/pred; peak_gpu; T0=$(date +%s); cd $E/EModelX
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 $E/emx/bin/python run.py --protocol=temp_flex --EM_map $MAP --fasta $FASTA --template_dir $T \
  --output_dir $OUT/pred --seed 2022 --mul_proc_num 4 --run_pulchra --pulchra_path modules/pulchra304/src/pulchra; RC=$?
F=$(ls $OUT/pred/*_EModelX\(+AF\)_all_atom_model.pdb $OUT/pred/*all_atom*.pdb 2>/dev/null | head -1); [ -n "$F" ] && cp "$F" $OUT/
record $RC "$OUT/$(basename "$F")" $(( $(date +%s)-T0 ))
