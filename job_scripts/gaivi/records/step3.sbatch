#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 23:00:00
#SBATCH -J mica_s3
#SBATCH -o /home/a/aravindav/mica_bench/logs/step3_%A_%a.out
# MICA step 3 = MICA_pipeline.sh step 3 verbatim (run.py ... --run_pulchra --run_phenix --device cuda) on Phenix-docked inputs.
# Recorded deviation (same as Kaggle): batch_threshold raised -> batch 1. Hardware: GAIVI Titan RTX 24 GB (Kaggle: T4 16 GB).
B=$HOME/mica_bench; W=$B/work; EMD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $B/${CASEFILE:-cases.txt})
read FASTA RES CASE < <($W/venv/bin/python -c "import json;m=json.load(open(\"$B/inputs/manifest.json\"))[\"$EMD\"];print(m[\"fasta\"],m[\"resolution_A\"],m[\"case\"])")
[ -f $B/logs/dock_record_$EMD.json ] && ! grep -q "\"step2_exit\":\"0\"" $B/logs/dock_record_$EMD.json && { echo "docking incomplete for $EMD; not running step 3"; exit 1; }
J=$B/jobtmp/$SLURM_JOB_ID; mkdir -p $J; cp -r $W/MICA $J/; cd $J/MICA; mkdir -p input  # node /tmp is a 2 GB shared tmpfs; use NFS home
tar -xzf $B/docked/mica_input_$EMD.tar.gz -C input
echo "case=$CASE emd=$EMD res=$RES fasta=$FASTA host=$(hostname) gpu=$CUDA_VISIBLE_DEVICES start=$(date -Is)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
export PATH=$W/venv/bin:$PATH OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled
( while :; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i ${CUDA_VISIBLE_DEVICES%%,*} >> $J/gpu_mem.txt 2>/dev/null; sleep 5; done ) & MON=$!
T0=$(date +%s)
python run.py -m input/$EMD/emd_$EMD.map -f input/$EMD/$FASTA -a input/$EMD/AF3_results --run_pulchra \
  --pulchra_path=modules/pulchra304/src/pulchra --run_phenix --phenix_act=$W/phenix/phenix-1.20.1-4487/phenix_env.sh \
  --resolution=$RES --device=cuda -o $J/out
RC=$?; T1=$(date +%s); kill $MON
O=$B/outputs/$EMD; mkdir -p $O; cp $J/out/* $O/ 2>/dev/null
PEAK=$(sort -n $J/gpu_mem.txt | tail -1)
F=$(ls $O/*_MICA_all_atom_model.pdb 2>/dev/null | head -1); A=0; [ -n "$F" ] && A=$(grep -cE "^(ATOM|HETATM)" $F)
echo "{\"case\":\"$CASE\",\"emd\":\"$EMD\",\"exit\":$RC,\"seconds\":$((T1-T0)),\"peak_gpu_mib\":${PEAK:-0},\"final_model\":\"$(basename "$F")\",\"final_atoms\":$A,\"host\":\"$(hostname)\"}" | tee $O/run_record.json
rm -rf $J
