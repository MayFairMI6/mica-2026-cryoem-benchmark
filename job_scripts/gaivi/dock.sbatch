#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 23:30:00
#SBATCH -J mica_dock
#SBATCH -o /home/a/aravindav/mica_bench/logs/dock_%j_%x.out
# MICA steps 1-2 (MICA_pipeline.sh verbatim): process_AF3_results.py, then dock_in_map.py (Phenix, domainwise).
# Usage: sbatch -J dock_<EMD> dock.sbatch <EMD> [resume_tarball]
# dock_in_map.py sets nproc = 0.6 x multiprocessing.cpu_count(). Under Slurm that call sees all 64 node cores (-> 38 Phenix
# workers on 16 allocated cores; fork() ran out of memory). shim/sitecustomize.py makes it affinity-aware, so nproc =
# 0.6 x 16 = 9: the authors' formula on the allocated cores (Kaggle: 0.6 x 4 -> 2). No MICA code is changed.
# Resume: dock_in_map.py reuses existing AF3_docked_models/*_docked.pdb (native). Domains that already ended
# "No solution found" against the same masked map are not retried (recorded deviation, as on Kaggle).
EMD=$1; PREV=$2; B=$HOME/mica_bench; W=$B/work
read FASTA RES CON CASE < <($W/venv/bin/python -c "import json;m=json.load(open('$B/inputs/manifest.json'))['$EMD'];print(m['fasta'],m['resolution_A'],m['contour_level'],m['case'])")
J=$B/jobtmp/dock_$SLURM_JOB_ID; mkdir -p $J; cp -r $W/MICA $J/; cd $J/MICA; mkdir -p input
if [ -n "$PREV" ]; then tar -xzf $PREV -C input; rm -rf input/$EMD/temp_maps
else cp -r $B/r3_inputs/$EMD input/; cp $HOME/bench/maps/$CASE/canonical.mrc input/$EMD/emd_$EMD.map; fi
FAILED=$(grep -l "No solution found" input/$EMD/docking_logs/phenix_fold*.log 2>/dev/null | sed 's|.*/phenix_||; s|\.pdb_.*|.pdb|')
mkdir -p input/$EMD/AF3_domains_failed; for f in $FAILED; do [ -f input/$EMD/AF3_domains/$f ] && mv input/$EMD/AF3_domains/$f input/$EMD/AF3_domains_failed/; done
echo "case=$CASE emd=$EMD res=$RES contour=$CON prev=${PREV:-none} skipped_failed=[$FAILED] host=$(hostname) start=$(date -Is)"
export PYTHONPATH=$B/shim PATH=$W/venv/bin:$PATH OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled
T0=$(date +%s)
if [ -z "$PREV" ]; then python utils/process_AF3_results.py -f input/$EMD/$FASTA -a input/$EMD/AF3_results; echo "step1 exit $?"; fi
T1=$(date +%s)
timeout 22h python utils/dock_in_map.py -m input/$EMD/emd_$EMD.map -c $CON -r $RES -f input/$EMD/$FASTA -a input/$EMD/AF3_results --phenix_act $W/phenix/phenix-1.20.1-4487/phenix_env.sh
RC=$?; T2=$(date +%s)
mv input/$EMD/AF3_domains_failed/* input/$EMD/AF3_domains/ 2>/dev/null; rmdir input/$EMD/AF3_domains_failed; rm -rf input/$EMD/temp_maps
tar -czf $B/docked/mica_input_$EMD.tar.gz -C input $EMD
N=$(ls input/$EMD/AF3_domains/*.pdb | wc -l); D=$(ls input/$EMD/AF3_docked_models/*.pdb 2>/dev/null | wc -l)
NOSOL=$(grep -l "No solution found" input/$EMD/docking_logs/phenix_fold*.log 2>/dev/null | sed 's|.*/phenix_||; s|\.pdb_.*||' | sort -u | wc -l)
UNEXPL=$((N-D-NOSOL))  # domains neither docked nor a completed "No solution" search (e.g. out of memory, timeout)
[ $RC -eq 0 ] && [ $UNEXPL -gt 0 ] && RC=99
echo "{\"case\":\"$CASE\",\"no_solution\":$NOSOL,\"unexplained_undocked\":$UNEXPL,\"step1_s\":$((T1-T0)),\"step2_exit\":\"$RC\",\"step2_s\":$((T2-T1)),\"domains\":$N,\"docked\":$D,\"combined\":$( [ -f input/$EMD/${EMD}_af3_docked.pdb ] && echo true || echo false)}" | tee $B/logs/dock_record_$EMD.json
rm -rf /home/a/aravindav/mica_bench/jobtmp/dock_$SLURM_JOB_ID
exit $RC
