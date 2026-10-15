#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 23:30:00
#SBATCH -J mica_dock
#SBATCH -o /home/a/aravindav/mica_bench/logs/dock_%j_%x.out
# MICA steps 1-2 (MICA_pipeline.sh verbatim): process_AF3_results.py, then dock_in_map.py (Phenix, domainwise).
# Usage: sbatch -J dock_<EMD> dock_ext.sbatch <EMD> [resume_tarball] [attempt]
# Same recorded deviations as the core run: affinity-aware cpu_count shim (nproc = 0.6 x allocated cores = 9);
# native resume from existing docked files; domains that ended "No solution found" on the same masked map not retried.
# Extension-only housekeeping (no MICA file changed): dock_in_map.py keeps one full masked-map copy per docked domain
# in temp_maps/ until it exits (tens of GB for 600^3 maps). A background loop deletes all but the two newest copies;
# dock_in_map.py only ever reads the newest one. On a 22 h timeout the job resubmits itself from its saved state;
# on success it submits step 3.
EMD=$1; PREV=$2; ATT=${3:-1}; B=$HOME/mica_bench; W=$B/work
read FASTA RES CON CASE < <($W/venv/bin/python -c "import json;m=json.load(open('$B/inputs/manifest.json'))['$EMD'];print(m['fasta'],m['resolution_A'],m['contour_level'],m['case'])")
[ -n "$CASE" ] || { echo "no manifest entry for $EMD"; exit 2; }
J=$B/jobtmp/dock_$SLURM_JOB_ID; mkdir -p $J; cp -r $W/MICA $J/; cd $J/MICA; mkdir -p input
if [ -n "$PREV" ]; then tar -xzf $PREV -C input; rm -rf input/$EMD/temp_maps
else cp -r $B/r3_inputs/$EMD input/; cp $HOME/bench/maps/$CASE/canonical.mrc input/$EMD/emd_$EMD.map; fi
FAILED=$(grep -l "No solution found" input/$EMD/docking_logs/phenix_fold*.log 2>/dev/null | sed 's|.*/phenix_||; s|\.pdb_.*|.pdb|')
mkdir -p input/$EMD/AF3_domains_failed; for f in $FAILED; do [ -f input/$EMD/AF3_domains/$f ] && mv input/$EMD/AF3_domains/$f input/$EMD/AF3_domains_failed/; done
echo "case=$CASE emd=$EMD res=$RES contour=$CON attempt=$ATT prev=${PREV:-none} skipped_failed=[$FAILED] host=$(hostname) start=$(date -Is)"
export PYTHONPATH=$B/shim PATH=$W/venv/bin:$PATH OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled
( while :; do ls -t input/$EMD/temp_maps/map_after_*.mrc 2>/dev/null | tail -n +3 | while read f; do rm -f "$f"; done; sleep 60; done ) & CLEAN=$!
T0=$(date +%s)
if [ -z "$PREV" ]; then python utils/process_AF3_results.py -f input/$EMD/$FASTA -a input/$EMD/AF3_results; echo "step1 exit $?"; fi
T1=$(date +%s)
timeout 22h python utils/dock_in_map.py -m input/$EMD/emd_$EMD.map -c $CON -r $RES -f input/$EMD/$FASTA -a input/$EMD/AF3_results --phenix_act $W/phenix/phenix-1.20.1-4487/phenix_env.sh
RC=$?; T2=$(date +%s); kill $CLEAN 2>/dev/null
mv input/$EMD/AF3_domains_failed/* input/$EMD/AF3_domains/ 2>/dev/null; rmdir input/$EMD/AF3_domains_failed; rm -rf input/$EMD/temp_maps
tar -czf $B/docked/mica_input_$EMD.tar.gz -C input $EMD
N=$(ls input/$EMD/AF3_domains/*.pdb | wc -l); D=$(ls input/$EMD/AF3_docked_models/*.pdb 2>/dev/null | wc -l)
NOSOL=$(grep -l "No solution found" input/$EMD/docking_logs/phenix_fold*.log 2>/dev/null | sed 's|.*/phenix_||; s|\.pdb_.*||' | sort -u | wc -l)
UNEXPL=$((N-D-NOSOL)); [ $RC -eq 0 ] && [ $UNEXPL -gt 0 ] && RC=99
PRIOR=$(python -c "import json;print(json.load(open('$B/logs/dock_record_$EMD.json')).get('step2_s_total',0))" 2>/dev/null || echo 0)
echo "{\"case\":\"$CASE\",\"attempt\":$ATT,\"no_solution\":$NOSOL,\"unexplained_undocked\":$UNEXPL,\"step1_s\":$((T1-T0)),\"step2_exit\":\"$RC\",\"step2_s\":$((T2-T1)),\"step2_s_total\":$((PRIOR+T2-T1)),\"domains\":$N,\"docked\":$D,\"combined\":$( [ -f input/$EMD/${EMD}_af3_docked.pdb ] && echo true || echo false),\"host\":\"$(hostname)\",\"cores\":16}" | tee $B/logs/dock_record_$EMD.json
cd $B; rm -rf /home/a/aravindav/mica_bench/jobtmp/dock_$SLURM_JOB_ID
if [ "$RC" = "0" ]; then MEM=$(awk -v e=$EMD '$1==e{print $2}' $B/step3_mem.txt); sbatch -J s3_$EMD --mem=${MEM:-48}G $B/step3_ext.sbatch $EMD
elif [ "$RC" = "124" ] && [ $ATT -lt 4 ]; then sbatch -J dock_$EMD $B/dock_ext.sbatch $EMD $B/docked/mica_input_$EMD.tar.gz $((ATT+1))
else echo "docking for $EMD stopped with exit $RC after attempt $ATT; not continuing"; fi
exit 0
