#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 23:00:00
#SBATCH -J finalize
#SBATCH -o /home/a/aravindav/bench/logs/finalize_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aravindav@usf.edu
# Waits for the EMProt setup + array (polling squeue, so it also ends if setup fails and the array can never start),
# bundles every result, verifies the bundle, then deletes only re-creatable bulk (environments, weights, scratch).
SETUP=$1; ARRAY=$2; B=$HOME/bench; MB=$HOME/mica_bench; STAMP=$(date +%Y%m%d_%H%M)
while :; do
  S=$(squeue -h -j $SETUP 2>/dev/null | wc -l); A=$(squeue -h -r -j $ARRAY -o "%T %r" 2>/dev/null)
  [ "$S" = "0" ] || { sleep 300; continue; }
  if echo "$A" | grep -q "DependencyNeverSatisfied"; then scancel $ARRAY; echo "EMProt setup failed; array cancelled"; break; fi
  [ -z "$A" ] && break; sleep 300
done
echo "== jobs finished $(date)"
# summary of every run record
python3 - > $B/RESULTS_SUMMARY.txt <<'PY'
import glob, json, os
H = os.path.expanduser("~")
rows = []
for f in sorted(glob.glob(f"{H}/bench/out/*/*/run_record.json")) + sorted(glob.glob(f"{H}/mica_bench/outputs/*/run_record.json")):
    try: r = json.load(open(f))
    except Exception as e: rows.append(f"{f}: unreadable ({e})"); continue
    rows.append(f"{r.get('method','MICA'):14} {r.get('case'):16} exit={r.get('exit')} atoms={r.get('final_atoms')} s={r.get('seconds')}")
for f in sorted(glob.glob(f"{H}/mica_bench/outputs/*/refine_record.json")):
    r = json.load(open(f)); rows.append(f"MICA-refine    {r.get('case'):16} exit={r.get('refine_exit')} atoms={r.get('refined_atoms')} s={r.get('seconds')}")
print("\n".join(rows))
PY
cat $B/RESULTS_SUMMARY.txt
cd $HOME && tar -czf $B/RESULTS_$STAMP.tar.gz --exclude="*/see_alpha_output" --exclude="*/cf_*" \
  bench/out bench/logs bench/RESULTS_SUMMARY.txt bench/chain_ids.txt bench/maps/integrity_report_*.json \
  mica_bench/outputs mica_bench/logs mica_bench/docked 2>/dev/null
N_OUT=$(find bench/out mica_bench/outputs -name "*.cif" -o -name "*.pdb" | wc -l)
N_TAR=$(tar -tzf $B/RESULTS_$STAMP.tar.gz | grep -cE "\.(cif|pdb)$")
echo "bundle RESULTS_$STAMP.tar.gz: $(du -h $B/RESULTS_$STAMP.tar.gz | cut -f1), structure files on disk=$N_OUT in bundle=$N_TAR"
if [ "$N_TAR" -ge "$N_OUT" ] && [ "$N_OUT" -gt 0 ]; then
  # optional off-cluster copy: only if the user configured an rclone remote named "gdrive" themselves
  if $HOME/.local/bin/rclone listremotes 2>/dev/null | grep -qx "gdrive:"; then
    $HOME/.local/bin/rclone copy $B/RESULTS_$STAMP.tar.gz gdrive:MICA_benchmark_results/ && echo "uploaded bundle to Google Drive (gdrive:MICA_benchmark_results)" || echo "Google Drive upload failed (bundle kept on GAIVI)"
  else echo "no rclone remote 'gdrive' configured; bundle stays on GAIVI"; fi
  EMPROT_OK=$(grep -c "^EMProt .*exit=0 atoms=[1-9]" $B/RESULTS_SUMMARY.txt)
  echo "EMProt non-empty models: $EMPROT_OK/12"
  # re-creatable bulk only; models, records, logs, maps and docked inputs are kept
  rm -rf /home/a/aravindav/bench/envs/ca /home/a/aravindav/bench/envs/CryoAtom /home/a/aravindav/bench/envs/EModelX \
         /home/a/aravindav/bench/weights /home/a/aravindav/bench/tmp /home/a/aravindav/bench/raw /home/a/aravindav/bench/conda_pkgs \
         /home/a/aravindav/mica_bench/jobtmp /home/a/aravindav/mica_bench/work
  if [ "$EMPROT_OK" -eq 12 ]; then rm -rf /home/a/aravindav/bench/envs/emprot_env /home/a/aravindav/bench/envs/EMProt; else echo "kept EMProt env for debugging"; fi
  echo "cleanup done; home now: $(du -sh $HOME 2>/dev/null | cut -f1)"
else
  echo "BUNDLE CHECK FAILED - nothing deleted"
fi
echo FINALIZE_DONE
