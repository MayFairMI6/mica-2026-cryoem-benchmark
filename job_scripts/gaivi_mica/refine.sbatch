#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH -c 4
#SBATCH -t 08:00:00
#SBATCH -J mica_rsr
#SBATCH -o /home/a/aravindav/mica_bench/logs/refine_%A_%a.out
# Completes MICA's final step. run.py calls modules/phenix.sh, which runs
#   phenix.real_space_refine <MICA_all_atom_model.pdb> <map> modules/phenix.eff resolution=<res>
# in the output dir with a trailing "&"; batch schedulers kill it when run.py exits. This job runs that
# exact command (same Phenix 1.20.1-4487, same phenix.eff, same model/map/resolution) to completion.
B=$HOME/mica_bench; W=$B/work; EMD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $B/cases.txt)
read RES CASE < <($W/venv/bin/python -c "import json;m=json.load(open('$B/inputs/manifest.json'))['$EMD'];print(m['resolution_A'],m['case'])")
O=$B/outputs/$EMD; M=$(ls $O/*_MICA_all_atom_model.pdb 2>/dev/null | head -1)
if [ -z "$M" ] || ! grep -q "^ATOM" "$M"; then echo "{\"case\":\"$CASE\",\"refine\":\"skipped_no_model\"}" > $O/refine_record.json; exit 0; fi
R=$B/refine/$EMD; rm -rf $R; mkdir -p $R/out; tar -xzf $B/docked/mica_input_$EMD.tar.gz -C $R $EMD/emd_$EMD.map
cp "$M" $R/out/; cd $R/out
source $W/phenix/phenix-1.20.1-4487/phenix_env.sh
T0=$(date +%s)
phenix.real_space_refine "$R/out/$(basename "$M")" $R/$EMD/emd_$EMD.map $W/MICA/modules/phenix.eff resolution=$RES > $R/out/rsr.log 2>&1
RC=$?; T1=$(date +%s)
F=$(ls $R/out/*_real_space_refined_000.pdb 2>/dev/null | head -1)
[ -n "$F" ] && cp "$F" $O/; cp $R/out/rsr.log $O/rsr.log
A=0; [ -n "$F" ] && A=$(grep -cE "^(ATOM|HETATM)" "$F")
echo "{\"case\":\"$CASE\",\"emd\":\"$EMD\",\"refine_exit\":$RC,\"seconds\":$((T1-T0)),\"refined_model\":\"$(basename "$F")\",\"refined_atoms\":$A}" | tee $O/refine_record.json
rm -rf $R
