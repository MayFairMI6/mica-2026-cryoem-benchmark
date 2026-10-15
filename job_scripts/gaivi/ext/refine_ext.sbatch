#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -J mica_rsr
#SBATCH -o /home/a/aravindav/mica_bench/logs/refine_%j_%x.out
# Completes MICA's final step (see core-run note): the exact phenix.real_space_refine command run.py launches in the
# background, run to completion (same Phenix 1.20.1-4487, modules/phenix.eff, model, map, resolution).
EMD=$1; B=$HOME/mica_bench; W=$B/work
read RES CASE < <($W/venv/bin/python -c "import json;m=json.load(open('$B/inputs/manifest.json'))['$EMD'];print(m['resolution_A'],m['case'])")
O=$B/outputs/$EMD; M=$(ls $O/*_MICA_all_atom_model.pdb 2>/dev/null | head -1)
if [ -z "$M" ] || ! grep -q "^ATOM" "$M"; then echo "{\"case\":\"$CASE\",\"refine\":\"skipped_no_model\"}" > $O/refine_record.json; exit 0; fi
R=$B/refine/$EMD; mkdir -p $R/out; tar -xzf $B/docked/mica_input_$EMD.tar.gz -C $R $EMD/emd_$EMD.map
cp "$M" $R/out/; cd $R/out; source $W/phenix/phenix-1.20.1-4487/phenix_env.sh
T0=$(date +%s)
phenix.real_space_refine "$R/out/$(basename "$M")" $R/$EMD/emd_$EMD.map $W/MICA/modules/phenix.eff resolution=$RES > $R/out/rsr.log 2>&1
RC=$?; T1=$(date +%s)
F=$(ls $R/out/*_real_space_refined_000.pdb 2>/dev/null | head -1); [ -n "$F" ] && cp "$F" $O/; cp $R/out/rsr.log $O/rsr.log
A=0; [ -n "$F" ] && A=$(grep -cE "^(ATOM|HETATM)" "$F")
echo "{\"case\":\"$CASE\",\"emd\":\"$EMD\",\"refine_exit\":$RC,\"seconds\":$((T1-T0)),\"refined_model\":\"$(basename "$F")\",\"refined_atoms\":$A}" | tee $O/refine_record.json
rm -f $R/$EMD/emd_$EMD.map $R/out/*.pdb $R/out/*.cif $R/out/*.eff $R/out/*.geo 2>/dev/null
