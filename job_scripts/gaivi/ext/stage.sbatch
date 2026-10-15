#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 23:30:00
#SBATCH -J stage
#SBATCH -o /home/a/aravindav/bench/logs/stage_%j_%x.out
# One method per job on the extension cases: install (same pinned versions as the core run) -> run every case that has
# no recorded outcome yet -> uninstall -> submit the next stage. Keeps one method's software on disk at a time.
# Usage: sbatch -J st_<METHOD> stage.sbatch <METHOD> [NEXT_METHOD ...]     METHOD in: af2 emx ma ca emprot
# After 20 h the job resubmits itself (finished cases are skipped) instead of running into the 23.5 h limit.
M=$1; shift; NEXT="$@"; B=$HOME/bench; E=$B/envs; mkdir -p $E $B/weights $B/tmp; export PATH=$HOME/.local/bin:$PATH TMPDIR=$B/tmp
export SSL_CERT_FILE=/etc/pki/tls/certs/ca-bundle.crt OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
START=$(date +%s); CASES=$(cat $B/${CASEFILE:-cases_ext.txt}); GPU=${CUDA_VISIBLE_DEVICES%%,*}
echo "stage=$M next=[$NEXT] host=$(hostname) gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader -i $GPU) start=$(date -Is)"
case $M in af2) NAME=af2_templates;; emx) NAME="EModelX(+AF)";; ma) NAME=ModelAngelo;; ca) NAME=CryoAtom;; emprot) NAME=EMProt;; *) echo "unknown stage $M"; exit 2;; esac
install() { set -e; cd $E; case $M in
  af2) [ -d cf ] || uv venv -q -p 3.11 cf; uv pip install -q -p cf/bin/python "colabfold[alphafold]==1.6.3" "jax[cuda12]>=0.6.2,<0.12"
       cf/bin/python -c "import jax;assert jax.default_backend()=='gpu';print('CF jax',jax.__version__,jax.devices())";;
  emx) [ -d emx ] || uv venv -q -p 3.8 emx
       uv pip install -q -p emx/bin/python torch==1.8.1+cu111 numpy==1.19.1 scipy==1.5.2 pandas==1.1.3 biopython==1.78 mrcfile==1.3.0 open3d==0.13.0 superpose3d==1.1.1 tqdm==4.61.2 scikit-image==0.18.2 networkx==2.6.3 argparse==1.4.0 requests gdown --find-links https://download.pytorch.org/whl/torch_stable.html
       [ -d EModelX ] || { git clone -q https://github.com/biomed-AI/EModelX.git; git -C EModelX checkout -q 53713cf820d10d78b1ca6669ef197850cf907843; }
       [ -f EModelX/models/best_model_ckpts_done ] || { emx/bin/gdown -q 13BKzEBfL0uubYgcJJTGSZ-PkO9hXAB9X -O emx_models.zip; mkdir -p EModelX/models emx_unzip; unzip -o -q emx_models.zip -d emx_unzip; find emx_unzip -name 'best_*_model.ckpt' -exec cp {} EModelX/models/ \; ; touch EModelX/models/best_model_ckpts_done; rm -rf $E/emx_models.zip $E/emx_unzip; }
       chmod 755 EModelX/modules/pulchra304/src/pulchra; emx/bin/python -c "import torch;assert torch.cuda.is_available();print('EMX torch',torch.__version__)";;
  ma)  [ -d ma ] || uv venv -q -p 3.11 ma; uv pip install -q -p ma/bin/python torch==2.9.1 torchvision
       uv pip install -q -p ma/bin/python git+https://github.com/3dem/model-angelo@994945bdfa6e5368e0d62349a47792f4864eebc3
       TORCH_HOME=$B/weights/modelangelo ma/bin/model_angelo setup_weights --bundle-name nucleotides
       ma/bin/python -c "import torch;assert torch.cuda.is_available();print('MA torch',torch.__version__)";;
  ca)  [ -d CryoAtom ] || { git clone -q https://github.com/YangLab-SDU/CryoAtom.git; git -C CryoAtom checkout -q v1.0.0; }; git -C CryoAtom rev-parse HEAD
       [ -d ca ] || uv venv -q -p 3.8 ca
       uv pip install -q -p ca/bin/python torch==1.11.0+cu113 torchvision==0.12.0+cu113 torchaudio==0.11.0+cu113 numpy==1.21.5 tqdm einops==0.6.0 scipy==1.7.3 biopython==1.81 matplotlib==3.5.3 mrcfile==1.4.3 pandas==1.3.5 fair-esm==2.0.0 pyhmmer==0.7.1 openpyxl setuptools --extra-index-url https://download.pytorch.org/whl/cu113
       [ -f CryoAtom/CryoAtom/checkpoint/CryNet.pth ] || ( cd CryoAtom && wget -q --no-check-certificate https://yanglab.qd.sdu.edu.cn/CryoAtom/download/checkpoints_v1.0.zip && unzip -o -q checkpoints_v1.0.zip && rm -f checkpoints_v1.0.zip )
       ( cd CryoAtom && ../ca/bin/python setup.py -q install ); ca/bin/python -c "import torch,esm,pyhmmer;assert torch.cuda.is_available();print('CA torch',torch.__version__)";;
  emprot) [ -d EMProt ] || { git clone -q https://github.com/huang-laboratory/EMProt.git; git -C EMProt checkout -q 9489d52e7fdc81c85ba82e386abcfa05bac2b9d1; }
       source /apps/anaconda3/etc/profile.d/conda.sh; export CONDA_PKGS_DIRS=$B/conda_pkgs PIP_NO_CACHE_DIR=1
       [ -d $E/emprot_env ] || conda env create -q -p $E/emprot_env -f EMProt/environment.yml
       conda activate $E/emprot_env; ( cd EMProt && pip install -q -e . && chmod +x emprot/bin/* )
       ( cd EMProt/emprot && { [ -d weights ] || { curl -sfL -o weights.tgz http://huanglab.phys.hust.edu.cn/EMProt/pretrained_weights/weights.tgz && tar -zxf weights.tgz && rm -f weights.tgz; }; } )
       rm -rf /home/a/aravindav/bench/conda_pkgs; python -c "import torch;assert torch.cuda.is_available();print('EMProt torch',torch.__version__)"; conda deactivate;;
  esac; uv cache clean -q 2>/dev/null; set +e; }
uninstall() { case $M in
  af2) rm -rf /home/a/aravindav/bench/envs/cf /home/a/aravindav/bench/weights/cf_cache;;
  emx) rm -rf /home/a/aravindav/bench/envs/emx /home/a/aravindav/bench/envs/EModelX;;
  ma)  rm -rf /home/a/aravindav/bench/envs/ma /home/a/aravindav/bench/weights/modelangelo;;
  ca)  rm -rf /home/a/aravindav/bench/envs/ca /home/a/aravindav/bench/envs/CryoAtom /home/a/aravindav/bench/weights/torch_hub;;
  emprot) rm -rf /home/a/aravindav/bench/envs/emprot_env /home/a/aravindav/bench/envs/EMProt;; esac; }
record() { # exit model seconds
  P=$(sort -n $OUT/.gpu 2>/dev/null | tail -1); rm -f $OUT/.gpu; A=0; [ -f "$2" ] && A=$(grep -cE "^(ATOM|HETATM)" "$2")
  echo "{\"method\":\"$NAME\",\"case\":\"$CASE\",\"exit\":$1,\"seconds\":$3,\"peak_gpu_mib\":${P:-0},\"final_model\":\"$(basename "$2")\",\"final_atoms\":$A,\"host\":\"$(hostname)\",\"gpu\":\"$(nvidia-smi --query-gpu=name --format=csv,noheader -i $GPU)\"}" | tee $OUT/run_record.json; }
run_case() { PDB=${CASE##*_}; FASTA=$B/frozen_inputs/$CASE/$PDB.fasta; MAP=$B/maps/$CASE/canonical.mrc; CUBE=$B/maps/$CASE/canonical_cube.mrc; OUT=$B/out/$NAME/$CASE; mkdir -p "$OUT"
  [ -f $MAP ] || { echo "map missing for $CASE"; return 9; }
  ( while :; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i $GPU >> "$OUT/.gpu" 2>/dev/null; sleep 5; done ) & MON=$!; T0=$(date +%s)
  case $M in
  af2) export XDG_CACHE_HOME=$B/weights/cf_cache TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0; mkdir -p $OUT/in; RC=0
       $E/cf/bin/python - $FASTA $OUT/in <<'PY'
import sys
for b in open(sys.argv[1]).read().strip().split(">")[1:]:
    h, *s = b.strip().split("\n"); e = h.split("|")[0]; open(f"{sys.argv[2]}/{e}.fasta", "w").write(f">{e}\n{''.join(s).strip()}\n")
PY
       for f in $OUT/in/*.fasta; do e=$(basename $f .fasta); [ -f $OUT/$e/ranked_0.pdb ] && continue
         $E/cf/bin/colabfold_batch --num-models 5 --num-recycle 3 --random-seed 0 $f $OUT/cf_$e || RC=1
         P=$(ls $OUT/cf_$e/${e}*_rank_001_*.pdb 2>/dev/null | head -1); [ -n "$P" ] && mkdir -p $OUT/$e && cp $P $OUT/$e/ranked_0.pdb && cp $OUT/cf_$e/${e}*scores_rank_001*.json $OUT/$e/ 2>/dev/null
         [ -n "$P" ] && rm -rf "$OUT/cf_$e"; done
       $E/cf/bin/python - $OUT <<'PY'
import glob, json, os, sys
three = dict(ALA="A",ARG="R",ASN="N",ASP="D",CYS="C",GLN="Q",GLU="E",GLY="G",HIS="H",ILE="I",LEU="L",LYS="K",MET="M",PHE="F",PRO="P",SER="S",THR="T",TRP="W",TYR="Y",VAL="V")
o = sys.argv[1]; rec = {}
for f in glob.glob(f"{o}/in/*.fasta"):
    e = os.path.basename(f)[:-6]; seq = open(f).read().split("\n", 1)[1].replace("\n", ""); p = f"{o}/{e}/ranked_0.pdb"
    if not os.path.exists(p): rec[e] = {"length": len(seq), "template": None}; continue
    ca = [l for l in open(p) if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    rec[e] = {"length": len(seq), "sequence_matches_fasta": "".join(three.get(l[17:20], "X") for l in ca) == seq, "mean_plddt": round(sum(float(l[60:66]) for l in ca) / len(ca), 1)}
json.dump(rec, open(f"{o}/templates_record.json", "w"), indent=1); print(rec)
PY
       kill $MON; rm -f $OUT/.gpu; echo "{\"method\":\"af2_templates\",\"case\":\"$CASE\",\"exit\":$RC,\"seconds\":$(( $(date +%s)-T0 ))}" | tee $OUT/run_record.json; return;;
  emx) T=$B/out/af2_templates/$CASE; for e in $(grep ">" $FASTA | sed 's/>//; s/|.*//'); do [ -f $T/$e/ranked_0.pdb ] || { kill $MON; echo "missing template $e"; record 2 none 0; return; }; done
       mkdir -p "$OUT/pred"; ( cd $E/EModelX && $E/emx/bin/python run.py --protocol=temp_flex --EM_map $MAP --fasta $FASTA --template_dir $T --output_dir "$OUT/pred" --seed 2022 --mul_proc_num 4 --run_pulchra --pulchra_path modules/pulchra304/src/pulchra ); RC=$?
       F=$(ls "$OUT"/pred/*all_atom*.pdb 2>/dev/null | head -1); [ -n "$F" ] && cp "$F" "$OUT/"; kill $MON; record $RC "$OUT/$(basename "$F")" $(( $(date +%s)-T0 ));;
  ma)  TORCH_HOME=$B/weights/modelangelo $E/ma/bin/model_angelo build -v $CUBE -pf $FASTA -o $OUT/$CASE --device 0; RC=$?
       cp $OUT/$CASE/$CASE.cif $OUT/$CASE/${CASE}_raw.cif $OUT/ 2>/dev/null; kill $MON; record $RC $OUT/$CASE.cif $(( $(date +%s)-T0 ));;
  ca)  ( cd $OUT && TORCH_HOME=$B/weights/torch_hub $E/ca/bin/cryoatom build -s $FASTA -v $MAP -o $OUT/$CASE -d cuda:0 ); RC=$?
       cp $OUT/$CASE/$CASE.cif $OUT/$CASE/${CASE}_raw.cif $OUT/ 2>/dev/null; kill $MON; record $RC $OUT/$CASE.cif $(( $(date +%s)-T0 ));;
  emprot) source /apps/anaconda3/etc/profile.d/conda.sh; conda activate $E/emprot_env
       LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH emprot build --map $MAP --seq $FASTA --device 0 --output $OUT/run > $OUT/emprot.log 2>&1; RC=$?; tail -n 5 $OUT/emprot.log; conda deactivate
       cp $OUT/run/output.cif $OUT/$CASE.cif 2>/dev/null; cp $OUT/run/output_denovo.cif $OUT/${CASE}_denovo.cif 2>/dev/null; rm -rf "$OUT/run/temp"
       kill $MON; record $RC $OUT/$CASE.cif $(( $(date +%s)-T0 ));;
  esac; }
( install ) || { echo "INSTALL FAILED for $M — stage stopped (later stages not started)"; exit 1; }
LEFT=0
for CASE in $CASES; do
  [ -f "$B/out/$NAME/$CASE/run_record.json" ] && { echo "skip $CASE (recorded)"; continue; }
  if [ $(( $(date +%s)-START )) -gt 72000 ]; then LEFT=1; break; fi
  echo "=== $M $CASE $(date +%T)"; run_case
done
cd $B
if [ $LEFT -eq 1 ]; then echo "time budget reached; resubmitting stage $M"; sbatch -J st_$M $B/stage.sbatch $M $NEXT; exit 0; fi
MISSING=0; for CASE in $CASES; do [ -f "$B/out/$NAME/$CASE/run_record.json" ] || MISSING=$((MISSING+1)); done
echo "stage $M finished; cases without a record: $MISSING"
[ $MISSING -eq 0 ] && uninstall
if [ -n "$NEXT" ]; then set -- $NEXT; N1=$1; shift; sbatch -J st_$N1 $B/stage.sbatch $N1 "$@"; fi
echo STAGE_DONE
