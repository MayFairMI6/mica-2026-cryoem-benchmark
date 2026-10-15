#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J bench_env
#SBATCH -o /home/a/aravindav/bench/logs/env_setup_%j.out
# Environments for GAIVI runs (same versions as the Kaggle runs; uv venvs instead of conda):
#  ModelAngelo 1.0.18 @994945bd (py3.11, torch 2.9.1) | ColabFold 1.6.3 (py3.11, jax[cuda12])
#  EModelX @53713cf8 (py3.8, torch 1.8.1+cu111)       | CryoAtom v1.0.0 @e2a31f2b (py3.8, torch 1.11.0+cu113, linux.yml pins)
export PATH=$HOME/.local/bin:$PATH; E=$HOME/bench/envs; mkdir -p $E $HOME/bench/weights; cd $E
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
step(){ echo "=== $1 $(date +%T)"; shift; "$@"; echo "=== exit $? $(date +%T)"; }
# ModelAngelo
[ -d ma ] || step ma_venv uv venv -q -p 3.11 ma
step ma_torch uv pip install -q -p ma/bin/python torch==2.9.1 torchvision
step ma_pkg uv pip install -q -p ma/bin/python git+https://github.com/3dem/model-angelo@994945bdfa6e5368e0d62349a47792f4864eebc3
step ma_weights env TORCH_HOME=$HOME/bench/weights/modelangelo ma/bin/model_angelo setup_weights --bundle-name nucleotides
step ma_gpu ma/bin/python -c "import torch;print('MA torch',torch.__version__,torch.version.cuda,torch.cuda.is_available());x=torch.randn(512,512,device='cuda');print((x@x).sum().item())"
# ColabFold
[ -d cf ] || step cf_venv uv venv -q -p 3.11 cf
step cf_pkg uv pip install -q -p cf/bin/python "colabfold[alphafold]==1.6.3" "jax[cuda12]>=0.6.2,<0.12"
step cf_gpu cf/bin/python -c "import jax;print('CF jax',jax.__version__,jax.default_backend(),jax.devices())"
# EModelX
[ -d emx ] || step emx_venv uv venv -q -p 3.8 emx
step emx_pkg uv pip install -q -p emx/bin/python torch==1.8.1+cu111 numpy==1.19.1 scipy==1.5.2 pandas==1.1.3 biopython==1.78 mrcfile==1.3.0 open3d==0.13.0 superpose3d==1.1.1 tqdm==4.61.2 scikit-image==0.18.2 networkx==2.6.3 argparse==1.4.0 requests gdown --find-links https://download.pytorch.org/whl/torch_stable.html
if [ ! -d EModelX ]; then git clone -q https://github.com/biomed-AI/EModelX.git && git -C EModelX checkout -q 53713cf820d10d78b1ca6669ef197850cf907843; fi
if [ ! -f EModelX/models/best_model_ckpts_done ]; then
  step emx_weights bash -c "emx/bin/gdown -q 13BKzEBfL0uubYgcJJTGSZ-PkO9hXAB9X -O emx_models.zip && mkdir -p EModelX/models emx_unzip && unzip -o -q emx_models.zip -d emx_unzip && find emx_unzip -name 'best_*_model.ckpt' -exec cp {} EModelX/models/ \; && ls EModelX/models && touch EModelX/models/best_model_ckpts_done && rm -rf emx_models.zip emx_unzip"
fi
chmod 755 EModelX/modules/pulchra304/src/pulchra
step emx_gpu emx/bin/python -c "import torch;print('EMX torch',torch.__version__,torch.cuda.is_available(),torch.cuda.get_device_name(0))"
# CryoAtom
if [ ! -d CryoAtom ]; then git clone -q https://github.com/YangLab-SDU/CryoAtom.git && git -C CryoAtom checkout -q v1.0.0; fi
git -C CryoAtom rev-parse HEAD
[ -d ca ] || step ca_venv uv venv -q -p 3.8 ca
step ca_pkg uv pip install -q -p ca/bin/python torch==1.11.0+cu113 torchvision==0.12.0+cu113 torchaudio==0.11.0+cu113 numpy==1.21.5 tqdm einops==0.6.0 scipy==1.7.3 biopython==1.81 matplotlib==3.5.3 mrcfile==1.4.3 pandas==1.3.5 fair-esm==2.0.0 pyhmmer==0.7.1 openpyxl setuptools --extra-index-url https://download.pytorch.org/whl/cu113
if [ ! -f CryoAtom/CryoAtom/checkpoint/CryNet.pth ]; then
  step ca_weights bash -c "cd CryoAtom && wget -q --no-check-certificate https://yanglab.qd.sdu.edu.cn/CryoAtom/download/checkpoints_v1.0.zip && unzip -o -q checkpoints_v1.0.zip && rm -f checkpoints_v1.0.zip && ls -la CryoAtom/checkpoint"
fi
step ca_install bash -c "cd CryoAtom && ../ca/bin/python setup.py -q install"
step ca_help ca/bin/cryoatom build -h
step ca_gpu ca/bin/python -c "import torch,esm,pyhmmer;print('CA torch',torch.__version__,torch.cuda.is_available())"
echo ENV_SETUP_DONE
