#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 03:00:00
#SBATCH -J mica_setup
#SBATCH -o /home/a/aravindav/mica_bench/logs/setup_%j.out
# MICA @ 726d1bb1 env (same pins as Kaggle runs), Phenix 1.20.1-4487 (private, licensee-only), trained models.
set -x
B=$HOME/mica_bench; W=$B/work; mkdir -p $W; cd $W
hostname; nvidia-smi; nproc; free -g | head -2; curl -sI https://pypi.org | head -1
export PATH=$HOME/.local/bin:$PATH
command -v uv || curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH=$HOME/.local/bin:$PATH
[ -d venv ] || uv venv -q -p 3.8 venv
uv pip install -q -p venv/bin/python torch==1.8.1+cu111 numpy==1.19.1 scipy==1.5.2 pandas==1.1.3 biopython==1.78 mrcfile==1.3.0 open3d==0.13.0 superpose3d==1.1.1 tqdm==4.61.2 scikit-image==0.18.2 networkx==2.6.3 psutil einops rotary-embedding-torch natsort wandb matplotlib scikit-learn requests --find-links https://download.pytorch.org/whl/torch_stable.html
[ -d MICA ] || { git clone -q https://github.com/jianlin-cheng/MICA.git MICA && git -C MICA checkout -q 726d1bb14be04887a43b7f5aefcff0d843761d01; }
cd MICA
grep -q "batch_threshold = 10\*\*9" utils/predict.py || sed -i "s/self.batch_threshold = 200/self.batch_threshold = 10**9  # per-sample path (batch 1); result-neutral in eval\/no_grad/" utils/predict.py
grep -n batch_threshold utils/predict.py | head -2
[ -f trained_models/MICA_best_model.pth ] || { curl -sfL "https://zenodo.org/records/15756654/files/trained_models.tar.gz?download=1" -o tm.tgz && tar -xzf tm.tgz && rm tm.tgz; }
chmod 755 modules/pulchra304/src/pulchra
cd $W
if [ ! -f phenix/phenix-1.20.1-4487/phenix_env.sh ]; then
  mkdir -p phx_src && tar -xzf $B/phenix/phenix-installer-1.20.1-4487-*.tar.gz -C phx_src
  (cd phx_src/phenix-installer-1.20.1-4487-* && ./install --prefix=$W/phenix --nproc=4) && rm -rf phx_src
fi
ls $W/phenix/phenix-1.20.1-4487/phenix_env.sh
$W/venv/bin/python -c "import torch;print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0))"
echo SETUP_OK
