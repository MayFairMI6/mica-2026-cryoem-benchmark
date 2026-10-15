#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 03:00:00
#SBATCH -J emprot_env
#SBATCH -o /home/a/aravindav/bench/logs/emprot_setup_%j.out
# EMProt v1.2 @9489d52e (README install: conda env from environment.yml, pip install -e ., chmod bin, weights.tgz).
E=$HOME/bench/envs; cd $E
[ -d EMProt ] || { git clone -q https://github.com/huang-laboratory/EMProt.git && git -C EMProt checkout -q 9489d52e7fdc81c85ba82e386abcfa05bac2b9d1; }
git -C EMProt rev-parse HEAD
source /apps/anaconda3/etc/profile.d/conda.sh
export CONDA_PKGS_DIRS=$HOME/bench/conda_pkgs
[ -d $E/emprot_env ] || conda env create -q -p $E/emprot_env -f EMProt/environment.yml
conda activate $E/emprot_env
cd EMProt && pip install -q -e . && chmod +x emprot/bin/*
cd emprot && [ -d weights ] || { curl -sfL -o weights.tgz http://huanglab.phys.hust.edu.cn/EMProt/pretrained_weights/weights.tgz && tar -zxf weights.tgz && rm -f weights.tgz; }
ls; cd $E
rm -rf $HOME/bench/conda_pkgs
emprot --version; emprot build --help | head -30
python -c "import torch;print('EMProt torch',torch.__version__,torch.cuda.is_available(),torch.cuda.get_device_name(0))"
echo EMPROT_SETUP_DONE
