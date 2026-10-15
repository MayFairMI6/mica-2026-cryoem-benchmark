#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH -J bench_maps
#SBATCH -o /home/a/aravindav/bench/logs/prep_maps_%j.out
cd $HOME/bench && $HOME/mica_bench/work/venv/bin/python prep_maps_gaivi.py "$@"
