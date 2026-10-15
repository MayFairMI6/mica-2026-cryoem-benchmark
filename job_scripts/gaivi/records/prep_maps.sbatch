#!/bin/bash
#SBATCH -p ClsGenAI
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH -J bench_maps
#SBATCH -o /home/a/aravindav/bench/logs/prep_maps_%j.out
cd $HOME/bench && $HOME/bench/envs/maps/bin/python prep_maps_gaivi.py "$@"
