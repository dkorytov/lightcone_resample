#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=0.4.8_1_2
#SBATCH -o params_lc_resamp/protoDC2_v0.4.8/logs/v0.4.8_z_1_2_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v0.4.8/logs/v0.4.8_z_1_2_%j.err



./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.8/cosmoDC2_v0.4.8_z_1_2_hp:a.param
