#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=0.4.10_0_1
#SBATCH -o params_lc_resamp/protoDC2_v0.4.10/logs/v0.4.10_z_0_1_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v0.4.10/logs/v0.4.10_z_0_1_%j.err



./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.10/cosmoDC2_v0.4.10_z_0_1_hp:a.param
