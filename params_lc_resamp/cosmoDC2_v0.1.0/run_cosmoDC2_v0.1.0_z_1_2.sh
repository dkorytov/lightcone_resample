#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=0.1.0_1_2
#SBATCH -o params_lc_resamp/protoDC2_v0.1.0/logs/v0.1.0_z_1_2_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v0.1.0/logs/v0.1.0_z_1_2_%j.err



./lc_resample.py params_lc_resamp/protoDC2_v0.1.0/protoDC2_v0.1.0_z_1_2.param
