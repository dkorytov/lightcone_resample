#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=4.15.2_1_2
#SBATCH -o params_lc_resamp/protoDC2_v4.15.2/logs/v4.15.2_z_1_2_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v4.15.2/logs/v4.15.2_z_1_2_%j.err



./lc_resample.py params_lc_resamp/protoDC2_v4.15.2/protoDC2_v4.15.2_z_1_2.param
