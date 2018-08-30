#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=0.4.6_0_1
#SBATCH -o params_lc_resamp/protoDC2_v0.4.6/logs/v0.4.6_z_0_1_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v0.4.6/logs/v0.4.6_z_0_1_%j.err



./lc_resample.py params_lc_resamp/protoDC2_v0.4.6/protoDC2_v0.4.6_z_0_1.param
