#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=4.15.3_0_1
#SBATCH -o params_lc_resamp/protoDC2_v4.15.3/logs/v4.15.3_z_0_1_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v4.15.3/logs/v4.15.3_z_0_1_%j.err



./lc_resample.py params_lc_resamp/protoDC2_v4.15.3/protoDC2_v4.15.3_z_0_1.param
