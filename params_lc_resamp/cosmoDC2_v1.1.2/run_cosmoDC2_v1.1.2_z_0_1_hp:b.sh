#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=1.1.2_0_1
#SBATCH -o params_lc_resamp/protoDC2_v1.1.2/logs/v1.1.2_z_0_1_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v1.1.2/logs/v1.1.2_z_0_1_%j.err



./lc_resample.py params_lc_resamp/cosmoDC2_v1.1.2/cosmoDC2_v1.1.2_z_0_1_hp:b.param
