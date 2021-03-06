#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=1.1.1_2_3
#SBATCH -o params_lc_resamp/protoDC2_v1.1.1/logs/v1.1.1_z_2_3_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v1.1.1/logs/v1.1.1_z_2_3_%j.err



./lc_resample.py params_lc_resamp/cosmoDC2_v1.1.1/cosmoDC2_v1.1.1_z_2_3_hp:b.param
