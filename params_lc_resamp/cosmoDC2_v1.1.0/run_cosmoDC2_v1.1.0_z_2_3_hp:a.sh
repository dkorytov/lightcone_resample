#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=1.1.0_2_3
#SBATCH -o params_lc_resamp/protoDC2_v1.1.0/logs/v1.1.0_z_2_3_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v1.1.0/logs/v1.1.0_z_2_3_%j.err



./lc_resample.py params_lc_resamp/cosmoDC2_v1.1.0/cosmoDC2_v1.1.0_z_2_3_hp:a.param
