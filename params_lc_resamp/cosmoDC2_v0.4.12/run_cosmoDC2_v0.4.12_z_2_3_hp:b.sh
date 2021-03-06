#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=0.4.12_2_3
#SBATCH -o params_lc_resamp/protoDC2_v0.4.12/logs/v0.4.12_z_2_3_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v0.4.12/logs/v0.4.12_z_2_3_%j.err



./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.12/cosmoDC2_v0.4.12_z_2_3_hp:b.param
