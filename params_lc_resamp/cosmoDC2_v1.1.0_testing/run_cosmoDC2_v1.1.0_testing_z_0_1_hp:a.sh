#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=1.1.0_0_1
#SBATCH -o params_lc_resamp/protoDC2_v1.1.0/logs/v1.1.0_z_0_1_%j.out
#SBATCH -e params_lc_resamp/protoDC2_v1.1.0/logs/v1.1.0_z_0_1_%j.err



./lc_resample.py params_lc_resamp/cosmoDC2_v1.1.0_testing/cosmoDC2_v1.1.0_testing_z_0_1_hp:a.param
