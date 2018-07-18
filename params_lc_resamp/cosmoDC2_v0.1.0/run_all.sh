mv logs/*.err logs/old/.
mv logs/*.out logs/old/.
sbatch params_lc_resamp/cosmoDC2_v0.1.0/run_cosmoDC2_v0.1.0_z_0_1.sh
sbatch params_lc_resamp/cosmoDC2_v0.1.0/run_cosmoDC2_v0.1.0_z_1_2.sh
sbatch params_lc_resamp/cosmoDC2_v0.1.0/run_cosmoDC2_v0.1.0_z_2_3.sh
