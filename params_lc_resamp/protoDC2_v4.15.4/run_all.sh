mv logs/*.err logs/old/.
mv logs/*.out logs/old/.
sbatch params_lc_resamp/protoDC2_v4.15.4/run_protoDC2_v4.15.4_z_0_1.sh
sbatch params_lc_resamp/protoDC2_v4.15.4/run_protoDC2_v4.15.4_z_1_2.sh
sbatch params_lc_resamp/protoDC2_v4.15.4/run_protoDC2_v4.15.4_z_2_3.sh
