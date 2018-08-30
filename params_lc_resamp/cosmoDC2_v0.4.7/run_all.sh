mv logs/*.err logs/old/.
mv logs/*.out logs/old/.
qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v0.4.7/logs/0_1.out -e params_lc_resamp/cosmoDC2_v0.4.7/logs/0_1.err  --debuglog=params_lc_resamp/cosmoDC2_v0.4.7/logs/0_1.colbalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.7/cosmoDC2_v0.4.7_z_0_1.param
