#!/bin/bash
mv logs/*.err logs/old/.
mv logs/*.out logs/old/.
qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v0.4.8/logs/0_1_a.out -e params_lc_resamp/cosmoDC2_v0.4.8/logs/0_1_a.err  --debuglog=params_lc_resamp/cosmoDC2_v0.4.8/logs/0_1_a.cobalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.8/cosmoDC2_v0.4.8_z_0_1_hp:a.param
qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v0.4.8/logs/0_1_b.out -e params_lc_resamp/cosmoDC2_v0.4.8/logs/0_1_b.err  --debuglog=params_lc_resamp/cosmoDC2_v0.4.8/logs/0_1_b.cobalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.8/cosmoDC2_v0.4.8_z_0_1_hp:b.param
qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v0.4.8/logs/1_2_a.out -e params_lc_resamp/cosmoDC2_v0.4.8/logs/1_2_a.err  --debuglog=params_lc_resamp/cosmoDC2_v0.4.8/logs/1_2_a.cobalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.8/cosmoDC2_v0.4.8_z_1_2_hp:a.param
qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v0.4.8/logs/1_2_b.out -e params_lc_resamp/cosmoDC2_v0.4.8/logs/1_2_b.err  --debuglog=params_lc_resamp/cosmoDC2_v0.4.8/logs/1_2_b.cobalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.8/cosmoDC2_v0.4.8_z_1_2_hp:b.param
qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v0.4.8/logs/2_3_a.out -e params_lc_resamp/cosmoDC2_v0.4.8/logs/2_3_a.err  --debuglog=params_lc_resamp/cosmoDC2_v0.4.8/logs/2_3_a.cobalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.8/cosmoDC2_v0.4.8_z_2_3_hp:a.param
qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v0.4.8/logs/2_3_b.out -e params_lc_resamp/cosmoDC2_v0.4.8/logs/2_3_b.err  --debuglog=params_lc_resamp/cosmoDC2_v0.4.8/logs/2_3_b.cobalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v0.4.8/cosmoDC2_v0.4.8_z_2_3_hp:b.param
