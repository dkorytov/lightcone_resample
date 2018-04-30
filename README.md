There is one main script to the galacticus matchup: lc_resample.py (light cone resample)
There are several precompute scritps:
All scirpts run with one argrument: their parameter file. 
i.e: "./lc_resample.py param_lc_resamp/v4.1.dust136.param"
      
1) lightcone_matchup.py param_lc_match/*.param:
  This script takes in parameters from param_lc_matchup. 
  This code takes in the binary outputs from the lightcone cutout cod and the
  UMachine + Color snapshot catalog and creates a lightcone catalog with all 
  snapshot catalog properties. If the positions don't change, then all one has
  to do is change the "snapshot_galaxy_fname" parameter to point to the new 
  UMachine + color snapshot catalog. This is a fast script. 

2) k_cor_step.py param_kcorr/*.param: 
   This scirpt generates the index to relate galaxies from step B to step A, where
    step B is later in time. 

3) precompute_mask.py param_lc_resamp/*.param:
   This script precomputes the masks on snaps shot (eg: color cuts, mstar cuts, etc)
   It should run on teh same parameter file as lc_resamp

=====================
4) Main Script:
     ./lc_resample.py param_lc_resampe/*.param
     This function matches up the lightcone UMachine galaxies to Galacticus galaxies. 

