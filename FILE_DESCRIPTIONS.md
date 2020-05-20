File Descriptions
=================

This file provides a quick description of each python file in the
repo. The files are separated into eight generic groups. There are
additional utility scripts for generating configuration files and
cooley submission scripts described in >link<


Main Executable
--------------
The main body of the pipeline. 

[lc_resample.py](lc_resample.py)
* Performs the matchup between the Galacticus library and baseDC2
  catalog.  The script takes in a single configuration files found in
  `param_lc_resample`. A description of the parameter file can be
  found here >link<


Precomputations
---------------

These script precompute quantities used in `lc_resample.py`. Usually
only ran once.

[k_corr_step.py](k_corr_step.py)
* Links library galaxies across an earlier and later time step. The
  links are recorded as indexes to rearrange the later time step to
  match the earlier timestep. Only needs to run once per library.

[precompute_mask.py](precompute_mask.py)
* Precomputes which galaxies pass Eve's Galmatcher color and other
  cuts. Stored as a boolean array for each timestep. 


Function Library & Definitions
----------------------------

These scripts contain definitions and functions used by other scripts.

[ellipticity_model.py](ellipticity_model.py)
* Generates the galaxy ellipticity distribution.

[pecZ.py](pecZ.py)
* Calculates observed redshift from peculiar velocity. Joe's script.

Diagnostics & Plots
-------------------

These scripts plot various validations and sanity plots. Most of these
scripts take in the same config file that `lc_resample.py` used to
produce the catalog.

[plot_color_match.py](plot_color_match.py)

[plot_color_redshift_diagnostics.py](plot_color_redshift_diagnostics.py)

[plot_colors.py](plot_colors.py)

[plot_dust_effect.py](plot_dust_effect.py)

[plot_ellip.py](plot_ellip.py)

[plot_gr_ri.py](plot_gr_ri.py)


One-shot Tests
--------------

Scripts used figure out found bugs & issues. They were not used again
after their purpose was done.

[find_galaxies.py](find_galaxies.py)
* Finds galaxies by hard coded ra/dec halo id. Saves basic halo and
  central galaxy info from those galaxies
  
[find_galaxy.py](find_galaxy.py)
* Finds a galaxy by ra/dec + other info. Was used to investigate a
  strange galaxy with a strange Av/Rv values.

[test_ellip.py](test_ellip.py)

[test_k_corr_step.py](test_k_corr_step.py)

[test_linear_interpolate.py](test_linear_interpolate.py)

[test_nan_totals.py](test_nan_totals.py)

[test_v4.py](test_v4.py)

[ellipticity_model_testing.py](ellipticity_model_testing.py)
* I don't recall.

Quick Fixes 
-----------
[correct_bulge_ellip.py](correct_bulge_ellip.py)
[correct_dec_true.py](correct_dec_true.py)
[move_ra_dec.py](move_ra_dec.py)
[shear_insert.py](shear_insert.py)
[fix_version.py](fix_version.py)
[recalculate_ellipticity.py](recalculate_ellipticity.py)

Deprecated 
----------------------
[lightcone_fuzzy.py](lightcone_fuzzy.py)
[lightcone_matchup.py](lightcone_matchup.py)
[k_corr_step.sl](k_corr_step.sl)
[k_corr_step_append.py](k_corr_step_append.py)
[precompute_interpolation_index.py](precompute_interpolation_index.py)


Parameter Folders
---------------
[params_lc_matchup](params_lc_matchup)

[params_kcorr](params_kcorr)

[params_lc_resamp](params_lc_resamp)












