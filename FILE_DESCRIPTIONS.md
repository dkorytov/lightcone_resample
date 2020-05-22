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
  found in [PARAM_DESCRIPTIONS.md](PARAM_DESCRIPTIONS.md)


Precomputations
---------------

These script precompute quantities used in `lc_resample.py`. Usually
only ran once.

[k_corr_step.py](k_corr_step.py)
* Links library galaxies across an earlier and later time step. The
  links are recorded as indexes to rearrange the later time step to
  match the earlier timestep. Only needs to run once per library. The
  parameter file for this script can be found in
  [PARAM_DESCRIPTIONS.md](PARAM_DESCRIPTIONS.md).

[precompute_mask.py](precompute_mask.py)
* Precomputes which galaxies pass Eve's Galmatcher color and other
  cuts. Stored as a boolean array for each timestep. This script takes
  in a `lc_resample.py` parameter file, overwriting the file that
  `lc_resample.py` expects to load.


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
* Plots the color match between baseDC2 and library galaxies for
  step323.

[plot_color_redshift_diagnostics.py](plot_color_redshift_diagnostics.py)
* Plots color vs redshift for a sanity checks and color bands. This
  was adapted into DESCQA test called `color_z` (I think)

[plot_colors.py](plot_colors.py)
* Plots many different color diagnostics plots. 

[plot_dust_effect.py](plot_dust_effect.py)
* Plots the strength of reddening from dust

[plot_ellip.py](plot_ellip.py)
* Plots the distribution of galaxy ellipticity. 

[plot_gr_ri.py](plot_gr_ri.py)
* Plots the gr-ri color-color for protoDC2

[test_ellip.py](test_ellip.py)
* Plots many ellipticity diagnostics. Also plots distribution of
  ellipticities for different B/T ratio galaxy bins.

One-shot Tests
--------------

Scripts used figure out found bugs & issues. They were not used again
after their purpose was completed.

[find_galaxies.py](find_galaxies.py)
* Finds galaxies by hard coded ra/dec halo id. Saves basic halo and
  central galaxy info from those galaxies
  
[find_galaxy.py](find_galaxy.py)
* Finds a galaxy by ra/dec + other info. Was used to investigate a
  strange galaxy with a strange Av/Rv values.

[test_linear_interpolate.py](test_linear_interpolate.py)
* Tests my implementation for linear interpolating Galacticus galaxy
properties between time steps. 

[test_nan_totals.py](test_nan_totals.py)
* Gives number of nan in various columns.

[test_v4.py](test_v4.py)
* Many sanity plots for protoDC2 v4.X. 

[ellipticity_model_testing.py](ellipticity_model_testing.py)
* Testing imports.

Quick Fixes 
-----------

These script were used to fix something directly in the catalog. Most
of these scripts take a `lc_resample.py` parameter file

[correct_bulge_ellip.py](correct_bulge_ellip.py)
* Recalculates ellipticity from major and minor galaxy axis lengths
  and overwrites the new values in the catalog hdf5 files. Takes in a
  `lc_resample.py` param file.
  
[correct_dec_true.py](correct_dec_true.py)
* Subtracts 85 deg from dec_true and writes new value into hdf5
  catalog files. Takes in an `lc_resample.py` param file.

[move_ra_dec.py](move_ra_dec.py) 
* Swaps ra and dec in hdf5 files. Takes in an `lc_resample.py` param file.

[shear_insert.py](shear_insert.py)
* Overwrites stored shear values with new values. The new shear values
  are taken from the shear file specified in the given
  `lc_resample.py` param file.

[fix_version.py](fix_version.py)
* Rewrites the current catalog version number in the hdf5 files and
  overwrites them with the version number in the param file.

[recalculate_ellipticity.py](recalculate_ellipticity.py)
* Recalculates ellipticity using the methods defined in
  `lc_resample.py` and overwrites new values into the catalog hdf5
  files. Takes in a `lc_resample.py` param file.
  
Deprecated 
----------------------
[lightcone_fuzzy.py](lightcone_fuzzy.py)
* Shuffles some baseDC2 galaxies into earlier and later time steps in
  attempt to solve color-z banding. It's superseded by galaxy
  interpolation. 

[lightcone_matchup.py](lightcone_matchup.py)
* Combines baseDC2 lightcone with calculated shears. The process is now 
done in `lc_resample.py`

[k_corr_step.sl](k_corr_step.sl)
* submission script for k_corr_step--not needed.

[k_corr_step_append.py](k_corr_step_append.py)
* Not sure what the functional difference is with `k_corr_step.py`. I
  believe it resumes from `k_corr_step.py` run instead of starting
  from 0. 
* TODO

[test_k_corr_step.py](test_k_corr_step.py)
* makes sure that `k_corr_step.py` computed slopes for every column

[precompute_interpolation_index.py](precompute_interpolation_index.py)
* Never implemented. `k_corr_step.py` preforms it's function.


Parameter Folders
---------------
[params_lc_matchup](params_lc_matchup)
* configuration/parameter files for the deprecated `lightcone_matchup.py` script.

[params_kcorr](params_kcorr)
* configuration/parameter files for `k_corr_step.py` script.

[params_lc_resamp](params_lc_resamp)
* configuration/parameter files for `lc_resample.py`











