Parameter/Configure File Description
====================================

Describes all the variables in parameter files. The order is the same
as found in saved param files. These parameter files follow a similar
HACC formating. Some strings have either "${step}" or "${healpix}"
that get replace for the specified files by step/healpix pixel.

Param file for `lc_resample.py`, in `param_lc_resamp/`
---------------

* `lightcone_fname` 
  - The location of baseDC2 light cone. 
* `gltcs_fname`
  - Galaticus library file location. 
* `gltcs_metadata_ref`
  - A single Galacitcus library file from which the metadata is
    pulled. The step doesn't matter since the metadata is the same for
    all files.
* `gltcs_slope_fname`
  - **deprecated** Used for older method of linear redshift interpolation of galaxy library properties 
* `halo_shape_fname`
  - **deprecated** Used for halo shape information. It's wrapped into baseDC2. 
* `halo_shape_red_fname`
  - **deprecated** Used for halo shape information. It's wrapped into baseDC2. 
* `sod_fname`
  - **deprecated** Used for sod halo information. It's wrapped into baseDC2. 
* `healpix_shear_file`
  - file location of calculated lightcone galaxy shears.
* `healpix_file`
  - boolean if we are running on multiple healpixels or not.
* `healpix_pixels`
  - the list of healpix pixels numbers to run together
* `fake_lensing`
  - Ignore the shear file and insert zero shear into the catalog
* `steps`
  - A list of time steps to run on. Must be in descending order. 


* `output_fname`
  - The output of for the catalog. As an intermediate step, one file
    will be generated for every time step interval with "${step}"
    being replaced with a time step. The final catalog will
    concatenate all the "step" files into a final "all" file. The
    intermediate step files can be deleted.

* `load_mask`
  - Boolean if to load precomputed masks (fast) or computing on the go (slow)
* `mask_loc`
  - The file location for precomputed galmatcher color cuts
* `index_loc`
  - The location of index file to match the same galaxies across time
    steps in the galaxy library
* `use_slope`
  - Must be set to `true`. Sets the pipeline to linearly interpolate
    galaxy properties--the only option as other options are deprecated.

* `substeps`
  - How many substeps to use between time steps for the match up. The
    lightcone galaxies are split into redshift groups/substeps and the
    library galaxies are interpolated to the middle red shift. The
    matching is only done within substeps. More substeps give a better
    matchup, but take longer. 5 substeps perform quite well and
    doesn't take too long.
	

* `use_substep_redshift`
  - boolean if to use the substep redshift or the galaxy redshift for
    galaxy property interpolation. This should always be set to
    `false` in production to avoid color-z banding. For quick test
    runs with few substeps, setting to `true` will give a better
    matchup but introduce color-z banding.

* `recolor`
  - should be set `false`. It recolors lightcone baseDC2 galaxies by
    calling the galaxy color functions in the cosmoDC2 package. Used
    to test new color distribution in baseDC2 withtout generating a
    new catalog.
	
* `short`
  - reduces run time by excluding SED columns from the output
    catalog. Should be set to `false` for production, of course.
  
* `supershort`
  - reduces run time even more by excluding almost everything but SDSS
    and LSST colors. Should be set to `false` for production, of
    course.
	
* `cut_small_galaxies`
  - should be set to `false`. Removes small galaxies from the catalog
    to reduce run time. 
	
* `cut_small_galaxies_mass`
  - mass limit for the above 

* `use_dust_factor`
  - must be set to `true` (all other options deprecated) even if no
    additional dust factors are used. Turns on the ability to multiply
    the effect of dust on galaxy properties (i.e. twice the reddening, etc)
  
* `dust_factors`
  - should be empty. The list of dust factors to apply, there's an
    implicit factor of 1 that is automatically added. Adding dust
    factors other than the implicit 1, will introduce color-z banding
    in red sequence cluster galaxies. 

* `ignore_mstar`
  - must be set to `true`. For the matchup between baseDC2 and the
    library, don't try to match the stellar mass, only the colors and
    magnitudes. Galacitcus's stellar mass is overwritten by baseDC2
    stellar mass.

* `match_obs_color_red_seq`
  - must be set to `true`. Forces red sequence galaxies to match
    expected observed colors in auditing to rest frame colors. 
  
* `red_sequence_transition_mass_start`
  - should be set to `12.5`. All red sequence galaxies in halos below this mass cut don't have
    to match observed colors. In between lower and upper limits, a
    linear proportion from 0% to 100% of red sequence galaxies are
    required to match observed colors
	
* `red_sequence_transition_mass_end`
  - should be set to `13.5` All red sequence galaxies above this cut have to match
    observed colors. In between lower and upper limits, a
    linear proportion from 0% to 100% of red sequence galaxies are
    required to match observed colors

* `red_sequence_scatter_query_gr`
  - should be set to `0.09`. Scatter added to baseDC2 expected red
    sequence colors. Tuned to get the right amount of scatter in
    cluster red sequence.
  
* `red_sequence_scatter_query_ri`
  - should be set to `0.114`. Scatter added to baseDC2 expected red
    sequence colors. Tuned to get the right amount of scatter in
    cluster red sequence.
  
* `red_sequence_scatter_query_iz`
  - should be set to `0.10`. Scatter added to baseDC2 expected red
    sequence colors. Tuned to get the right amount of scatter in
    cluster red sequence.
  
* `red_sequence_scatter_tree_gr`
  - should be set to `0.09`. Scatter added to Galacticus red sequence
    colors before matchup. Tuned to get the right amount of scatter in
    cluster red sequence.
* `red_sequence_scatter_tree_ri`
  - should be set to `0.114`. Scatter added to Galacticus red sequence
    colors before matchup. Tuned to get the right amount of scatter in
    cluster red sequence.
* `red_sequence_scatter_tree_iz`
  - should be set to `0.10`. Scatter added to Galacticus red sequence
    colors before matchup. Tuned to get the right amount of scatter in
    cluster red sequence.
  
* `rescale_bright_luminosity`
  - should be set to `true`.  boolean if to rescale luminosity of
    matched galacticus galaxy to the same level of the matched baseDC2
    galaxy.
* `rescale_bright_luminosity_threshold`
  - should be set to `0`. The magnitude threshold to start rescaling
    the magnitude. With a `0` value all galaxies should be rescaled.

* `ignore_bright_luminosity`
  - should be set to `true`. For bright galaxies, don't try to match
    the luminosity between baseDC2 and Galacticus. This is done by
    compressing the bright range into a narrower range. 
* `ignore_bright_luminosity_threshold`
  - should be set to `-21`. The maximum brightness value in the compressed range. 
* `ignore_bright_luminosity_softness`
  - should be set to `2`. The softness of the transition of not compressed to compressed
    Values below threshold minus softness are not compressed.
* `plot` 
  - boolean if to display various diagnostic plots to display during
    runtime
* `plot_substep`
  - boolean if to display various diagnostic plots related to
    substepping to display during runtime

* `version_major`
  - major version number of catalog
* `version_minor`
  - minor version number of catalog
* `version_minor_minor`
  - revision/minor-minor version number of catalog


* `concatenate_only`
  - should be set to `false` for production runs. If set to true, it
    will skip the matchup process and work from the intermediate
    "steps" files to produce the "all" file. (see description for
    `output_fname`). Allows for fixing errors or bugs that occurred
    between the intimated "step" file generation and final "all" file
    generation.
