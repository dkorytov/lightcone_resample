Light Cone Resampling of Galacticus Galaxies
============================================

This package combines baseDC2 and Galacticus catalogs to produce
cosmoDC2. 

There is a description of each file in
[FILE_DESCRIPTION.md](FILE_DESCRIPTION.md). Most files are diagnostics
or plotting routines.

Blackbox Overview
=================

![Blackbox Overview](doc_figures/Matchup%20Pipeline.png)

The main matchup pipeline is `lc_resample.py` which produces
the cosmoDC2 catalogs. The matchup pipeline requires two preproccess data sets to run. `precompute_masks.py` precomputes which galaxies pass galmatcher's color cuts. `k_corr_step.py` 

Each script takes in a single configuration/parameter file that
specifies input/output file locations and various run
settings. `lc_reample.py` and `k_corr_step.py` runs off .param files
located in `param_lc_resamp/..` and `param_k_corr/...`,
respectively. `precompute_mask.py` runs off the same .param file as
`lc_resample.py`. A full description of the parameter files and the settings in them can be
found in [PARAM_DESCRIPTION.md](PARAM_DESCRIPTION.md)

lc_resample.py
--------------

precompute_masks.py
-------------------

k_corr_step.py
--------------

Internal Overview
=================

This section reviews the many internal processes in
`lc_resample.py`. Below is a workflow chart. The purple cylinders are
files on disk, the orange boxes are processes, the red oval is an
internal arrays in memory and the teal wavy-shape is the configuration
file that controls. 

![Whitebox Overview](doc_figures/Internal%20Matchup%20Pipeline.png)

In short, the scripts loads all run time specifications and
input/output files locations from a parameter file. The script runs
mainly on time step intervals, producing almost complete cosmoDC2
catalogs for each time step (which I'm calling "step" files). To
produce these intermediate files, it loads in the time step slice of
the baseDC2 lightcone galaxy and the two Galacticus catalogs that
straddle the redshift range of the lightcone section and finds matches
between between Galacticus and baseDC2 galaxies. The matching is done
by splitting up the time step interval into even smaller redshift
ranges (called "substeps"), interpolating Galacitcus galaxies to the
substep redshift and preforming a nearest neighbor search by restframe
Mr, g-r, r-i colors for field galaxies and additional observer frame
g-r, r-i and i-z colors for cluster red sequence galaxies. Once all
the substeps are done and we have matches for all the galaxies within
the time step intervals, properties from matched Galacticus galaxies
are interpolated to the exact redshift of the lightcone baseDC2 galaxy
and additional properties are added on top (such as size, shear,
stellar mass, etc) of what Galacticus provides. All these properties
are written to disk as a "step" file. Once all the step files are
produced, a new file for all the steps is created (called the "all"
file because the step number is replace by "all") which is effectively
a simple concatenation. After the concatenation, metadata and column
units are added to produce the final cosmoDC2 catalog.


Workflow Walk Through 
---------------------

The `lc_resample.py` is call with a single argument which is the path
to the configuration or parameter file. These parameter files control
all run time options and point to all input and output files. The
parameter files for `lc_resample.py` are stored in the
`params_lc_resamp/..` directory. A full description of the param file
variables can be found in [PARAM_DESCRIPTIONS.md](PARAM_DESCRIPTIONS.md).

The main function in the script is `lightcone_resample(param_file)`.
