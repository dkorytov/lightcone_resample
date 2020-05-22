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

There are bunch of internal processes 
![Whitebox Overview](doc_figures/Internal%20Matchup%20Pipeline.png)

