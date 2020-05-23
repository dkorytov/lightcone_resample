Light Cone Resampling of Galacticus Galaxies
============================================

This package combines baseDC2 lightcone and Galacticus galaxy catalogs
to produce the cosmoDC2 galaxy catalog.

There is a description of each file in
[FILE_DESCRIPTIONS.md](FILE_DESCRIPTIONS.md). Most files are
diagnostics or plotting routines.

Broad Overview
=================

![Blackbox Overview](doc_figures/Matchup%20Pipeline.png)

The main matchup pipeline is `lc_resample.py` which produces the
cosmoDC2 catalogs. The matchup pipeline requires two preproccess data
sets to run: precomputed masks produced by `precompute_masks.py` and
indexes to identify the same galaxies in Galacticus across time steps
produced by `k_corr_step.py`.

Each script takes in one argument which is a single
configuration/parameter file that specifies input/output file
locations and various run settings. `lc_reample.py` and
`k_corr_step.py` runs off .param files located in `param_lc_resamp/..`
and `param_k_corr/...`, respectively. `precompute_mask.py` runs off
the same .param file as `lc_resample.py`. A full description of the
parameter files and the settings in them can be found in
[PARAM_DESCRIPTIONS.md](PARAM_DESCRIPTIONS.md)

Precomputed Galaxy Masking, `precompute_masks.py`
-------------------

Not all Galacticus galaxy have reasonable colors or other
properties. Galaxies with unreasonable properties are excluded from
being used in the match up process. Eve's Galmatcher package
determines which galaxies are acceptable. Unacceptable galaxies are
excluded by applying a mask to the Galacticus galaxy array. Unless the
Galacticus catalog or Galmatcher's requirements change, the masks only
needs to be generated once.

`precompute_masks.py` runs on the same .param file as for
`lc_resample.py`. The param file specifies Galacticus files locations,
time steps and output mask location (the mask location is an input for
`lc_resample.py`)

Precomputed Galacticus Interpoliation Indexes, `k_corr_step.py`
--------------

A core feature of the match up pipeline is the redshift interpolation
of Galacticus galaxy properties. The Galacticus catalog only has
descrete redshifts while we have continous redshifts in the baseDC2
lightcone. If galaxy properties are not interpolated, there will be
descrete bands in color-z plots, see figures below.

![](doc_figures/redshift_no_interpolation.png)![](doc_figures/redshift_interpolation.png)

The redshift interpolation requires to find the same galaxy in two
adjacent time step Galacticus catalogs. The order of galaxies are not
the same, so an indexing array is created to reorder step_i+1 catalog
order to match step_i, see table below for an example.
| step_i | step_i+1 | indexing | reorder step_i+1|
|--------------------| -- | -- | -- |
| galaxy_0 | galaxy_3 | 5 | galaxy_0 |
| galaxy_1 | galaxy_1 | 1 | galaxy_1 |
| galaxy_2 | galaxy_4 | -1| - |
| galaxy_3 | galaxy_5 | 0 | galaxy_3 |
| galaxy_4 | galaxy_6 | 2 | galaxy_4 |
|      | galaxy_0 |  | |

There are a couple complications:
* Central galaxies in Galacticus don't have unique ids--only satellite
  galaxies keep the same id between timesteps. To track centrals we
  have to load in the halo merger trees to figure out where central
  galaxies go.
* Not each galaxy exists the following time step, such as
  galaxy_2. This galaxy is simply excluded from the match up pipeline.
* A galaxy may fail Galmatcher's color cuts for either step_i or
  step_i+1 but not both. To make sure a galaxy isn't interpolated into
  cut region, that galaxy is excluded from the matchup.

Matchup Pipeline, `lc_resample.py`
--------------

 This is main body of the pipeline that does the matchup and generates
 the cosmoDC2 catalog. Like the other scripts, it takes in a single
 argument which is a paramter file. The .param file describes run time
 options and all input/output file locations. A full description of
 the param file can be found in
 [PARAM_DESCRIPTIONS.md](PARAM_DESCRIPTIONS.md). The pipeline takes in
 lighcone baseDC2, the Galaticus library, lightcone shears,
 precomputed masks and interpolation indexes. The lightcone shears
 need to have exactly the same file structure and galaxy order as
 baseDC2. The pipeline produces intermediate files for each time step
 interval and the final catalog. The intermediate files can be deleted
 afterwards or kept to skip a part of the pipeline during a rerun. A
 full description of what goes on inside the matchup pipe is below.


Internal Overview of the Matchup Pipeline, `lc_resample.py`
=================

This section reviews the many internal processes in
`lc_resample.py`. Below is a workflow chart. The purple cylinders are
files on disk, the orange boxes are processes, the red oval is an
internal arrays in memory and the teal wavy-shape is the configuration
file that controls. 

![Whitebox Overview](doc_figures/Internal%20Matchup%20Pipeline.png)


The scripts loads all run time specifications and input/output files
locations from a parameter file. The script runs mainly on time step
intervals, producing almost complete cosmoDC2 catalogs for each time
step, which I'm calling "step" files. The step files are combined to
produce the final catalog.

To produce these intermediate files, it loads in the time step slice
of the baseDC2 lightcone galaxy and the two Galacticus catalogs that
straddle the redshift range of the lightcone section and finds matches
between between Galacticus and baseDC2 galaxies for the time step. The
matching is done by further splitting up the time step interval into even
smaller redshift ranges (called "substeps"), interpolating Galacitcus
galaxies to the substep redshift and preforming a nearest neighbor
search by restframe Mr, g-r, r-i colors for field galaxies and
additional observer frame g-r, r-i and i-z colors for cluster red
sequence galaxies.

Once all the substeps are done and we have matches for all the
galaxies within the time step intervals, properties from matched
Galacticus galaxies are interpolated to the exact redshift of the
lightcone baseDC2 galaxy and additional properties are added on top
(such as size, shear, stellar mass, etc) of what Galacticus
provides. All these properties are written to disk as a "step"
file. 

After all the step files are produced, a new file for all the steps is
created (called the "all" file because the step number in the file
name is replace by "all") which is effectively a simple
concatenation. After the concatenation, metadata and column units are
added to produce the final cosmoDC2 catalog.

Workflow Walk Through 
---------------------

This is more detailed walk through the script. This section is written
as if I'm going over the code line-by-line-ish. 

The `lc_resample.py` is call with a single argument which is the path
to the configuration or parameter file. These parameter files control
all run time options and point to all the input and output files. The
parameter files for `lc_resample.py` are stored in the
`params_lc_resamp/..` directory. A full description of the param file
variables can be found in
[PARAM_DESCRIPTIONS.md](PARAM_DESCRIPTIONS.md).

The pipeline run starts with a function call to
`lightcone_resample(param_file)`, which I believe should be able to be
called by another script if `lc_resample.py` is imported. The function
loads the specified parameter file and pull outs all expected
parameters. Some parameters have default values that don't required to
be defined in the .param file. That's to ensure a current version of
the script will be able to run older parameter files before those
options were added. Immediately after extracting the variables, the
script processes and checks the variables.

Afterwards the script prepares how to handle Galmatchers masks. If
there are cached precomputed masks, it will load those masks. If there
aren't, it will create objects that calculate the mask on the
fly. It's faster to use the precomputed masks of course. The
precomputed masks are made by `precompute_masks.py`. 

Around this area, the script sets options depending if this is a
lightcone catalog (which is the typical run) or a snapshot
catalog. The main run time difference between lightcone and snapshot
catalogs is how the steps are handled. For lightcones, the pipeline
interpolates Galacticus galaxy properties between discrete time
steps. For snapshots, the pipeline doesn't need to interpolate since
all the Galacticus galaxies are already at the correct
redshift. Technically, the pipeline *still* interpolates properties
but between two of the same time step in Galacticus. This was an
simple solution to avoid writing large chunks of new code and
introducing bugs.

Now we start on making on intermediate "step" files. We will iterate
over the adjacent time steps pairs for lightcone or single time steps
for snapshots. First we check if we need to skip producing this "step"
file because of some option set in the paramter file. Then we get the
correct file names for this time step. Load in the cached galmatcher's
masks or calculate them. 

Hmm...this verbose parameter...I might just pull that into the param file
