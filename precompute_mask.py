#!/usr/bin/env python2.7
from __future__ import print_function, division

import numpy as np
import scipy as sp
import pdb
import dtk
import h5py
import time
import sys
import datetime

import galmatcher




if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    steps = param.get_int_list("steps")
    gltcs_fname = param.get_string("gltcs_fname")
    mask_loc = param.get_string("mask_loc")
    selection1 = galmatcher.read_selections(yamlfile='galmatcher/yaml/vet_protoDC2.yaml')
    selection2 = galmatcher.read_selections(yamlfile='galmatcher/yaml/colors_protoDC2.yaml')
    hfile = h5py.File(mask_loc, 'w')
    for step in steps:
        print("\n\n=============")
        print("STEP {}".format(str(step)))
        gltcs_step_fname = gltcs_fname.replace('${step}',str(step))
        mask1 = galmatcher.mask_cat(h5py.File(gltcs_step_fname, 'r'), selections=selection1)
        mask2 = galmatcher.mask_cat(h5py.File(gltcs_step_fname, 'r'), selections=selection2)
        mask_all = mask1 & mask2
        hfile['{}'.format(str(step))] = mask_all
        
