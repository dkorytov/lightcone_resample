#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
import h5py 
import dtk
import sys
import time
import datetime


def load_gio(lc_gio_fname):
    lc = {}
    lc['x'] = dtk.gio_read(lc_gio_fname,'x')
    lc['y'] = dtk.gio_read(lc_gio_fname,'y')
    lc['z'] = dtk.gio_read(lc_gio_fname,'z')
    lc['vx'] = dtk.gio_read(lc_gio_fname,'vx')
    lc['vy'] = dtk.gio_read(lc_gio_fname,'vy')
    lc['vz'] = dtk.gio_read(lc_gio_fname,'vz')
    lc['id'] = dtk.gio_read(lc_gio_fname,'id')
    return lc

def load_hdf5(ss_gal_fname):
    ss = {}
    hfile = h5py.File(ss_gal_fname,'r')
    keys = hfile.keys()
    for key in keys:
        print(key)
        ss[key]=hfile[key].value
    return ss



def match_up(lc, ss, output):
    t1 = time.time()
    ss_id = ss['lightcone_id']
    srt = np.argsort(ss_id)
    indx = dtk.search_sorted(ss_id, lc['id'], sorter=srt)
    num_not_found = np.sum(indx == -1)
    print("num not found: ",num_not_found)
    hfile = h5py.File(output,'w')
    keys = ss.keys()
    for key in keys:
        print("\t",key,end="")
        t2 = time.time()
        hfile[key]=ss[key][indx]
        print("--",time.time()-t2)
    print('done. ',time.time()-t1)
    hfile.close()


if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    lightcone_gio_fname = param.get_string("lightcone_gio_fname")
    snapshot_galaxy_fname = param.get_string("snapshot_galaxy_fname")
    lightcone_hdf5_fname = param.get_string("lightcone_hdf5_fname")
    steps = param.get_int_list("steps")
    for step in steps:
        match_up(lightcon_gio_fname.replace("${step}",str(step)),
                 snapshot_galaxy_fname.replace("${step}",str(step)),
                 lightcone_hdf5_fname.replace("${step}",str(step)))
