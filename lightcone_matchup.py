#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
import h5py 
import dtk
import sys
import time
import datetime
from astropy.table import Table

def load_lightcone(lc_fname):
    print("loading lightcone...",end='')
    t1 = time.time()
    lc = {}
    lc['x'] = np.fromfile(lc_fname.replace("${var}","x"),dtype='f4')
    lc['y'] = np.fromfile(lc_fname.replace("${var}","y"),dtype='f4')
    lc['z'] = np.fromfile(lc_fname.replace("${var}","z"),dtype='f4')
    lc['vx'] = np.fromfile(lc_fname.replace("${var}","vx"),dtype='f4')
    lc['vy'] = np.fromfile(lc_fname.replace("${var}","vy"),dtype='f4')
    lc['vz'] = np.fromfile(lc_fname.replace("${var}","vz"),dtype='f4')
    lc['id'] = np.fromfile(lc_fname.replace("${var}","id"),dtype='i8')
    lc['ra'] = np.fromfile(lc_fname.replace("${var}","phi"),dtype='f4')
    lc['dec'] = np.fromfile(lc_fname.replace("${var}","theta"),dtype='f4')
    lc['redshift'] = np.fromfile(lc_fname.replace("${var}","redshift"),dtype='f4')
    print("done {}".format(time.time()-t1))
    return lc

def load_snapshot(ss_fname):
    print("loading snapshot...",end='')
    t1 = time.time()
    ss = {}
    tbl = Table.read(ss_fname,path='data')
    # hfile = h5py.File(ss_fname,'r')
    keys = tbl.keys()
    for key in keys:
         print("\t",key)
         ss[key]=tbl[key].quantity
    print("done. {}".format(time.time()-t1))
    return ss



def match_up(lc, ss, output):
    t1 = time.time()
    ss_id = ss['lightcone_id']
    srt = np.argsort(ss_id)
    indx = dtk.search_sorted(ss_id, lc['id'], sorter=srt)
    num_not_found = np.sum(indx == -1)
    print("num not found: ",num_not_found)
    assert(num_not_found == 0)
    hfile = h5py.File(output,'w')
    lc_keys = lc.keys()
    for key in lc_keys:
        if(key != 'id'):
            hfile[key] = lc[key]
    ss_keys = ss.keys()
    for key in ss_keys:
        print("\t",key,end="")
        if(key not in lc_keys):
            t2 = time.time()
            hfile[key]=ss[key][indx]
            print("--",time.time()-t2)
        else:
            print(" not copied")
    print('done. ',time.time()-t1)
    hfile.close()


if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    lightcone_bin_fname = param.get_string("lightcone_bin_fname")
    snapshot_galaxy_fname = param.get_string("snapshot_galaxy_fname")
    lightcone_output_fname = param.get_string("lightcone_output_fname")
    steps = param.get_int_list("steps")
    t0 =time.time()
    for step in steps:
        t1 = time.time()
        print("\n\n=====================\n STEP: {}".format(step))
        lc = load_lightcone(lightcone_bin_fname.replace("${step}",str(step)))
        ss = load_snapshot(snapshot_galaxy_fname.replace("${step}",str(step)))
        output_fname = lightcone_output_fname.replace("${step}",str(step))
        match_up(lc, ss, output_fname)
        print("\n=== done: {}".format(time.time()-t1))

    print("\n\n=======================\n=========================")
    print("All done: ",time.time()-t0)
