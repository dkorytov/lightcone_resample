#!/usr/bin/env python2.7
from __future__ import print_function, division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk
import h5py
import time
import sys
import datetime
from astropy.table import Table
from scipy.spatial import cKDTree
from pecZ import pecZ
from astropy.cosmology import WMAP7 as cosmo
from scipy.interpolate import interp1d 


def load_mags(fname,frame,use_dust=False):
    dic ={}
    hgroup = h5py.File(fname,'r')
    if use_dust:
        dust = ':dustAtlas'
    else:
        dust = ''
    dic['mag_g'] = hgroup['galaxyProperties/SDSS_filters/magnitude:SDSS_g:'+frame+dust].value
    dic['mag_r'] = hgroup['galaxyProperties/SDSS_filters/magnitude:SDSS_r:'+frame+dust].value
    dic['mag_i'] = hgroup['galaxyProperties/SDSS_filters/magnitude:SDSS_i:'+frame+dust].value
    dic['clr_gr'] = dic['mag_g']-dic['mag_r']
    dic['clr_ri'] = dic['mag_r']-dic['mag_i']
    return dic

def multiply(dic,val):
    new_dic = {}
    for key in dic.keys():
        new_dic[key] = dic[key]*val
    return new_dic

def add(dic1, dic2):
    new_dic = {}
    for key in dic1.keys():
        new_dic[key] = dic1[key]+dic2[key]
    return new_dic;

def plot_colors(mags,title,frame):
    plt.figure()
    ybins = np.linspace(-1,2,100)
    if frame == 'rest':
        xbins = np.linspace(-26,-12,100)
    elif frame == 'observed':
        xbins = np.linspace(12,26,100) 
    elif frame == 'slope':
        xbins = np.linspace(-5,5,100)
        ybins = np.linspace(-5,5,100)
    else:
        print("frame \"{}\" is not allowed".format(frame))
        raise KeyError
    h,xbins,ybins = np.histogram2d(mags['mag_r'],mags['clr_gr'],bins=(xbins,ybins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.grid()
    plt.xlabel('mag_r');plt.ylabel('g-r')
    plt.title(title)

if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    gltcs_fname = param.get_string('gltcs_fname')
    gltcs_slope_fname = param.get_string('gltcs_slope_fname')
    steps = param.get_int_list('steps')
    stepz = dtk.StepZ(sim_name='AlphaQ')
    frame ='rest'
    frame ='observed'
    for i in range(0,len(steps)-1):
        step1 = steps[i+1]
        step2 = steps[i]
        a_1 = stepz.get_a(step1)
        a_2 = stepz.get_a(step2)
        print(a_1,a_2)
        del_a = a_2 - a_1
        mag1 = load_mags(gltcs_fname.replace('${step}',str(step1)),frame,use_dust=True)
        mag2 = load_mags(gltcs_fname.replace('${step}',str(step2)),frame,use_dust=True)
        slope1 = load_mags(gltcs_slope_fname.replace('${step}',str(step1)),frame, use_dust=True)
        plot_colors(mag1,str(step1),frame)

        for factor in np.linspace(0,1,25):
            mag_slope = add(mag1,multiply(slope1,del_a*factor))
            plot_colors(mag_slope,"{:.2f} of the way from {} to {}.".format(factor,step1,step2),frame)
        plot_colors(mag2,str(step2),frame)
        for factor in np.linspace(1,0,25):
            mag_slope = add(mag1,multiply(slope1,del_a*factor))
            plot_colors(mag_slope,"{:.2f} of the way from {} to {}.".format(factor,step1,step2),frame)


        plot_colors(multiply(slope1,del_a),str(step1),'slope')
        dtk.save_figs('figs/'+sys.argv[1]+'/'+__file__+"/")

