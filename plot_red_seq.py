#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk 
import h5py
import sys
import time
from numpy.random import normal


def plot_colors(fname,z_min,z_max,title,mass_cut):
    hgroup = h5py.File(fname,'r')['galaxyProperties']
    redshift = hgroup['redshift'].value
    mag_g = hgroup['SDSS_filters/magnitude:SDSS_g:observed:dustAtlas'].value
    mag_r = hgroup['SDSS_filters/magnitude:SDSS_r:observed:dustAtlas'].value
    #mag_i = hgroup['SDSS_filters/magnitude:SDSS_i:rest:dustAtlas'].value
    host_mass = hgroup['hostHaloMass'].value
    slct = (z_min < redshift) & (redshift < z_max) & (mass_cut < host_mass)
    r = mag_r[slct]
    gr = mag_g[slct] - mag_r[slct]
    ybins = np.linspace(-0.5,2,250)
    xbins = np.linspace(14,30,250)
    h,xbins,ybins = np.histogram2d(r,gr,bins=(xbins,ybins))
    plt.figure()
    plt.title(title+"\n{:.2f} < z < {:.2f})".format(z_min,z_max))
    plt.plot(r,gr,'.',alpha=0.3)
    # plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    # plt.colorbar()
    plt.grid()
    plt.xlabel('mag r');plt.ylabel('g-r color')
    plt.tight_layout()
    
def plot_color_z(fname,title,mass_cut):
    hgroup = h5py.File(fname,'r')['galaxyProperties']
    redshift = hgroup['redshift'].value
    host_mass = hgroup['hostHaloMass'].value
    mag_g = hgroup['SDSS_filters/magnitude:SDSS_g:observed:dustAtlas'].value
    mag_r = hgroup['SDSS_filters/magnitude:SDSS_r:observed:dustAtlas'].value
    gr = mag_g - mag_r
    ybins = np.linspace(-0.5,2,250)
    xbins = np.linspace(0,1,250)
    slct = (mass_cut < host_mass)
    print("color_z num: ", np.sum(slct))

    h,xbins,ybins = np.histogram2d(redshift[slct],gr[slct],bins=(xbins,ybins))
    plt.figure()
    #plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    #plt.colorbar()
    plt.plot(redshift[slct],gr[slct],'.',alpha=0.3)
    
    plt.grid()
    plt.xlabel('redshift');plt.ylabel('g-r')
    plt.tight_layout()

if __name__ == "__main__":
    title = sys.argv[1]
    fname = sys.argv[2]
    mass_cut = 1e14
    plot_color_z(fname,title,mass_cut)
    zbins = np.linspace(0.28,0.3,4)
    for i in range(0,len(zbins)-1):
        plot_colors(fname,zbins[i],zbins[i+1],title,mass_cut)
    plt.show()
    
