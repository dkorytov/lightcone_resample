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
    #plt.plot(r,gr,'.',alpha=0.3)
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.colorbar()
    plt.grid()
    plt.xlabel('mag r');plt.ylabel('g-r color')
    plt.tight_layout()
    
mag_trans = {'g':'LSST_filters/magnitude:LSST_g:observed:dustAtlas',
             'r':'LSST_filters/magnitude:LSST_r:observed:dustAtlas',
             'i':'LSST_filters/magnitude:LSST_i:observed:dustAtlas',
             'z':'LSST_filters/magnitude:LSST_z:observed:dustAtlas',
             'y':'LSST_filters/magnitude:LSST_y:observed:dustAtlas',}

def plot_color_z(fname,title,mass_cut, mag1, mag2,central_cut=False):
    hgroup = h5py.File(fname,'r')['galaxyProperties']
    print(hgroup['LSST_filters'].keys())
    redshift = hgroup['redshift'].value
    host_mass = hgroup['hostHaloMass'].value
    mag1_val = hgroup[mag_trans[mag1]].value
    mag2_val = hgroup[mag_trans[mag2]].value
    clr_mag = mag1_val - mag2_val
    ybins = np.linspace(-0.5,2,250)
    xbins = np.linspace(0,1,250)
    slct = (mass_cut < host_mass)
    if central_cut:
        central = hgroup['isCentral'].value
        slct = slct & (central == 1)
    print("color_z num: ", np.sum(slct))

    h,xbins,ybins = np.histogram2d(redshift[slct],clr_mag[slct],bins=(xbins,ybins))
    plt.figure(figsize=(7,5))
    # plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    # plt.colorbar()
    plt.plot(redshift[slct],clr_mag[slct],'+',alpha=0.3)
    
    plt.grid()
    plt.xlabel('redshift');plt.ylabel('LSST {}-{}'.format(mag1,mag2))
    plt.title(title)
    plt.tight_layout()


if __name__ == "__main__":
    title = sys.argv[1]
    fname = sys.argv[2]
    mass_cut = 3e13
    # plot_color_z(fname,title,mass_cut,'g','r', True)
    plot_color_z(fname,title,mass_cut,'r','i', True)
    # plot_color_z(fname,title,mass_cut,'i','z', True)
    # plot_color_z(fname,title,mass_cut,'z','y', True)
    # for i in range(0,len(zbins)-1):
    #     plot_colors(fname,zbins[i],zbins[i+1],title,mass_cut)
    dtk.save_figs(path='figs/'+__file__+"/"+title+"/"+fname+"/")
    plt.show()
    
