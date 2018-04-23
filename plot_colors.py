#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk 
import h5py
import sys
import time
from numpy.random import normal

def plot_clr_clr(v3_data, key1, key2, title):
    plt.figure()
    bins = np.linspace(-1,2,100)
    h,xbins,ybins = np.histogram2d(gr,ri,bins=(bins,bins))
    plt.pcolor(xbins, ybins, h.T, cmap='PuBu', norm=clr.LogNorm())
    plt.grid()
    plt.xlabel('g-r')
    plt.ylabel('r-i')
    plt.title(title)

def load_protoDC2(fname):
    t1 = time.time()
    print "loading protodc2...",
    hgroup = h5py.File(fname,'r')['galaxyProperties']
    redshift = hgroup['redshiftHubble'].value
    result  = {}
    mag_g = hgroup['SDSS_filters/magnitude:SDSS_g:rest:dustAtlas'].value
    mag_r = hgroup['SDSS_filters/magnitude:SDSS_r:rest:dustAtlas'].value
    mag_i = hgroup['SDSS_filters/magnitude:SDSS_i:rest:dustAtlas'].value
    mag_z = hgroup['SDSS_filters/magnitude:SDSS_z:rest:dustAtlas'].value
    result['g-r rest'] = mag_g - mag_r
    result['r-i rest'] = mag_r - mag_i
    result['i-z rest'] = mag_i - mag_z
    result['mag g rest'] = mag_g
    result['mag r rest'] = mag_r 
    result['mag i rest'] = mag_i
    result['mag z rest'] = mag_z

    mag_g = hgroup['SDSS_filters/magnitude:SDSS_g:observed:dustAtlas'].value
    mag_r = hgroup['SDSS_filters/magnitude:SDSS_r:observed:dustAtlas'].value
    mag_i = hgroup['SDSS_filters/magnitude:SDSS_i:observed:dustAtlas'].value
    mag_z = hgroup['SDSS_filters/magnitude:SDSS_z:observed:dustAtlas'].value
    result['mag g obs'] = mag_g
    result['mag r obs'] = mag_r
    result['mag i obs'] = mag_i
    result['mag z obs'] = mag_z
    result['g-r obs'] = mag_g - mag_r
    result['r-i obs'] = mag_r - mag_i
    result['i-z obs'] = mag_i - mag_z
    result['redshift'] = redshift
    result['dust factor'] = hgroup['dustFactor'].value
    print(result['dust factor'])
    print "\n\ttime: ",time.time() - t1
    return result

def load_umachine(fname):
    hfile = h5py.File(fname,'r')
    result = {}
    result['mag r rest']= hfile['restframe_extincted_sdss_abs_magr'].value
    result['g-r'] = hfile['restframe_extincted_sdss_gr'].value
    result['r-i'] = hfile['restframe_extincted_sdss_ri'].value
    result['redshift'] = hfile['redshift'].value
    return result

def append_dics(dics):
    result = {}
    keys = dics[0].keys()
    for key in keys:
        result[key] = []
    for dic in dics:
        for key in keys:
            result[key].append(dic[key])
    for key in keys:
        result[key] = np.concatenate(result[key])
    return result



if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    lightcone_fname = param.get_string("lightcone_fname")
    output_fname = param.get_string("output_fname")
    steps = param.get_int_list("steps")
    v1 = param.get_int("version_major")
    v2 = param.get_int("version_minor")
    v3 = param.get_int("version_minor_minor")
    lc_dics = []
    for i,step in enumerate(steps):
        if i == 0:
            continue
        lc_dics.append(load_umachine(lightcone_fname.replace("${step}",str(step))))
    lc_data = append_dics(lc_dics)
    
    plt.figure()
    h,xbins,ybins = np.histogram2d(lc_data['redshift'],lc_data['g-r'],bins=(250,250))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('g-r rest'); plt.xlabel('redshift')
    plt.title("UMachine+SDSS Light Cone")
    plt.grid()

    plt.figure()
    h,xbins,ybins = np.histogram2d(lc_data['redshift'],lc_data['r-i'],bins=(250,250))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('r-i rest'); plt.xlabel('redshift')
    plt.title("UMachine+SDSS Light Cone")
    plt.grid()
    clrclrbins = np.linspace(-.5,1,250)

    plt.figure()
    h,xbins,ybins = np.histogram2d(lc_data['g-r'],lc_data['r-i'],bins=(clrclrbins,clrclrbins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.xlabel('g-r rest'); plt.ylabel('r-i rest')
    plt.title("UMachine+SDSS Light Cone")
    plt.grid()

    
    protoDC2 = load_protoDC2(output_fname.replace("${step}","all"))
    bins = np.linspace(-1,3,250)
    zbins = np.linspace(0,1,250)
    magobsbins = np.linspace(14,39,100)
    magrestbins = np.linspace(-30,-12,100)
    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['redshift'],protoDC2['mag r rest'],bins=(zbins,bins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('Mag r rest'); plt.xlabel('redshift')
    plt.title("ProtoDC2 v{}.{}.{}".format(v1,v2,v3))
    plt.grid()

    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['redshift'],protoDC2['mag g rest'],bins=(zbins,bins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('Mag g rest'); plt.xlabel('redshift')
    plt.grid()
    plt.title("ProtoDC2 v{}.{}.{}".format(v1,v2,v3))

    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['redshift'],protoDC2['g-r obs'],bins=(zbins,bins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('g-r observed'); plt.xlabel('redshift')
    plt.title("ProtoDC2 v{}.{}.{}".format(v1,v2,v3))
    plt.grid()

    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['redshift'],protoDC2['r-i obs'],bins=(zbins,bins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('r-i observed'); plt.xlabel('redshift')
    plt.title("ProtoDC2 v{}.{}.{}".format(v1,v2,v3))
    plt.grid()

    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['redshift'],protoDC2['i-z obs'],bins=(zbins,bins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('i-z observed'); plt.xlabel('redshift')
    plt.title("ProtoDC2 v{}.{}.{}".format(v1,v2,v3))
    plt.grid()


    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['g-r rest'],protoDC2['r-i rest'],bins=(clrclrbins,clrclrbins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.xlabel('g-r rest'); plt.ylabel('r-i rest')
    plt.title("ProtoDC2 v{}.{}.{}".format(v1,v2,v3))
    plt.grid()


    clrclrbins= np.linspace(-1,2,250)

    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['g-r obs'],protoDC2['r-i obs'],bins=(clrclrbins,clrclrbins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.xlabel('g-r observed'); plt.ylabel('r-i observed')
    plt.title("ProtoDC2 v{}.{}.{}".format(v1,v2,v3))
    plt.grid()


    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['mag r rest'], protoDC2['g-r rest'],bins=(magrestbins,clrclrbins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.xlabel('r rest'); plt.ylabel('g-r rest')
    plt.title("ProtoDC2 v{}.{}.{}".format(v1,v2,v3))
    plt.grid()

    slct = protoDC2['mag r rest']<-18
    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['g-r obs'][slct],protoDC2['r-i obs'][slct],bins=(clrclrbins,clrclrbins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.xlabel('g-r observed'); plt.ylabel('r-i observed')
    plt.title("ProtoDC2 v{}.{}.{}\nMr<-18".format(v1,v2,v3))
    plt.grid()

    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['redshift'][slct],protoDC2['g-r obs'][slct],bins=(zbins,bins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('g-r observed'); plt.xlabel('redshift')
    plt.title("ProtoDC2 v{}.{}.{}\nMr<-18".format(v1,v2,v3))
    plt.grid()

    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['redshift'][slct],protoDC2['r-i obs'][slct],bins=(zbins,bins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('r-i observed'); plt.xlabel('redshift')
    plt.title("ProtoDC2 v{}.{}.{}\nMr<-18".format(v1,v2,v3))
    plt.grid()

    plt.figure()
    h,xbins,ybins = np.histogram2d(protoDC2['redshift'][slct],protoDC2['i-z obs'][slct],bins=(zbins,bins))
    plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    plt.ylabel('i-z observed'); plt.xlabel('redshift')
    plt.title("ProtoDC2 v{}.{}.{}\nMr<-18".format(v1,v2,v3))
    plt.grid()


    # plt.figure()
    # h,xbins,ybins = np.histogram2d(protoDC2['redshift'],protoDC2['mag i obs'],bins=(zbins,magobsbins))
    # plt.pcolor(xbins,ybins,h.T,cmap="PuBu",norm=clr.LogNorm())
    # plt.colorbar()
    # plt.plot([0,1],[25.3,25.3],'r',label='gold sample')
    # plt.plot([0,1],[26.8, 26.8], 'r--',label='i<26.8')
    # plt.legend(loc='best',framealpha=0.5)
    # plt.xlabel('redshift');plt.ylabel('mag i observed')
    # plt.grid()
    # plt.xlim([0,1])
    # plt.tight_layout()

    # plt.figure()
    # h,xbins,ybins = np.histogram2d(protoDC2['redshift'],protoDC2['mag r obs'],bins=(zbins,magobsbins))
    # plt.pcolor(xbins,ybins,h.T,cmap="PuBu",norm=clr.LogNorm())
    # plt.colorbar()
    # plt.legend(loc='best',framealpha=0.5)
    # plt.xlabel('redshift');plt.ylabel('mag r observed')
    # plt.grid()
    # plt.xlim([0,1])
    # plt.tight_layout()

    # plt.figure()
    # h,xbins,ybins = np.histogram2d(protoDC2['redshift'],protoDC2['mag g obs'],bins=(zbins,magobsbins))
    # plt.pcolor(xbins,ybins,h.T,cmap="PuBu",norm=clr.LogNorm())
    # plt.colorbar()
    # plt.legend(loc='best',framealpha=0.5)
    # plt.xlabel('redshift');plt.ylabel('mag g observed')
    # plt.grid()
    # plt.xlim([0,1])
    # plt.tight_layout()

    plt.figure()
    h, xbins = np.histogram(protoDC2['redshift'],bins=250)
    plt.plot(dtk.bins_avg(xbins), h, label='total')

    slct = protoDC2['dust factor'] == 1.0
    h, xbins = np.histogram(protoDC2['redshift'][slct], bins=250)
    plt.plot(dtk.bins_avg(xbins), h, label='dust factor = 1')

    slct = protoDC2['dust factor'] == 3.0
    h, xbins = np.histogram(protoDC2['redshift'][slct], bins=250)
    plt.plot(dtk.bins_avg(xbins), h, label='dust factor = 3')

    slct = protoDC2['dust factor'] == 6.0
    h, xbins = np.histogram(protoDC2['redshift'][slct], bins=250)
    plt.plot(dtk.bins_avg(xbins), h, label='dust factor = 6')

    plt.xlabel('redshift')
    plt.ylabel('count')
    plt.legend(loc='best')
    plt.yscale('log')
    plt.grid()

    dtk.save_figs("figs/"+sys.argv[1]+"/"+__file__+"/")
    plt.show()
