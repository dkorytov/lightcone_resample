#!/usr/bin/env python2.7
from __future__ import print_function, division
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr


def plot_differences(lc_data, gal_prop, index):
    keys = ['Mag_r','clr_gr','clr_ri','m_star']
    dist = {}
    dist_all = None
    for key in keys:
        d = lc_data[key]-gal_prop[key][index]
        dist[key] = d
        if(dist_all is None):
            dist_all = d*d
        else:
            dist_all += d*d
    dist_all = np.sqrt(dist_all)
    plt.figure()
    for key in keys:
        slct_fnt = np.isfinite(dist[key])
        bins = np.linspace(np.min(dist[key][slct_fnt]), np.max(dist[key][slct_fnt]), 100)
        h,xbins = np.histogram(dist[key][slct_fnt],bins=bins)
        plt.plot(dtk.bins_avg(xbins),h,label=key)
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('original value - matched value')
    plt.ylabel('count')
    plt.figure()
    slct_fnt = np.isfinite(dist_all)
    bins = np.linspace(np.min(dist_all[slct_fnt]), np.max(dist_all[slct_fnt]), 100)
    h,xbins = np.histogram(dist_all[slct_fnt],bins=bins)
    plt.plot(dtk.bins_avg(xbins),h,label='all',lw=2.0)
    for key in keys:
        slct_fnt = np.isfinite(dist[key])
        h,xbins = np.histogram(dist[key][slct_fnt],bins=xbins)
        plt.plot(dtk.bins_avg(xbins),h,label=key)
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('distance in match')
    plt.ylabel('count')
    return


def plot_differences_obs_color(lc_data, gal_prop, index):
    slct = lc_data['is_cluster_red_sequence']
    keys = ['Mag_r','clr_gr','clr_ri', 'clr_gr_obs', 'clr_ri_obs', 'clr_iz_obs']
    dist = {}
    dist_all = None
    for key in keys:
        d = lc_data[key][slct]-gal_prop[key][index][slct]
        dist[key] = d
        if(dist_all is None):
            dist_all = d*d
        else:
            dist_all += d*d
    dist_all = np.sqrt(dist_all)
    plt.figure()
    for key in keys:
        slct_fnt = np.isfinite(dist[key])
        bins = np.linspace(np.min(dist[key][slct_fnt]), np.max(dist[key][slct_fnt]), 100)
        h,xbins = np.histogram(dist[key][slct_fnt],bins=bins)
        plt.plot(dtk.bins_avg(xbins),h,label=key)
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('original value - matched value')
    plt.ylabel('count')

    plt.figure()
    slct_fnt = np.isfinite(dist_all)
    bins = np.linspace(np.min(dist_all[slct_fnt]), np.max(dist_all[slct_fnt]), 100)
    for key in keys:
        slct_fnt = np.isfinite(dist[key])
        h,xbins = np.histogram(dist[key][slct_fnt],bins=xbins)
        plt.plot(dtk.bins_avg(xbins),h,label=key)
    h,xbins = np.histogram(dist_all[slct_fnt],bins=bins)
    plt.plot(dtk.bins_avg(xbins),h,label='all',lw=2.0)
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('distance in match')
    plt.ylabel('count')
    return


def plot_differences_2d(lc_data, gal_prop,index, x='Mag_r'):
    keys = ['Mag_r','clr_gr','clr_ri','m_star']
    for key in keys:
        # if key == x:
        #     continue
        plt.figure()
        h,xbins,ybins = np.histogram2d(lc_data[x],lc_data[key]-gal_prop[key][index],bins=(100,100))
        plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm = clr.LogNorm())
        plt.ylabel("diff {} (orginal-new)".format(key))
        plt.xlabel(x)
        plt.grid()
    return

 
def plot_side_by_side(lc_data, gal_prop, index, x='Mag_r'):
    keys =  ['Mag_r','clr_gr','clr_ri','m_star']
    for key in keys:
        if key == x:
            continue
        fig,axs = plt.subplots(1,3,sharey=True,sharex=True,figsize=(15,5))
        h,xbins,ybins = np.histogram2d(lc_data[x],lc_data[key],bins=(100,100))
        axs[0].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
        axs[0].grid()
        axs[0].set_title('UMachine + SDSS')
        axs[0].set_xlabel(x)
        axs[0].set_ylabel(key)

        h,xbins,ybins = np.histogram2d(gal_prop[x][index],gal_prop[key][index],bins=(xbins,ybins))
        axs[1].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
        axs[1].grid()
        axs[1].set_title('Matched Galacticus')
        axs[1].set_xlabel(x)
        axs[1].set_ylabel(key)

        h,xbins,ybins = np.histogram2d(gal_prop[x],gal_prop[key],bins=(xbins,ybins))
        axs[2].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
        axs[2].grid()
        axs[2].set_title('Galacticus ')
        axs[2].set_xlabel(x)
        axs[2].set_ylabel(key)
        
    fig,axs = plt.subplots(1,3,sharey=True,sharex=True,figsize=(15,5))
    h,xbins,ybins = np.histogram2d(lc_data['clr_gr'],lc_data['clr_ri'],bins=(100,100))
    axs[0].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    axs[0].grid()
    axs[0].set_title('UMachine + SDSS')
    axs[0].set_xlabel('clr_gr')
    axs[0].set_ylabel('clr_ri')
    
    h,xbins,ybins = np.histogram2d(gal_prop['clr_gr'][index],gal_prop['clr_ri'][index],bins=(xbins,ybins))
    axs[1].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    axs[1].grid()
    axs[1].set_title('Matched Galacticus')
    axs[1].set_xlabel('clr_gr')
    axs[1].set_ylabel('clr_ri')
    
    h,xbins,ybins = np.histogram2d(gal_prop['clr_gr'],gal_prop['clr_ri'],bins=(xbins,ybins))
    axs[2].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    axs[2].grid()
    axs[2].set_title('Galacticus ')
    axs[2].set_xlabel('clr_gr')
    axs[2].set_ylabel('clr_ri')
        
    return


def plot_mag_r(lc_data,gal_prop,index):
    plt.figure()
    slct_fnt_lc = np.isfinite(lc_data['Mag_r'])
    slct_fnt_gp = np.isfinite(gal_prop['Mag_r'][index])
    max_r = np.max((np.max(lc_data['Mag_r'][slct_fnt_lc]),np.max(gal_prop['Mag_r'][index][slct_fnt_gp])))
    min_r = np.min((np.min(lc_data['Mag_r'][slct_fnt_lc]),np.min(gal_prop['Mag_r'][index][slct_fnt_gp])))
    bins = np.linspace(min_r,max_r,100)
    h_lc,_ = np.histogram(lc_data['Mag_r'][slct_fnt_lc],bins=bins)
    h_mg,_ = np.histogram(gal_prop['Mag_r'][index][slct_fnt_gp],bins=bins)
    plt.plot(dtk.bins_avg(bins),h_lc, 'b', label='UMachine-SDSS')
    plt.plot(dtk.bins_avg(bins),h_mg, 'r', label='Matched Glctcs')
    plt.grid()
    plt.xlabel("Mr")
    plt.ylabel('Count')
    plt.legend(loc='best')
    return


def plot_m_star(lc_data,gal_prop,index):
    plt.figure()
    max_r = np.max((np.max(lc_data['m_star']),np.max(gal_prop['m_star'][index])))
    min_r = np.min((np.min(lc_data['m_star']),np.min(gal_prop['m_star'][index])))
    bins = np.linspace(min_r,max_r,100)
    h_lc,_ = np.histogram(lc_data['m_star'],bins=bins)
    h_mg,_ = np.histogram(gal_prop['m_star'][index],bins=bins)
    plt.plot(dtk.bins_avg(bins),h_lc, 'b', label='UMachine-SDSS')
    plt.plot(dtk.bins_avg(bins),h_mg, 'r', label='Matched Glctcs')
    plt.grid()
    plt.xlabel("log10(Stellar Mass)")
    plt.ylabel('Count')
    plt.legend(loc='best')
    return


def plot_single_dist(lc_data,gal_prop,index,key_name,key_label,bins = None):
    plt.figure()
    if bins is None:
        max_r = np.max((np.max(lc_data[key_name]),np.max(gal_prop[key_name][index])))
        min_r = np.min((np.min(lc_data[key_name]),np.min(gal_prop[key_name][index])))
        bins = np.linspace(min_r,max_r,100)
    h_lc,_ = np.histogram(lc_data[key_name],bins=bins)
    h_mg,_ = np.histogram(gal_prop[key_name][index],bins=bins)
    plt.plot(dtk.bins_avg(bins),h_lc, 'b', label='UMachine-SDSS')
    plt.plot(dtk.bins_avg(bins),h_mg, 'r', label='Matched Glctcs')
    plt.grid()
    plt.xlabel(key_label)
    plt.ylabel('Count')
    plt.legend(loc='best')
    return


def plot_clr_mag(lc_data,gal_prop,index,mag_bins,data_key, data_name):
    fig,ax = plt.subplots(1,len(mag_bins),figsize=(15,5))
    # max_gr = np.max((np.max(lc_data[data_key]),np.max(gal_prop[data_key])))
    # min_gr = np.min((np.min(lc_data[data_key]),np.min(gal_prop[data_key])))
    bins = np.linspace(0.0, 1.1, 50)
    for i in range(0, len(mag_bins)):
        if i == 0:
            slct_lc = lc_data['Mag_r']<mag_bins[i]
            slct_mg = gal_prop['Mag_r'][index]<mag_bins[i]
            ax[i].set_title('Mr < {}'.format(mag_bins[i]))
        else:
            slct_lc = (mag_bins[i-1] < lc_data['Mag_r']) & ( lc_data['Mag_r'] < mag_bins[i])
            slct_mg = (mag_bins[i-1] < gal_prop['Mag_r'][index]) & ( gal_prop['Mag_r'][index] < mag_bins[i])
            ax[i].set_title('{} < Mr < {}'.format(mag_bins[i-1], mag_bins[i]))
        h_lc, _ = np.histogram(lc_data[data_key][slct_lc],bins=bins)
        h_mg, _ = np.histogram(gal_prop[data_key][index][slct_mg],bins=bins)
        ax[i].plot(dtk.bins_avg(bins),h_lc,'b', label = 'UMachine-SDSS')
        ax[i].plot(dtk.bins_avg(bins),h_mg,'r', label = 'Matched Glctcs')
        if i ==0:
            ax[i].legend(loc='best', framealpha=0.3)
        ax[i].set_xlabel(data_name)
        ax[i].set_ylabel('Count')
        ax[i].grid()
        

def plot_ri_gr_mag(lc_data, gal_prop, index, mag_bins):
    fig,ax = plt.subplots(1,len(mag_bins),figsize=(15,5))
    for i in range(0, len(mag_bins)):
        if i == 0:
            slct_lc = lc_data['Mag_r']<mag_bins[i]
            slct_mg = gal_prop['Mag_r'][index]<mag_bins[i]
            ax[i].set_title('Mr < {}'.format(mag_bins[i]))
        else:
            slct_lc = (mag_bins[i-1] < lc_data['Mag_r']) & ( lc_data['Mag_r'] < mag_bins[i])
            slct_mg = (mag_bins[i-1] < gal_prop['Mag_r'][index]) & ( gal_prop['Mag_r'][index] < mag_bins[i])
            ax[i].set_title('{} < Mr < {}'.format(mag_bins[i-1], mag_bins[i]))
        # print('{} < Mr < {}'.format(mag_bins[i], mag_bins[i-1]))
        # print(np.average(lc_data['Mag_r'][slct_lc]))
        ax[i].plot(lc_data['clr_gr'][slct_lc], lc_data['clr_ri'][slct_lc],'.b',alpha=0.5,label='UMachine-SDSS',ms=4)
        ax[i].plot(gal_prop['clr_gr'][index][slct_mg], gal_prop['clr_ri'][index][slct_mg],'.r',alpha=0.5,label='Matched Glctcs', ms=4)
        if i ==0:
            ax[i].legend(loc='best', framealpha=0.3)
        ax[i].set_xlabel('g-r color')
        ax[i].set_ylabel('r-i color')
        ax[i].grid()


def plot_clr_z(lc_data, gal_prop, index,clr='clr_gr'):
    fig,axs = plt.subplots(1,3,sharey=True,sharex=True,figsize=(15,5))
    h,xbins,ybins = np.histogram2d(lc_data['redshift'],lc_data[clr],bins=(100,100))
    axs[0].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    axs[0].grid()
    axs[0].set_title('UMachine + SDSS')
    axs[0].set_xlabel('redshift')
    axs[0].set_ylabel(clr)
    h,xbins,ybins = np.histogram2d(gal_prop['redshift'][index],gal_prop[clr][index],bins=(xbins,ybins))
    axs[1].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    axs[1].grid()
    axs[1].set_title('Matched Galacticus')
    axs[1].set_xlabel('redshift')
    axs[1].set_ylabel(clr)
    h,xbins,ybins = np.histogram2d(gal_prop['redshift'],gal_prop[clr],bins=(xbins,ybins))
    axs[2].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    axs[2].grid()
    axs[2].set_title('Galacticus ')
    axs[2].set_xlabel('redshift')
    axs[2].set_ylabel(clr)
    

def plot_gal_prop_dist(gal_props, gal_props_names):
    num = len(gal_props)
    num_list = range(0, num)
    plt.figure()
    bins = np.linspace(5,14,100)
    for i in num_list:
        h,xbins = np.histogram(gal_props[i]['m_star'], bins = bins, normed=True)
        plt.plot(dtk.bins_avg(xbins), h, label=gal_props_names[i])
    plt.grid()
    plt.xlabel('m_star')
    plt.ylabel('count')
    plt.legend(loc='best', framealpha=0.3)

    plt.figure()
    bins = np.linspace(-25,-10, 100)
    for i in num_list:
        h,xbins = np.histogram(gal_props[i]['Mag_r'], bins = bins, normed=True)
        plt.plot(dtk.bins_avg(xbins), h, label=gal_props_names[i])
    plt.grid()
    plt.xlabel('Mag_r')
    plt.ylabel('count')
    plt.legend(loc='best', framealpha=0.3)

    plt.figure()
    bins = np.linspace(-.5,3, 100)
    for i in num_list:
        h,xbins = np.histogram(gal_props[i]['clr_gr'], bins = bins, normed=True)
        plt.plot(dtk.bins_avg(xbins), h, label=gal_props_names[i])
    plt.grid()
    plt.xlabel('clr_gr')
    plt.ylabel('count')
    plt.legend(loc='best', framealpha=0.3)

    plt.figure()
    bins = np.linspace(-.5,3, 100)
    for i in num_list:
        h,xbins = np.histogram(gal_props[i]['clr_ri'], bins = bins, normed=True)
        plt.plot(dtk.bins_avg(xbins), h, label=gal_props_names[i])
    plt.grid()
    plt.xlabel('clr_ri')
    plt.ylabel('count')
    plt.legend(loc='best', framealpha=0.3)
    
    xbins, ybins = (np.linspace(-25,-10,250), np.linspace(-.5,2,250))
    fig, axs = plt.subplots(1,num,figsize=(num*5,5), sharex = True)
    for i in num_list:
        h,xbins,ybins = np.histogram2d(gal_props[i]['Mag_r'], gal_props[i]['clr_gr'],bins=(xbins,ybins))
        axs[i].pcolor(xbins, ybins, h.T, cmap='PuBu', norm=clr.LogNorm())
        axs[i].grid()
        axs[i].set_title(gal_props_names[i])
        axs[i].set_xlabel('Mag_r')
        axs[i].set_ylabel('clr_gr')

    fig, axs = plt.subplots(1,num,figsize=(num*5,5), sharex = True)
    for i in num_list:
        h,xbins,ybins = np.histogram2d(gal_props[i]['Mag_r'], gal_props[i]['clr_ri'],bins=(xbins,ybins))
        axs[i].pcolor(xbins, ybins, h.T, cmap='PuBu', norm=clr.LogNorm())
        axs[i].grid()
        axs[i].set_title(gal_props_names[i])
        axs[i].set_xlabel('Mag_r')
        axs[i].set_ylabel('clr_ri')
    
