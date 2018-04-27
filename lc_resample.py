#!/usr/bin/env python2.7
from __future__ import print_function, division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import pdb
import dtk
import h5py
import time
import sys
import datetime
from astropy.table import Table
from scipy.spatial import cKDTree
from pecZ import pecZ
from astropy.cosmology import WMAP7 as cosmo
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from cosmodc2.black_hole_modeling import monte_carlo_bh_acc_rate, bh_mass_from_bulge_mass, monte_carlo_black_hole_mass
from cosmodc2.size_modeling import mc_size_vs_luminosity_late_type, mc_size_vs_luminosity_early_type
from cosmodc2.sdss_colors import assign_restframe_sdss_gri
import galmatcher

def construct_gal_prop(fname, verbose=False, mask = None, mag_r_cut = False):
    t1 = time.time()
    gal_prop = {}
    hfile = h5py.File(fname,'r')
    hgp = hfile['galaxyProperties']
    m_star = np.log10(hgp['totalMassStellar'].value)
    mag_g = hgp['SDSS_filters/magnitude:SDSS_g:rest:dustAtlas'].value
    mag_r = hgp['SDSS_filters/magnitude:SDSS_r:rest:dustAtlas'].value
    mag_i = hgp['SDSS_filters/magnitude:SDSS_i:rest:dustAtlas'].value
    if mask is None:
        mask = np.ones(mag_r.size,dtype=bool)
    if mag_r_cut:
        mask = (mag_r < -10) & mask
    gal_prop['m_star'] = m_star[mask]
    gal_prop['Mag_r']  = mag_r[mask]
    gal_prop['clr_gr'] = mag_g[mask]-mag_r[mask]
    gal_prop['clr_ri'] = mag_r[mask]-mag_i[mask]
    if verbose:
        print('done loading gal prop. {}'.format(time.time()-t1))
    return gal_prop,mask


def construct_gal_prop_dust_factor(fname, dust_factor, verbose=False, mask = None, mag_r_cut = False):
    t1 = time.time()
    gal_prop = {}
    hfile = h5py.File(fname,'r')
    hgp = hfile['galaxyProperties']
    m_star = np.log10(hgp['totalMassStellar'].value)
    mag_gd = hgp['SDSS_filters/magnitude:SDSS_g:rest:dustAtlas'].value
    mag_rd = hgp['SDSS_filters/magnitude:SDSS_r:rest:dustAtlas'].value
    mag_id = hgp['SDSS_filters/magnitude:SDSS_i:rest:dustAtlas'].value
    mag_gnd = hgp['SDSS_filters/magnitude:SDSS_g:rest'].value
    mag_rnd = hgp['SDSS_filters/magnitude:SDSS_r:rest'].value
    mag_ind = hgp['SDSS_filters/magnitude:SDSS_i:rest'].value
    mag_dgd = mag_gd - mag_gnd
    mag_drd = mag_rd - mag_rnd
    mag_did = mag_id - mag_ind
    mag_g = mag_gd + (dust_factor-1.0)*mag_dgd
    mag_r = mag_rd + (dust_factor-1.0)*mag_drd
    mag_i = mag_id + (dust_factor-1.0)*mag_did
    if verbose:
        print("dust factor: ", dust_factor)
        slct_fnt = np.isfinite(mag_g)
        print(np.sum(mag_g[slct_fnt]-mag_gd[slct_fnt]))
        print("num diff: ", np.sum( (mag_g[slct_fnt]-mag_gd[slct_fnt]) == 0), np.sum(slct_fnt))
        slct_fnt = np.isfinite(mag_r)
        print(np.sum(mag_r[slct_fnt]-mag_rd[slct_fnt]))
        print("num diff: ", np.sum( (mag_r[slct_fnt]-mag_rd[slct_fnt]) == 0), np.sum(slct_fnt))
        slct_fnt = np.isfinite(mag_i)
        print(np.sum(mag_i[slct_fnt]-mag_id[slct_fnt]))
        print("num diff: ", np.sum( (mag_i[slct_fnt]-mag_id[slct_fnt]) == 0), np.sum(slct_fnt))

    if mask is None:
        mask = np.ones(mag_r.size,dtype=bool)
    if mag_r_cut:
        mask = (mag_r < -10) & mask
    gal_prop['m_star'] = m_star[mask]
    gal_prop['Mag_r']  = mag_r[mask]
    gal_prop['clr_gr'] = mag_g[mask]-mag_r[mask]
    gal_prop['clr_ri'] = mag_r[mask]-mag_i[mask]
    gal_prop['Mag_g']  = mag_g[mask]
    gal_prop['Mag_i']  = mag_i[mask]
    #lum debuging
    lum_gd = np.log10(hgp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest:dustAtlas'].value)
    lum_rd = np.log10(hgp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest:dustAtlas'].value)
    # lum_id = np.log10(hgp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest:dustAtlas'].value)
    lum_gnd =np.log10( hgp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest'].value)
    lum_rnd =np.log10( hgp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest'].value)
    # lum_ind =np.log10( hgp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest'].value)
    lum_dgd = lum_gd - lum_gnd
    lum_g = 10**(lum_gd + (dust_factor - 1.0)* lum_dgd)
    
    lum_g = 10**(lum_gd + (dust_factor - 1.0)* lum_dgd)
    print("g lum vs mag: ", np.sum(-2.5*np.log10(lum_gd)!=mag_gd))
    print("dust_mag_g vs dust_lum_g: ", np.sum( -2.5*np.log10(lum_g) != mag_g), mag_g.shape)
    print("dust_mag_g vs dust_lum_g: ", np.nansum( -2.5*np.log10(lum_g) - mag_g))
    if verbose:
        print('done loading gal prop. {}'.format(time.time()-t1))
    return gal_prop,mask


def cat_dics(dics, keys = None):
    new_dic = {}
    if keys is None:
        keys = dics[0].keys()
    for key in keys:
        new_dic[key] = []
        for dic in dics:
            new_dic[key].append(dic[key])
        new_dic[key] = np.concatenate(new_dic[key])
    return new_dic


def construct_gal_prop_redshift(fname,slope_fname,snap_a,target_a,verbose=False,mask=None,mag_r_cut = False,index = None):
    t1 = time.time()
    del_a = target_a - snap_a
    gal_prop = {}
    gal_prop_slope = {}
    hfile = h5py.File(fname,'r')
    hfile_slp = h5py.File(slope_fname,'r')
    hgp = hfile['galaxyProperties']
    hgp_slp = hfile_slp['galaxyProperties']
    m_star =hgp['totalMassStellar'].value
    mag_g = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest:dustAtlas'].value
    mag_r = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest:dustAtlas'].value
    mag_i = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest:dustAtlas'].value
    if mask is None:
        mask = np.ones(mag_r.size,dtype=bool)
    if mag_r_cut:
        mask = (mag_r < -10) & mask
    if index is None:
        index = np.arange(0,np.sum(mask)-1,dtype='i4')
    m_star_slp = hgp_slp['totalMassStellar'].value
    mag_g_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest:dustAtlas'].value
    mag_r_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest:dustAtlas'].value
    mag_i_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest:dustAtlas'].value
    gal_prop['m_star'] = np.log10(m_star[mask][index] + m_star_slp[mask][index]*del_a)
    gal_prop['Mag_r'] = -2.5*np.log10(mag_r[mask][index] + mag_r_slp[mask][index]*del_a)
    gal_prop['clr_gr'] = -2.5*( np.log10(mag_g[mask][index] + mag_g_slp[mask][index]*del_a) - np.log10(mag_r[mask][index] + mag_r_slp[mask][index]*del_a) )
    gal_prop['clr_ri'] = -2.5*( np.log10(mag_r[mask][index] + mag_r_slp[mask][index]*del_a) - np.log10(mag_i[mask][index] + mag_i_slp[mask][index]*del_a) )
    
    #Debug
    # mag_r = hgp['SDSS_filters/magnitude:SDSS_r:rest:dustAtlas'].value
    # lum_r = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest:dustAtlas'].value
    # mag2_r = -2.5*np.log10(
    if verbose:
        print('done loading slope gal prop. {}'.format(time.time()-t1))
    return gal_prop,mask


def clean_up_gal_prop(gal_prop):
    """For each galaxy, if any property is not finite, set all other
    properties to some value (4) that will not be selected by the
    kdtree query.

    """
    print("Cleaning up gal prop: ",end="")
    slct_nfnt =  ~np.isfinite(gal_prop['m_star'])
    for key in gal_prop.keys():
        slct_nfnt = slct_nfnt | ~np.isfinite(gal_prop[key])
        print("bad vals: ", np.sum(slct_nfnt))
    for key in gal_prop.keys():
        gal_prop[key][slct_nfnt] = -4
    return gal_prop


def construct_gal_prop_redshift_dust(fname,slope_fname,snap_a,target_a,verbose=False,mask=None,mag_r_cut = False,index = None, dust_factor = 1.0 ):
    t1 = time.time()
    del_a = target_a - snap_a
    gal_prop = {}
    gal_prop_slope = {}
    hfile = h5py.File(fname,'r')
    hfile_slp = h5py.File(slope_fname,'r')
    hgp = hfile['galaxyProperties']
    hgp_slp = hfile_slp['galaxyProperties']
    m_star_const =hgp['totalMassStellar'].value
    # mag with dust (d= dust)
    mag_gd = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest:dustAtlas'].value
    mag_rd = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest:dustAtlas'].value
    mag_id = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest:dustAtlas'].value
    # mag witout dust (nd = no dust)
    mag_gnd = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest'].value
    mag_rnd = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest'].value
    mag_ind = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest'].value
    # dust extinction
    mag_g_delta = -2.5*np.log10(mag_gd) - -2.5*np.log10(mag_gnd)
    mag_r_delta = -2.5*np.log10(mag_rd) - -2.5*np.log10(mag_rnd)
    mag_i_delta = -2.5*np.log10(mag_id) - -2.5*np.log10(mag_ind)
    if mask is None:
        print("mask is none")
        mask = np.ones(mag_gd.size,dtype=bool)
    if mag_r_cut:
        mask = (mag_rd < -10) & mask
    if index is None:
        print("index is none")
        index = np.arange(0,np.sum(mask)-1,dtype='i4')
    m_star_slp = hgp_slp['totalMassStellar'].value
    # no dust interpolation slopes
    mag_gnd_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest'].value
    mag_rnd_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest'].value
    mag_ind_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest'].value
    # no dust interpolated values
    mag_g = -2.5*np.log10(mag_gnd[mask][index] + mag_gnd_slp[mask][index]*del_a) + mag_g_delta[mask][index]*dust_factor
    mag_r = -2.5*np.log10(mag_rnd[mask][index] + mag_rnd_slp[mask][index]*del_a) + mag_r_delta[mask][index]*dust_factor
    mag_i = -2.5*np.log10(mag_ind[mask][index] + mag_ind_slp[mask][index]*del_a) + mag_i_delta[mask][index]*dust_factor
    m_star = np.log10(m_star_const[mask][index] + m_star_slp[mask][index]*del_a)
    indx = [ 355112, 2744050, 3707872, 5329360, 5963433, 5984882]
    gal_prop['m_star'] = m_star
    gal_prop['Mag_r'] = mag_r
    gal_prop['clr_gr'] = mag_g - mag_r
    gal_prop['clr_ri'] = mag_r - mag_i
    gal_prop['dust_factor'] = np.ones(m_star.size,dtype='f4')*dust_factor
    gal_prop['index'] = np.arange(m_star.size,dtype='i8')
    gal_prop = clean_up_gal_prop(gal_prop)
    if verbose:
        print('done loading slope gal prop. {}'.format(time.time()-t1))
    return gal_prop,mask


def construct_gal_prop_redshift_dust_raw(fname, index, step1, step2, target_a, mask1, mask2, dust_factor = 1.0):
    """Constructs gal_prop using the interpolation scheme from the galacticus
    snapshots and index matching galaxies in step2 to galaxies in step1. 
    """
    h_in_gp1 = h5py.File(fname.replace("${step}", str(step1)), 'r')['galaxyProperties']
    h_in_gp2 = h5py.File(fname.replace("${step}", str(step2)), 'r')['galaxyProperties']
    stepz = dtk.StepZ(sim_name="AlphaQ")
    step1_a = stepz.get_a(step1)
    step2_a = stepz.get_a(step2)
    
    lum_g = get_column_interpolation_dust_raw(
        'SDSS_filters/totalLuminositiesStellar:SDSS_g:rest:dustAtlas',
        h_in_gp1, h_in_gp2, index, mask1, mask2, step1_a, step2_a,
        target_a, dust_factor)
    lum_r = get_column_interpolation_dust_raw(
        'SDSS_filters/totalLuminositiesStellar:SDSS_r:rest:dustAtlas',
        h_in_gp1, h_in_gp2, index, mask1, mask2, step1_a, step2_a,
        target_a, dust_factor)
    lum_i = get_column_interpolation_dust_raw(
        'SDSS_filters/totalLuminositiesStellar:SDSS_i:rest:dustAtlas',
        h_in_gp1, h_in_gp2, index, mask1, mask2, step1_a, step2_a,
        target_a, dust_factor)
    m_star = get_column_interpolation_dust_raw(
        'totalMassStellar',
        h_in_gp1, h_in_gp2, index, mask1, mask2, step1_a, step2_a,
        target_a, dust_factor)
    mag_g = -2.5*np.log10(lum_g)
    mag_r = -2.5*np.log10(lum_r)
    mag_i = -2.5*np.log10(lum_i)
    size = m_star.size
    gal_prop = {}
    gal_prop['m_star'] = np.log10(m_star)
    gal_prop['Mag_r']  = mag_r
    gal_prop['clr_gr'] = mag_g - mag_r
    gal_prop['clr_ri'] = mag_r - mag_i
    gal_prop['dust_factor'] = np.ones(size, dtype='f4')*dust_factor
    gal_prop['index']  = np.arange(size, dtype='i8')
    # print("nan test")
    # print(np.sum(np.isnan(gal_prop['m_star'])))
    # print("not finite test")
    # print(np.sum(~np.isfinite(gal_prop['m_star'])))
    # print(gal_prop['m_star'][np.isnan(gal_prop['m_star'])])
    return gal_prop


def construct_many_gal_prop_redshift_dust(fname, slope_fname, snap_a, target_a_list, verbose=False,mask=None,mag_r_cut = False,index = None, dust_factor_list = [1.0] ):
    t1 = time.time()
    #TODO target_a not defined
    del_a = target_a_list - snap_a
    gal_prop = {}
    gal_prop_slope = {}
    hfile = h5py.File(fname,'r')
    hfile_slp = h5py.File(slope_fname,'r')
    hgp = hfile['galaxyProperties']
    hgp_slp = hfile_slp['galaxyProperties']
    m_star =hgp['totalMassStellar'].value
    # mag with dust (d= dust)
    mag_gd = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest:dustAtlas'].value
    mag_rd = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest:dustAtlas'].value
    mag_id = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest:dustAtlas'].value
    # mag witout dust (nd = no dust)
    mag_gnd = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest'].value
    mag_rnd = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest'].value
    mag_ind = hgp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest'].value
    # dust extinction
    mag_g_delta = -2.5*np.log10(mag_gd) - -2.5*np.log10(mag_gnd)
    mag_r_delta = -2.5*np.log10(mag_rd) - -2.5*np.log10(mag_rnd)
    mag_i_delta = -2.5*np.log10(mag_id) - -2.5*np.log10(mag_ind)
    if mask is None:
        mask = np.ones(mag_rnd.size,dtype=bool)
    if mag_r_cut:
        mask = (mag_rnd < -10) & mask
    if index is None:
        index = np.arange(0,np.sum(mask)-1,dtype='i4')
    m_star_slp = hgp_slp['totalMassStellar'].value
    # no dust interpolation slopes
    mag_gnd_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_g:rest'].value
    mag_rnd_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_r:rest'].value
    mag_ind_slp = hgp_slp['SDSS_filters/totalLuminositiesStellar:SDSS_i:rest'].value
    gal_prop_list = []
    for target_a in target_a_list:
        result_a = []
        del_a = target_a - snap_a
        mag_g_list = []
        mag_r_list = []
        mag_i_list = []
        m_star_list = []
        gal_indx_list = []
        for dust_factor in dust_factor_list:
        # no dust interpolated values
            mag_g_list.append( -2.5*np.log10(mag_gnd[mask][index] + mag_gnd_slp[mask][index]*del_a) + mag_g_delta[mask][index]*dust_factor)
            mag_r_list.append( -2.5*np.log10(mag_rnd[mask][index] + mag_rnd_slp[mask][index]*del_a) + mag_r_delta[mask][index]*dust_factor)
            mag_i_list.append( -2.5*np.log10(mag_ind[mask][index] + mag_ind_slp[mask][index]*del_a) + mag_i_delta[mask][index]*dust_factor)
            m_star_list.append( np.log10(m_star[mask][index] + m_star_slp[mask][index]*del_a))
            gal_indx_list.append(np.arange(m_star[-1].size))
        m_star = np.concatenate(m_star_list)
        mag_g = np.concatenate(mag_g_list)
        mag_r = np.concatenate(mag_r_list)
        mag_i = np.concatenate(mag_i_list)
        gal_indx = np.concatenate(gal_indx_list)
        gal_prop['m_star'] = m_star
        gal_prop['Mag_r'] = mag_r
        gal_prop['clr_gr'] = mag_g - mag_r
        gal_prop['clr_ri'] = mag_r - mag_i
        gal_prop['dust_factor'] = np.ones(m_star.size,dtype='f4')*dust_factor
        gal_prop['index'] = gal_indx
        gal_prop_list.append(gal_prop)
    if verbose:
        print('done loading slope gal prop. {}'.format(time.time()-t1))
    return gal_prop_list,mask


def construct_lc_data(fname,verbose = False, recolor=False):
    t1 = time.time()
    lc_data = {}
    #tbl = Table.read(fname,path='data')
    hfile = h5py.File(fname,'r')
    # lc_data['m_star'] = np.log10(sdss['obs_sm'].quantity)
    # lc_data['Mag_r']  = sdss['rmag'].quantity
    # lc_data['clr_gr'] = sdss['sdss_petrosian_gr'].quantity
    # lc_data['clr_ri'] = sdss['sdss_petrosian_ri'].quantity
    lc_data['m_star'] = np.log10(hfile['obs_sm'].value)
    lc_data['Mag_r'] = hfile['restframe_extincted_sdss_abs_magr'].value
    lc_data['clr_gr'] = hfile['restframe_extincted_sdss_gr'].value
    lc_data['clr_ri'] = hfile['restframe_extincted_sdss_ri'].value
    lc_data['redshift'] = hfile['redshift'].value
    lc_data['sfr_percentile'] = hfile['sfr_percentile'].value
    if recolor:
        print(hfile.keys())
        upid_mock = hfile['upid'].value
        mstar_mock = hfile['obs_sm'].value
        sfr_percentile_mock = hfile['sfr_percentile'].value
        host_halo_mvir_mock = hfile['host_halo_mvir'].value
        redshift_mock = lc_data['redshift']
        a,b,c = assign_restframe_sdss_gri(upid_mock, mstar_mock, sfr_percentile_mock,
                                          host_halo_mvir_mock, redshift_mock)
        # plt.figure()
        # h,xbins,ybins = np.histogram2d(lc_data['Mag_r'], a, bins=250)
        # plt.pcolor(xbins,ybins, h.T, cmap='PuBu', norm =clr.LogNorm())
        # plt.grid()

        # plt.figure()
        # h,xbins,ybins = np.histogram2d(lc_data['clr_gr'], b, bins=250)
        # plt.pcolor(xbins,ybins, h.T, cmap='PuBu', norm =clr.LogNorm())
        # plt.grid()

        # plt.figure()
        # h,xbins,ybins = np.histogram2d(lc_data['clr_ri'], c, bins=250)
        # plt.pcolor(xbins,ybins, h.T, cmap='PuBu', norm =clr.LogNorm())
        # plt.grid()
        # plt.show()
        lc_data['Mag_r'] = a 
        lc_data['clr_gr'] = b
        lc_data['clr_ri'] = c
        #lc_data['Mag_r'], lc_data['clr_gr'], lc_data['clr_ri'] = [a,b,c]
        
    if verbose:
        print('done loading lc data. {}'.format(time.time()-t1))
    return lc_data


def dic_select(dic, slct):
    new_dic = {}
    for key in dic.keys():
        new_dic[key] = dic[key][slct]
    return new_dic


def select_by_index(data,index):
    new_data = {}
    for key in data.keys():
        new_data[key] = data[key][index]
    return new_data


def resample_index(lc_data, gal_prop, ignore_mstar = False, nnk = 10, verbose = False):
    if verbose:
        t1 = time.time()
        print("Starting kdtree resampling")
    # print("lc_data size: {}".format(lc_data['m_star'].size))
    # print("gal_prop size: {}".format(gal_prop['m_star'].size))
    m_star = lc_data['m_star']
    mag_r  = lc_data['Mag_r']
    clr_gr = lc_data['clr_gr']
    clr_ri = lc_data['clr_ri']
    print("ignore_mstar: ", ignore_mstar)
    if ignore_mstar:
        print("Ignoring Mstar!")
        lc_mat = np.stack((mag_r,clr_gr,clr_ri),axis=1)
        gal_mat = np.stack((gal_prop['Mag_r'],
                            gal_prop['clr_gr'],
                            gal_prop['clr_ri']),axis=1)
    else:
        lc_mat = np.stack((m_star,mag_r,clr_gr,clr_ri),axis=1)
        gal_mat = np.stack((gal_prop['m_star'],
                            gal_prop['Mag_r'],
                            gal_prop['clr_gr'],
                            gal_prop['clr_ri']),axis=1)

    if verbose:
        t2 = time.time()
        print('\tdone formating data. {}'.format(t2-t1))

    ckdtree = cKDTree(gal_mat, balanced_tree = False, compact_nodes = False)
    if verbose:
        t3 = time.time()
        print('\tdone making tree. {}'.format(t3-t2))
    dist, index = ckdtree.query(lc_mat, nnk, n_jobs=10)
    if verbose:
        t4= time.time()
        print('\tdone querying. {}'.format(t4-t3))
    if nnk > 1:
        rand = np.random.randint(nnk,size=dist.shape[0])
        aa = np.arange(dist.shape[0])
        #dist = dist[aa,rand]
        index = index[aa,rand]
    # print("distance: ",dist)
    # print("gal_mat")
    # print("{:20} {:20} {:20} {:20}".format("mstar", "MagR","g-r", "r-i"))
    # for i in range(0,20):
    #     print("{:20} {:20} {:20} {:20}".format(gal_mat[i,0], gal_mat[i,1], gal_mat[i,2], gal_mat[i,3]))
    # print("lc_mat")
    # print("{:20} {:20} {:20} {:20}".format("mstar", "MagR","g-r", "r-i"))
    # for i in range(0,20):
    #     print("{:20} {:20} {:20} {:20}".format(lc_mat[i,0], lc_mat[i,1], lc_mat[i,2], lc_mat[i,3]))
    # print("index: ")
    # for i in range(0,20):
    #     print(index[i])
    return index
                                

def get_keys(hgroup):
    keys = []
    def _collect_keys(name, obj):
        if isinstance(obj, h5py.Dataset): 
            keys.append(name)
    hgroup.visititems(_collect_keys)
    return keys


copy_avoids = ('x','y','z','vx','vy','vz', 'peculiarVelocity','galaxyID','redshift',
               'redshiftHubble','placementType','isCentral','hostIndex', 
               'blackHoleAccretionRate','blackHoleMass')
copy_avoids_ptrn = ('hostHalo','magnitude','ageStatistics','Radius','Axis','Ellipticity','positionAngle')
no_slope_var = ('x','y','z','vx','vy','vz', 'peculiarVelocity','galaxyID','redshift','redshiftHubble','inclination','positionAngle')
no_slope_ptrn  =('morphology','hostHalo','infall')

def to_copy(key, short, supershort):
    if short:
        if "LSST" in key or "SED" in key or "other" in key or "Lines" in key:
            print("short: ",key)
            return False
    if supershort:
        if "SDSS" not in key and "total" not in key and ":rest" not in key and "MassStellar" not in key and "infallIndex" != key and "inclination" not in key:
            print("supershort: ", key)
            return False
    if any([ ca == key for ca in copy_avoids]) or any([ cap in key for cap in copy_avoids_ptrn ]):
        print("{} isn't copied".format(key))
        return False
    return True


def copy_columns(input_fname, output_fname, index, 
                 verbose = False, mask = None, 
                 short = False, supershort = False, step = -1):
    h_in = h5py.File(input_fname,'r')
    h_out = h5py.File(output_fname,'w')
    h_in_gp = h_in['galaxyProperties']
    h_out_gp = h_out.create_group('galaxyProperties')
    keys = get_keys(h_in_gp)
    for i in range(0,len(keys)):
        key = keys[i]
        if "LSST" in key or "SED" in key or "other" in key or "Lines" in key or "morphology" in key:
            if short:
                continue
        if "SDSS" not in key and "total" not in key and "MassStellar" not in key and "infallIndex" != key:
            if supershort:
                continue
        if any([ ca == key for ca in copy_avoids]) or any([ cap in key for cap in copy_avoids_ptrn ]):
            print("{} isn't copied".format(key))
            continue
        print('{}/{},{} {} {}'.format(i,len(keys),step,float(i)/float(len(keys)), key))
        data = h_in_gp[key].value
        if mask is not None:
            data = data[mask]
        h_out_gp[key]=data[index]
    return


def copy_columns_dust(input_fname, output_fname, raw_index, dust_factor, 
                      verbose = False,mask = None,
                      short = False, supershort = False, step = -1, dust_factors = 1.0):
    h_in = h5py.File(input_fname,'r')
    h_out = h5py.File(output_fname,'w')
    h_in_gp = h_in['galaxyProperties']
    h_out_gp = h_out.create_group('galaxyProperties')
    h_out_gp['dustFactor'] = dust_factor
    keys = get_keys(h_in_gp)
    for i in range(0,len(keys)):
        t1 = time.time()
        key = keys[i]
        if "LSST" in key or "SED" in key or "other" in key or "Lines" in key or "morphology" in key:
            if short:
                continue
        if "SDSS" not in key and "total" not in key and ":rest" not in key and "MassStellar" not in key and "infallIndex" != key:
            if supershort:
                continue
        if any([ ca == key for ca in copy_avoids]) or any([ cap in key for cap in copy_avoids_ptrn ]):
            print("{} isn't copied".format(key))
            continue
        print('{}/{},{} {} {}'.format(i,len(keys),step,float(i)/float(len(keys)), key), end='')
        if ":dustAtlas" in key:
            org_data = h_in_gp[key].value[mask][raw_index]
            nd_key = key.replace(":dustAtlas","")
            dust_data =  np.log10(h_in_gp[key].value[mask][raw_index])
            no_dust_data = np.log10(h_in_gp[nd_key].value[mask][raw_index])
            delta_dust = dust_data - no_dust_data
            data = 10**(dust_data + (dust_factor -1.0)*delta_dust)
            # nd_key = key.replace(":dustAtlas","")
            # dust_data =  np.log10(h_in_gp[key].value[mask][index])
            # no_dust_data = np.log10(h_in_gp[nd_key].value[mask][index])
            # delta_dust = dust_data - no_dust_data 
            # data = 10**(dust_data + (dust_factor -1.0 ) * delta_dust)
            # print(h_in_gp[key][mask][index])
            # print(data)
            # print("sum diff: ", np.sum(h_in_gp[key][mask][index] - data))
        else:
            data = h_in_gp[key].value
            if mask is not None:
                data = data[mask]
            data = data[raw_index]
        h_out_gp[key]=data
        print("\n\ttime: {}".format(time.time()-t1))
    return
    

def copy_columns_slope(input_fname, input_slope_fname, 
                       output_fname, index,  
                       input_a, lc_a, 
                       verbose = False, mask = None,
                       short = False, supershort = False, step = -1):
    # lc_a = 1.0/(1.0+lc_redshift)
    # input_a = 1.0/(1.0 + input_redshift)
    del_a = lc_a-input_a
    h_in = h5py.File(input_fname,'r')
    h_in_slope = h5py.File(input_slope_fname,'r')
    h_out = h5py.File(output_fname,'w')
    h_in_gp = h_in['galaxyProperties']
    h_in_slope_gp = h_in_slope['galaxyProperties']
    h_out_gp = h_out.create_group('galaxyProperties')
    keys = get_keys(h_in_gp)
    max_float = np.finfo(np.float32).max #The max float size
    for i in range(0,len(keys)):
        key = keys[i]
        if "LSST" in key or "SED" in key or "other" in key or "Lines" in key or "morphology" in key:
            if short:
                continue
        if "SDSS" not in key and "total" not in key and ":rest" not in key and "MassStellar" not in key and 'infallIndex' != key:
            if supershort:
                continue
        if any([ ca == key for ca in copy_avoids]) or any([ cap in key for cap in copy_avoids_ptrn ]):
            print("{} isn't copied".format(key))
            continue
   
        if verbose:
            print('{}/{},{} {} {}'.format(i,len(keys),step,float(i)/float(len(keys)), key))
        data = h_in_gp[key].value
        slope = h_in_slope_gp[key].value
        if mask is not None:
            data = data[mask]
            slope = slope[mask]
        no_slope = key in no_slope_var or any(s in key for s in no_slope_ptrn)
        if (data.dtype == np.float64 or data.dtype == np.float32) and not no_slope:
            #print("\tslope")
            new_data = data[index] + slope[index]*del_a
        else:
            #print("\tno slope")
            new_data = data[index]
        slct_finite = np.isfinite(new_data)
        if(new_data.dtype == np.float64 and np.sum(new_data[slct_finite]>max_float) == 0):
            h_out_gp[key]= new_data.astype(np.float32)
        else:
            h_out_gp[key] = new_data
    return


def get_column_slope_dust(key, h_in_gp, h_in_slope_gp, index, del_a, mask = None, dust_factors = 1.0):
    if ":dustAtlas" in key:
        print("\thas dust")
        key_nd = key.replace(":dustAtlas","")
        data = h_in_gp[key_nd].value
        slope = h_in_slope_gp[key_nd].value
        # multiplicitive effect of dust
        data_dust = h_in_gp[key].value
        dust_effect = data_dust/data 
        if mask is not None:
            dust_effect = dust_effect[mask]
    else:
        data = h_in_gp[key].value
        slope = h_in_slope_gp[key].value
    if mask is not None:
        data = data[mask]
        slope = slope[mask]

    no_slope = key in no_slope_var or any(s in key for s in no_slope_ptrn)
    if (data.dtype == np.float64 or data.dtype == np.float32) and not no_slope:
        print("\t float & slope")
        if ":dustAtlas" in key:
            print("\t\t has dustAtlas")
            # after interpolating on the undusted luminosity, apply the effect of dust
            new_data = (data[index] + slope[index]*del_a)*(dust_effect[index]**dust_factors)
        else:
            print("\t\t does not have dustAtlas")
            new_data = data[index] + slope[index]*del_a
    else:
        new_data = data[index]
    return new_data


def get_column_interpolation_dust_raw(key, h_in_gp1, h_in_gp2, index, mask1, mask2, step1_a, step2_a, target_a, dust_factors, kdtree_index=None):
    """This function returns the interpolated quantity between two
    timesteps, from step1 to step2. Some galaxies are masked out: Any
    galaxy that doesn't pass the mask in step1 (mask1), any galaxy
    that doesn't a decendent in step2, or any galaxy that whose
    descendent doesn't pass the step2 mask (mask2).

    """
    print("\tLoading key: {}".format(key), end="")
    #print("dust_factors: ", dust_factors)
    t1 = time.time()
    step_del_a = step2_a - step1_a
    target_del_a = target_a - step1_a
    # The masking all galaxies that fail galmatcher's requirements at
    # step1, galaxies that don't have a descndent, or if the
    # descendent galaxy at step2 doesn't pass galmatcher requirements.
    mask_tot = mask1 & (index != -1) & mask2[index]
    if ":dustAtlas" in key:
        key_no_dust = key.replace(":dustAtlas","")
        val1_no_dust = h_in_gp1[key_no_dust].value[mask_tot]
        val1_dust = h_in_gp1[key].value[mask_tot]
        val2_no_dust = h_in_gp2[key].value[index][mask_tot]
        dust_effect = val1_dust/val1_no_dust
        #dust_effect = val1_dust - val1_no_dust
        slope = (val2_no_dust - val1_no_dust)/step_del_a
        if kdtree_index is None:
            val_out = (val1_no_dust + slope*target_del_a)*(dust_effect**dust_factors)
            #val_out = 10**(val1_no_dust + slope*target_del_a + dust_effect*dust_factor)
        else:
            val_out = (val1_no_dust[kdtree_index] + slope[kdtree_index]*target_del_a)*(dust_effect[kdtree_index]**dust_factors)
            #val_out = 10**(val1_no_dust[kdtree_index] + slope[kdtree_index]*target_del_a + dust_effect[kdtree_index]*dust_factor)
    else:
        val1_data = h_in_gp1[key].value[mask_tot]
        val2_data = h_in_gp2[key].value[index][mask_tot]
        slope = (val2_data - val1_data)/step_del_a
        if kdtree_index is None:
            val_out = val1_data + slope*target_del_a
        else:
            val_out = val1_data[kdtree_index] + slope[kdtree_index]*target_del_a
    if(val_out.dtype == np.float64):
        val_out = val_out.astype(np.float32)
    print("\t\tread + format time: {}".format(time.time()-t1))
    return val_out


def get_func_interpolation_dust_raw(key, h_in_gp1, h_in_gp2, index, mask1, mask2, step1_a, step2_a, target_a, dust_factors=1.0):
    """Returns the constants required to reconstruct the column at any redshift between 
    the interpolation steps. The return values are val0, slope0 and dust_effect0 for the
    function val(z, dust) = (val0 + slope0*del_a)*(dust_effect0**dust_factor)
    """
    step_del_a = step2_a - step1_a
    target_del_a = target_a - step1_a
    # The masking all galaxies that fail galmatcher's requirements at
    # step1, galaxies that don't have a descndent, or if the
    # descendent galaxy at step2 doesn't pass galmatcher requirements.
    mask_tot = mask1 & (index != -1) & mask2[index]
    if ":dustAtlas" in key:
        key_no_dust = key.replace(":dustAtlas","")
        val1_no_dust = h_in_gp1[key_no_dust].value[mask_tot]
        val1_dust = h_in_gp1[key].value[mask_tot]
        val2_no_dust = h_in_gp2[key].value[index][mask_tot]
        dust_effect0 = val1_dust/val1_no_dust
        val0 = val1_no_dust
        slope0 = (val2_no_dust - val1_no_dust)/step_del_a
        #val_out = (val1_no_dust + slope*target_del_a)*(dust_effect**dust_effect)
    else:
        val1_data = h_in_gp1[key].value[mask_tot]
        val2_data = h_in_gp2[key].value[index][mask_tot]
        slope0 = (val2_data - val1_data)/step_del_a
        val0 = val1_data
        dust_effect0 = np.ones(val0.size, dtype='f4')
        #val_out = val1_data + slope*target_del_a
    if(val0.dtype == np.float64):
        val0 = val0.astype(np.float32)
    if(slope0.dtype == np.float64):
        slope0 = slope0.astype(np.float32)
    if(dust_effect0.dtype == np.float64):
        dust_effect0 = dust_effect0.astype(np.float32)
    return val0, slope0, dust_effect0
    

def copy_columns_slope_dust(input_fname, input_slope_fname, 
                            output_fname, index,  
                            input_a, lc_a, 
                            verbose = False, mask = None,
                            short = False, supershort = False, 
                            step = -1, dust_factors = 1.0):
    # lc_a = 1.0/(1.0+lc_redshift)
    # input_a = 1.0/(1.0 + input_redshift)
    del_a = lc_a-input_a
    print("del_a: ", del_a)
    h_in = h5py.File(input_fname,'r')
    h_in_slope = h5py.File(input_slope_fname,'r')
    h_out = h5py.File(output_fname,'w')
    h_in_gp = h_in['galaxyProperties']
    h_in_slope_gp = h_in_slope['galaxyProperties']
    h_out_gp = h_out.create_group('galaxyProperties')
    keys = get_keys(h_in_gp)
    max_float = np.finfo(np.float32).max #The max float size
    for i in range(0,len(keys)):
        key = keys[i]
        if not to_copy(key, short, supershort):
            continue
        if verbose:
            print('{}/{} [{}] {}'.format(i,len(keys),step, key))
        new_data = get_column_slope_dust(key, h_in_gp, h_in_slope_gp, index, del_a, mask, dust_factors)
        slct_finite = np.isfinite(new_data)
        #If the data is a double, record it as a float to save on disk space
        if(new_data.dtype == np.float64 and np.sum(new_data[slct_finite]>max_float) == 0):
            h_out_gp[key]= new_data.astype(np.float32)
        else:
            h_out_gp[key] = new_data
    return


def copy_columns_interpolation_dust_raw(input_fname, output_fname,
                                        kdtree_index, step1, step2,
                                        step1_a, step2_a, mask1, mask2, 
                                        index_2to1, lc_a, 
                                        verbose = False, 
                                        short = False, supershort = False, 
                                        step = -1, dust_factors = 1.0):
    print("===================================")
    print("copy columns interpolation dust raw")
    # lc_a = 1.0/(1.0+lc_redshift)
    # input_a = 1.0/(1.0 + input_redshift)
    del_a = lc_a-step1_a
    print("del_a: ", del_a)
    h_out = h5py.File(output_fname,'w')
    h_out_gp = h_out.create_group('galaxyProperties')
    h_out_gp['dustFactor'] = dust_factors
    h_in_gp1 = h5py.File(input_fname.replace("${step}",str(step1)),'r')['galaxyProperties']
    h_in_gp2 = h5py.File(input_fname.replace("${step}",str(step2)),'r')['galaxyProperties']

    keys = get_keys(h_in_gp1)
    max_float = np.finfo(np.float32).max #The max float size
    for i in range(0,len(keys)):
        t1 = time.time()
        key = keys[i]
        if not to_copy(key, short, supershort):
            continue
        if verbose:
            print('{}/{} [{}] {}'.format(i,len(keys),step, key))
        new_data = get_column_interpolation_dust_raw(
            key, h_in_gp1, h_in_gp2, index_2to1, mask1, mask2, step1_a, step2_a, lc_a, dust_factors, 
            kdtree_index = kdtree_index)
        slct_finite = np.isfinite(new_data)
        #If the data is a double, record it as a float to save on disk space
        if(new_data.dtype == np.float64 and np.sum(new_data[slct_finite]>max_float) == 0):
            h_out_gp[key]= new_data.astype(np.float32)
        else:
            h_out_gp[key] = new_data
        print("\t\tDone writing. read+format+write: {}".format(time.time()-t1))
    return


def overwrite_columns(input_fname, output_fname, ignore_mstar = False, verbose=False):
    t1 = time.time()
    if verbose:
        print("Overwriting columns.")
        #sdss = Table.read(input_fname,path='data')
    h_in = h5py.File(input_fname,'r')
    #redshift = np.ones(sdss['x'].quantity.size)*0.1
    h_out = h5py.File(output_fname, 'a')
    h_out_gp = h_out['galaxyProperties']
    t2 = time.time()
    if verbose:
        print("\t done reading in data", t2-t1)
    #xyz,v(xyz)
    x = h_in['x'].value
    y = h_in['y'].value
    z = h_in['z'].value
    vx = h_in['vx'].value
    vy = h_in['vy'].value
    vz = h_in['vz'].value
    redshift  =h_in['redshift'].value
    h_out_gp['x']=x
    h_out_gp['y']=y
    h_out_gp['z']=z
    h_out_gp['vx']=vx
    h_out_gp['vy']=vy
    h_out_gp['vz']=vz
    h_out_gp['lightcone_rotation'] = h_in['lightcone_rotation'].value
    h_out_gp['lightcone_replication'] = h_in['lightcone_replication'].value
    if ignore_mstar:
        print("Ignoring M* in stellar mass assignment!")
        # m*_delta = M*_new/M*_old
        mstar_delta = h_in['obs_sm'].value/h_out_gp['totalMassStellar'].value
        h_out_gp['totalMassStellar'][:] = h_out_gp['totalMassStellar'].value*mstar_delta
        h_out_gp['diskMassStellar'][:] = h_out_gp['diskMassStellar'].value*mstar_delta
        h_out_gp['spheroidMassStellar'][:] = h_out_gp['spheroidMassStellar'].value*mstar_delta
    t3 = time.time()
    if verbose:
        print("\t done overwriting xyz, v_(xyz)",t3-t2)
    #peculiar velocity
    _,z_obs,v_pec,_,_,_,_ = pecZ(x,y,z,vx,vy,vz,redshift)
    h_out_gp['peculiarVelocity'] = v_pec
    #obs mag
    #Calculate the oringal redshift 
    stepz = dtk.StepZ(200,0,500)
    zs = np.linspace(0,1.5,1000)
    z_to_dl = interp1d(zs,cosmo.luminosity_distance(zs))
    dl = z_to_dl(redshift)
    adjust_mag = -2.5*np.log10(1.0+redshift)+5*np.log10(dl)+25.0
    keys = get_keys(h_out_gp)
    t4 = time.time()
    for key in keys:
        # Calculating new observer frame magnitudes
        if("totalLuminositiesStellar" in key and  ":observed" in key and ("SDSS" in key or "LSST" in key)):
            new_key = key.replace("totalLuminositiesStellar",'magnitude',1)
            print("making: "+new_key+" from "+key)
            h_out_gp[new_key]=adjust_mag -2.5*np.log10(h_out_gp[key].value)
        # Calculating new rest frame magnitudes
        if("totalLuminositiesStellar" in key and  ":rest" in key and ("SDSS" in key or "LSST" in key)):
            new_key = key.replace("totalLuminositiesStellar","magnitude",1)
            print("making: "+new_key+" from "+key)
            h_out_gp[new_key]=-2.5*np.log10(h_out_gp[key].value)

    # redshift_org = h_out_gp['redshiftHubble'].value

    # dl_org = z_to_dl(redshift_org)
    # adjust_org = -2.5*np.log10(1.0+redshift_org) + 5*np.log10(dl_org)
    # dl_new = z_to_dl(redshift_new)
    # adjust_new = -2.5*np.log10(1.0+redshift_new) + 5*np.log10(dl_new)
    # t4 = time.time()
    # if verbose:
    #     print("\t done calculating mag adjustment")
    # adjust_tot = -adjust_org + adjust_new
    # keys = get_keys(h_out_gp)
    # for key in keys:
    #     if "magnitude" in key and "observed" in key:
    #         print("\t\tadjusting mag: {}".format(key))
    #         h_out_gp[key][:] = h_out_gp[key].value + adjust_tot
    t5 = time.time()
    if verbose:
        print("\t done rewriting mags",t5-t4)
    #redshift
    h_out_gp['redshift'] = z_obs
    h_out_gp['redshiftHubble'] = redshift
    h_out_gp['ra_true'] = h_in['ra'].value
    h_out_gp['dec_true'] = h_in['dec'].value
    h_out_gp['ra'] = h_in['ra_lensed'].value
    h_out_gp['dec'] = h_in['dec_lensed'].value
    h_out_gp['shear1'] = h_in['shear1'].value
    h_out_gp['shear2'] = h_in['shear2'].value
    h_out_gp['magnification'] = h_in['magnification'].value
    h_out_gp['convergence'] = h_in['convergence'].value
    central = (h_in['host_centric_x'].value ==0) & (h_in['host_centric_y'].value ==0) & (h_in['host_centric_z'].value == 0)
    h_out_gp['isCentral'] = central
    h_out_gp['hostHaloTag'] = h_in['target_halo_id'].value
    h_out_gp['hostHaloMass'] = h_in['target_halo_mass'].value
    unq, indx, cnt = np.unique(h_out_gp['infallIndex'].value, return_inverse=True, return_counts = True)
    h_out_gp['NumberSelected'] = cnt[indx]
    tf = time.time()
    if verbose:
        print("\tDone overwrite columns", tf-t1)


def swap(slct, x1, x2):
    xa = x1[slct]
    x1[slct] = x2[slct]
    x2[slct]=xa


def rotate_host_halo(rot, x,y,z):
    """"
    From the documentation:
    
    // 0 = no rotation
    // 1 = swap(x, y) rotation applied
    // 2 = swap(y, z) rotation applied
    // 3 = swap(x, y) then swap(y, z) rotations applied 
    // 4 = swap(z, x) rotation applied
    // 5 = swap(x, y) then swap(z, x) rotations applied 

    // There are also two other possibilities:
    // 6 = swap(y, z) then swap(z, x) rotations applied 
    // 7 = swap(x, y), then swap(y, z), then swap(z, x) rotations applied 
    // which are equivalent to rotations 3 and 2 above, respectively

    """
    slct1 = rot == 1
    swap(slct1, x, y)

    slct2 = (rot == 2) | (rot == 7)
    swap(slct2, y, z)

    slct3 = (rot == 3) | (rot == 6)
    swap(slct3, x, y )
    swap(slct3, y, z)

    slct4 = rot == 4
    swap(slct4, z, x)

    slct5 = rot == 5  
    swap(slct5, x, y)
    swap(slct5, z, x)
    return


def overwrite_host_halo(output_fname, sod_loc, halo_shape_loc, halo_shape_red_loc, verbose=False):
    hgroup = h5py.File(output_fname,'r+')['galaxyProperties']

    halo_tag = hgroup['hostHaloTag'].value
    size = halo_tag.size
    halo_rot = hgroup['lightcone_rotation'].value
    # fof_tag =  dtk.read_gio(fof_loc,'fof_halo_tag')
    # fof_mass = dtk.read_gio(fof_loc,'fof_mass')
    # fof_srt = np.argsort(fof_tag)

    sod_cat_tag = dtk.gio_read(sod_loc,'fof_halo_tag')
    sod_cat_mass = dtk.gio_read(sod_loc,'sod_halo_mass')
    sod_cat_srt = np.argsort(sod_cat_tag)
    sod_mass = -1*np.ones(halo_tag.size)
    indx = dtk.search_sorted(sod_cat_mass,halo_tag,sorter=sod_cat_srt)
    slct = indx != -1
    sod_mass[slct] = sod_cat_mass[indx[slct]]
    hgroup['hostHaloSODMass']=sod_mass

    print("Num of galaxies: ", halo_tag.size)
    eg_cat_htag = dtk.gio_read(halo_shape_loc,'halo_id')
    srt = np.argsort(eg_cat_htag)
    indx = dtk.search_sorted(eg_cat_htag,halo_tag,sorter=srt)
    slct_indx = indx != -1
    indx_slct = indx[slct_indx]
    print("num selected: ",np.sum(slct_indx))
    

    eg_cat_eg1 = np.zeros(size)
    eg_cat_eg2 = np.zeros(size)
    eg_cat_eg3 = np.zeros(size)
    eg_cat_eg1_x =np.zeros(size)
    eg_cat_eg1_y =np.zeros(size)
    eg_cat_eg1_z =np.zeros(size)
    eg_cat_eg2_x =np.zeros(size)
    eg_cat_eg2_y =np.zeros(size)
    eg_cat_eg2_z =np.zeros(size)
    eg_cat_eg3_x = np.zeros(size)
    eg_cat_eg3_y =np.zeros(size)
    eg_cat_eg3_z =np.zeros(size)

    eg_cat_eg1[slct_indx] = dtk.gio_read(halo_shape_loc,'eval1')[indx_slct]
    eg_cat_eg2[slct_indx] = dtk.gio_read(halo_shape_loc,'eval2')[indx_slct]
    eg_cat_eg3[slct_indx] = dtk.gio_read(halo_shape_loc,'eval3')[indx_slct]
    eg_cat_eg1_x[slct_indx] = dtk.gio_read(halo_shape_loc,'evec1x')[indx_slct]
    eg_cat_eg1_y[slct_indx] = dtk.gio_read(halo_shape_loc,'evec1y')[indx_slct]
    eg_cat_eg1_z[slct_indx] = dtk.gio_read(halo_shape_loc,'evec1z')[indx_slct]
    eg_cat_eg2_x[slct_indx] = dtk.gio_read(halo_shape_loc,'evec2x')[indx_slct]
    eg_cat_eg2_y[slct_indx] = dtk.gio_read(halo_shape_loc,'evec2y')[indx_slct]
    eg_cat_eg2_z[slct_indx] = dtk.gio_read(halo_shape_loc,'evec2z')[indx_slct]
    eg_cat_eg3_x[slct_indx] = dtk.gio_read(halo_shape_loc,'evec3x')[indx_slct]
    eg_cat_eg3_y[slct_indx] = dtk.gio_read(halo_shape_loc,'evec3y')[indx_slct]
    eg_cat_eg3_z[slct_indx] = dtk.gio_read(halo_shape_loc,'evec3z')[indx_slct]
    rotate_host_halo(halo_rot, eg_cat_eg1_x, eg_cat_eg1_y, eg_cat_eg1_z)
    rotate_host_halo(halo_rot, eg_cat_eg2_x, eg_cat_eg2_y, eg_cat_eg2_z)
    rotate_host_halo(halo_rot, eg_cat_eg3_x, eg_cat_eg3_y, eg_cat_eg3_z)
    hgroup['hostHaloEigenValue1'] = eg_cat_eg1
    hgroup['hostHaloEigenValue2'] = eg_cat_eg2
    hgroup['hostHaloEigenValue3'] = eg_cat_eg3
    hgroup['hostHaloEigenVector1X'] = eg_cat_eg1_x
    hgroup['hostHaloEigenVector1Y'] = eg_cat_eg1_y
    hgroup['hostHaloEigenVector1Z'] = eg_cat_eg1_z
    hgroup['hostHaloEigenVector2X'] = eg_cat_eg2_x
    hgroup['hostHaloEigenVector2Y'] = eg_cat_eg2_y
    hgroup['hostHaloEigenVector2Z'] = eg_cat_eg2_z
    hgroup['hostHaloEigenVector3X'] = eg_cat_eg3_x
    hgroup['hostHaloEigenVector3Y'] = eg_cat_eg3_y
    hgroup['hostHaloEigenVector3Z'] = eg_cat_eg3_z


    eg_cat_htag = dtk.gio_read(halo_shape_red_step_loc,'halo_id')[indx_slct]
    eg_cat_eg1[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'eval1')[indx_slct]
    eg_cat_eg2[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'eval2')[indx_slct]
    eg_cat_eg3[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'eval3')[indx_slct]
    eg_cat_eg1_x[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec1x')[indx_slct]
    eg_cat_eg1_y[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec1y')[indx_slct]
    eg_cat_eg1_z[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec1z')[indx_slct]
    eg_cat_eg2_x[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec2x')[indx_slct]
    eg_cat_eg2_y[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec2y')[indx_slct]
    eg_cat_eg2_z[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec2z')[indx_slct]
    eg_cat_eg3_x[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec3x')[indx_slct]
    eg_cat_eg3_y[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec3y')[indx_slct]
    eg_cat_eg3_z[slct_indx] = dtk.gio_read(halo_shape_red_step_loc,'evec3z')[indx_slct]
    rotate_host_halo(halo_rot, eg_cat_eg1_x, eg_cat_eg1_y, eg_cat_eg1_z)
    rotate_host_halo(halo_rot, eg_cat_eg2_x, eg_cat_eg2_y, eg_cat_eg2_z)
    rotate_host_halo(halo_rot, eg_cat_eg3_x, eg_cat_eg3_y, eg_cat_eg3_z)
    hgroup['hostHaloEigenValueReduced1'] = eg_cat_eg1
    hgroup['hostHaloEigenValueReduced2'] = eg_cat_eg2
    hgroup['hostHaloEigenValueReduced3'] = eg_cat_eg3
    hgroup['hostHaloEigenVectorReduced1X'] = eg_cat_eg1_x
    hgroup['hostHaloEigenVectorReduced1Y'] = eg_cat_eg1_y
    hgroup['hostHaloEigenVectorReduced1Z'] = eg_cat_eg1_z
    hgroup['hostHaloEigenVectorReduced2X'] = eg_cat_eg2_x
    hgroup['hostHaloEigenVectorReduced2Y'] = eg_cat_eg2_y
    hgroup['hostHaloEigenVectorReduced2Z'] = eg_cat_eg2_z
    hgroup['hostHaloEigenVectorReduced3X'] = eg_cat_eg3_x
    hgroup['hostHaloEigenVectorReduced3Y'] = eg_cat_eg3_y
    hgroup['hostHaloEigenVectorReduced3Z'] = eg_cat_eg3_z
    return

   
def add_native_umachine(output_fname, umachine_native):
    t1 = time.time()
    h_in = h5py.File(umachine_native,'r')
    hgroup = h5py.File(output_fname, 'r+')['galaxyProperties']
    for key in h_in.keys():
           hgroup['UMachineNative/'+key] = h_in[key].value
    print("done addign umachine quantities. time: {:.2f}".format(time.time()-t1))
    return
 

def add_blackhole_quantities(output_fname, redshift, percentile_sfr):
    hgroup = h5py.File(output_fname,'r+')['galaxyProperties']
    print(hgroup.keys())
    bhm = monte_carlo_black_hole_mass(hgroup['spheroidMassStellar'].value)
    eddington_ratio, bh_acc_rate = monte_carlo_bh_acc_rate(redshift, bhm, percentile_sfr)
    hgroup['blackHoleMass'] = bhm
    hgroup['blackHoleAccretionRate'] = bh_acc_rate*1e9
    hgroup['blackHoleEddingtonRatio'] = eddington_ratio


def add_size_quantities(output_fname):
    hgroup = h5py.File(output_fname,'r+')['galaxyProperties']
    mag_r = hgroup['SDSS_filters/magnitude:SDSS_r:rest'].value
    redshift = hgroup['redshift'].value
    arcsec_per_kpc = interp1d(redshift,cosmo.arcsec_per_kpc_proper(redshift).value)
    size_disk = mc_size_vs_luminosity_late_type(mag_r, redshift)
    size_sphere = mc_size_vs_luminosity_early_type(mag_r, redshift)
    f = arcsec_per_kpc(redshift)
    hgroup['morphology/spheroidHalfLightRadius'] =       size_sphere
    hgroup['morphology/spheroidHalfLightRadiusArcsec'] = size_sphere*f
    hgroup['morphology/diskHalfLightRadius'] =       size_disk
    hgroup['morphology/diskHalfLightRadiusArcsec'] = size_disk*f


def add_ellipticity_quantities(output_fname):
    def gal_zoo_dist(x):
        val = np.zeros_like(x)
        a = 2
        slct = x<0.2
        val[slct] = 0
        
        slct = (0.1<=x) & (x<0.7)
        val[slct] = np.tanh((x[slct]-.3)*np.pi*a) - np.tanh((-0.2)*np.pi*a)
        
        slct = (0.7<=x) & (x<1.0)
        val[slct] = np.tanh(-(x[slct]-.95)*np.pi*6.) - np.tanh((-0.2)*np.pi*a) -(np.tanh(-(0.7-0.95)*np.pi*6)-np.tanh(0.4*np.pi*a))
        
        slct = 1.0<=x
        val[slct] = 0
        return val
    hgroup = h5py.File(output_fname, 'r+')['galaxyProperties']
    inclination = hgroup['morphology/inclination'].value
    del hgroup['morphology/inclination']
    size = np.size(inclination)
    pos_angle = np.random.uniform(size=size)*np.pi

    spheroid_axis_ratio = dtk.clipped_gaussian(0.8, 0.2, size, max_val = 1.0, min_val=0.0)
    dist,lim = dtk.make_distribution(-inclination)
    resamp = dtk.resample_distribution(dist,gal_zoo_dist,lim,[0.0,1.0])


    ellip_spheroid = (1.0 - spheroid_axis_ratio)/(1.0 + spheroid_axis_ratio)

    hgroup['morphology/spheroidAxisRatio'] = spheroid_axis_ratio
    hgroup['morphology/spheroidMajorAxisArcsec'] = hgroup['morphology/spheroidHalfLightRadiusArcsec'].value
    hgroup['morphology/spheroidMinorAxisArcsec'] = hgroup['morphology/spheroidHalfLightRadiusArcsec'].value*spheroid_axis_ratio
    hgroup['morphology/spheroidEllipticity'] = ellip_spheroid
    hgroup['morphology/spheroidEllipticity1'] = np.cos(2.0*pos_angle)*ellip_spheroid
    hgroup['morphology/spheroidEllipticity2'] = np.sin(2.0*pos_angle)*ellip_spheroid

    disk_axis_ratio = resamp(-inclination)
    ellip_disk = (1.0 - disk_axis_ratio)/(1.0 + disk_axis_ratio)
    hgroup['morphology/diskAxisRatio'] = disk_axis_ratio
    hgroup['morphology/diskMajorAxisArcsec'] = hgroup['morphology/diskHalfLightRadiusArcsec'].value
    hgroup['morphology/diskMinorAxisArcsec'] = hgroup['morphology/diskHalfLightRadiusArcsec'].value*disk_axis_ratio
    hgroup['morphology/diskEllipticity'] = ellip_disk
    hgroup['morphology/diskEllipticity1'] = np.cos(2.0*pos_angle)*ellip_disk
    hgroup['morphology/diskEllipticity2'] = np.sin(2.0*pos_angle)*ellip_disk

    lum_disk = hgroup['SDSS_filters/diskLuminositiesStellar:SDSS_r:rest'].value
    lum_sphere = hgroup['SDSS_filters/spheroidLuminositiesStellar:SDSS_r:rest'].value
    lum_tot = lum_disk + lum_sphere
    tot_ellip =  (lum_disk*ellip_disk + lum_sphere*ellip_spheroid)/(lum_tot)
    hgroup['morphology/totalEllipticity']  = tot_ellip
    hgroup['morphology/totalAxisRatio']    = (1.0 - tot_ellip)/(1.0 + tot_ellip)
    hgroup['morphology/totalEllipticity1'] = np.cos(2.0*pos_angle)*tot_ellip
    hgroup['morphology/totalEllipticity2'] = np.sin(2.0*pos_angle)*tot_ellip
    
    hgroup['morphology/positionAngle'] = pos_angle*180.0/np.pi


def combine_step_lc_into_one(step_fname_list, out_fname):
    print("combining into one file")
    hfile_out = h5py.File(out_fname,'w')
    hfile_gp_out = hfile_out.create_group('galaxyProperties')
    hfile_steps = []
    hfile_steps_gp = []
    for fname in step_fname_list:
        hfile = h5py.File(fname,'r')
        gp = hfile['galaxyProperties']
        hfile_steps.append(hfile)
        hfile_steps_gp.append(gp)
    keys = get_keys(hfile_steps_gp[0])
    for i,key in enumerate(keys):
        t1 = time.time()
        print("{}/{} {}".format(i,len(keys),key))
        data_list = []
        #units = None
        for h_gp in hfile_steps_gp:
            data_list.append(h_gp[key].value)
            #units = h_gp[key].attrs['units']
        data = np.concatenate(data_list)
        hfile_gp_out[key]=data
        #hfile_gp_out[key].attrs['units']=units
        print("\t time: {:.2f}".format(time.time()-t1))
    hfile_gp_out['galaxyID'] = np.arange(hfile_gp_out['redshift'].size,dtype='i8')
    return 


def add_metadata(gal_ref_fname, out_fname, version_major, version_minor, version_minor_minor):
    """
    Takes the metadata group and copies it over the final output product. 
    Also for each data column, copies the units attribute. 
    """
    add_units(out_fname)
    hfile_gf = h5py.File(gal_ref_fname,'r')
    hfile_out = h5py.File(out_fname,'a')
    # keys_a = get_keys(hfile_gf['galaxyProperties'])
    # keys_b = get_keys(hfile_out['galaxyProperties'])
    # #assert(len(keys_a) == len(keys_b))
    # print(keys_a)
    # print(keys_b)
    # for key in keys_b:
    #     print(key)
    #     if( key not in keys_a):
    #         print("^^^this ========== this^^^")
    # for key in keys_b:
    #     hfile_out['galaxyProperties'][key].attrs['units'] = hfile_gf['galaxyProperties'][key].attrs['units']
    # #copy over metadata
    # del hfile_out['/metaData']
    hfile_out.copy(hfile_gf['metaData'],'metaData')
    del hfile_out['/metaData/catalogCreationDate']
    del hfile_out['/metaData/versionChangeNotes']
    del hfile_out['/metaData/versionMajor'] 
    del hfile_out['/metaData/versionMinor'] 
    del hfile_out['/metaData/versionMinorMinor'] 
    del hfile_out['/metaData/version'] 
    # del hgroup['versionChangeNotes']

    hfile_out['/metaData/versionMajor'] = version_major
    hfile_out['/metaData/versionMinor'] = version_minor
    hfile_out['/metaData/versionMinorMinor'] = version_minor_minor
    hfile_out['/metaData/version'] = "{}.{}.{}".format(version_major, version_minor, version_minor_minor)

    hfile_out['/metaData/catalogCreationDate']=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")


def add_units(out_fname):
    hfile = h5py.File(out_fname,'a')['galaxyProperties']
    #################################
    ###Add units to all fields
    #################################
    mag_list = ['magnitude']; mag_unit = 'AB magnitude'
    arcsec_list= ['Arcsec']; arcsec_unit = 'arcsecond'
    rad_list = []; rad_unit ='radians'
    deg_list = ['ra','dec','positionAngle','inclination']; deg_unit = 'degrees'
    phys_kpc_list = ['Radius']; phys_kpc_unit = 'physical kpc'
    phys_mpc_list = []; phys_mpc_unit = 'physical Mpc'
    reduced_dist_list =['Reduced','EigenVector', 'Eddington'];reduced_dist_unit = 'unitless'
    eigen_val_list = ['EigenValue'];eigen_val_unit = 'comoving Mpc/h'
    comv_mpc_list = ['x','y','z']; comv_mpc_unit = 'comoving Mpc/h'
    vel_list = ['vx','vy','vz','Velocity']; vel_unit = 'km/s'
    timeSFR_list =['TimeWeightedIntegratedSFR']; timeSFR_unit = 'Gyr*Msun'
    sfr_list =['SFR','blackHoleAccretionRate','StarFormationRate']; sfr_unit = 'Msun/Gyr'
    mass_list =['Mass','IntegratedSFR']; mass_unit = 'Msun'
    abundance_list =['Abundance'];abundance_unit = 'Msun'
    luminosity_list =['Luminosities','Luminosity']; luminosity_unit = 'AB luminosity (4.4659e13 W/Hz)'
    unitless_list = ['redshift','shear','magnification','convergence','Ellipticity','Sersic','AxisRatio','dustFactor']; unitless_unit ='unitless'
    id_list = ['Index','Tag','placementType','galaxyID','lightcone_replication','lightcone_rotation']; id_unit = 'id/index'
    angular_list = ['angularMomentum'];angular_unit = 'Msun*km/s*Mpc'
    bool_list =['nodeIsIsolated'];bool_unit = 'boolean'
    spinSpin_list =['spinSpin'];spinSpin_unit ='lambda'
    step_list = ['step'];step_unit = 'simluation step'
    umachine_list = ['UMachineNative'];umachine_unit = 'Unspecified'
    count_list =['NumberSelected'];count_unit = 'count'
    print("assigning units")
    keys = get_keys(hfile)
    print( keys)
    for key in keys:
        print(key)
        #add magnitude units
        if(any(l in key for l in mag_list)):
            hfile[key].attrs['units']=mag_unit
            print("\t mag")
            # umachine list
        elif(any(l in key for l in umachine_list)):
            hfile[key].attrs['units']=umachine_unit
            #add arcsec units
        elif(any(l in key for l in arcsec_list)):
            hfile[key].attrs['units']=arcsec_unit
            print( "\t ",arcsec_unit)
            #add rad units
        elif(any(l in key for l in rad_list)):
            hfile[key].attrs['units']=rad_unit
            print( "\t ",rad_unit)
            #add degree units
        elif(any(l in key for l in deg_list)):
            hfile[key].attrs['units']=deg_unit
            print( '\t',deg_unit)
            #add kpc units
        elif(any(l in key for l in phys_kpc_list)):
            hfile[key].attrs['units']=phys_kpc_unit
            print( "\t ",phys_kpc_unit)
            #add mpc units
        elif(any(l in key for l in phys_mpc_list)):
            hfile[key].attrs['units']=phys_mpc_unit
            print ("\t ",phys_mpc_unit)
            #reduced distances units
        elif(any(l in key for l in reduced_dist_list)):
            hfile[key].attrs['units']=reduced_dist_unit
            print ("\t ",reduced_dist_unit)
            #eigen val units
        elif(any(l in key for l in eigen_val_list)):
            hfile[key].attrs['units']=eigen_val_unit
            print ("\t ",reduced_dist_unit)
            #add comoving mpc units
        elif(any(l == key for l in comv_mpc_list)):
            hfile[key].attrs['units']=comv_mpc_unit
            print ("\t ",comv_mpc_unit)
            #add velocity units
        elif(any(l in key for l in vel_list)):
            hfile[key].attrs['units']=vel_unit
            print ("\t ",vel_unit)
            #add timesfr
        elif(any(l in key for l in timeSFR_list)):
            hfile[key].attrs['units']=timeSFR_unit
            print ("\t ",timeSFR_unit)
            #add sfr
        elif(any(l in key for l in sfr_list)):
            hfile[key].attrs['units']=sfr_unit
            print ("\t ",sfr_unit)
            #add mass
        elif(any(l in key for l in mass_list)):
            hfile[key].attrs['units']=mass_unit
            print ("\t ",mass_unit)
            #add abundance
        elif(any(l in key for l in abundance_list)):
            hfile[key].attrs['units']=abundance_unit
            print ("\t ",abundance_unit)
            #add luminosity units
        elif(any(l in key for l in luminosity_list)):
            hfile[key].attrs['units']=luminosity_unit
            print ("\t ",luminosity_unit)
            #add unit less
        elif(any(l in key for l in unitless_list)):
            hfile[key].attrs['units']=unitless_unit
            print ("\t ",unitless_unit)
            #add mass units
        elif(any(l in key for l in id_list)):
            hfile[key].attrs['units']=id_unit
            print ("\t ",id_unit)
            #angular momentum 
        elif(any(l in key for l in angular_list)):
            hfile[key].attrs['units']=angular_unit
            print ("\t ",angular_unit)
            #boolean
        elif(any(l in key for l in bool_list)):
            hfile[key].attrs['units']=bool_unit
            print ("\t", bool_unit)
            #spinSpin
        elif(any(l in key for l in spinSpin_list)):
            hfile[key].attrs['units']=spinSpin_unit
            # step
        elif(any(l in key for l in step_list)):
            hfile[key].attrs['units']=step_unit
            #counts
        elif(any(l in key for l in count_list)):
            hfile[key].attrs['units']=count_unit
            #Everything should have a unit!
        else:
            print("column", key, "was not assigned a unit :(")
            raise;


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
    


if __name__ == "__main__":
    t00 = time.time()
    param = dtk.Param(sys.argv[1])
    lightcone_fname = param.get_string('lightcone_fname')
    gltcs_fname = param.get_string('gltcs_fname')
    gltcs_metadata_ref = param.get_string('gltcs_metadata_ref')
    gltcs_slope_fname = param.get_string('gltcs_slope_fname')
    sod_fname = param.get_string("sod_fname")
    halo_shape_fname = param.get_string("halo_shape_fname")
    halo_shape_red_fname = param.get_string("halo_shape_red_fname")
    output_fname = param.get_string('output_fname')
    steps = param.get_int_list('steps')
    use_slope = param.get_bool('use_slope')
    short = param.get_bool('short')
    supershort = param.get_bool('supershort')
    substeps = param.get_int('substeps')
    use_substep_redshift = param.get_bool('use_substep_redshift')
    load_mask = param.get_bool("load_mask")
    mask_loc  = param.get_string("mask_loc")
    index_loc = param.get_string("index_loc")
    recolor = param.get_bool('recolor')
    plot = param.get_bool('plot')
    plot_substep = param.get_bool('plot_substep')
    use_dust_factor = param.get_bool('use_dust_factor')
    dust_factors = param.get_float_list('dust_factors')
    ignore_mstar = param.get_bool('ignore_mstar')
    version_major = param.get_int('version_major')
    version_minor = param.get_int('version_minor')
    version_minor_minor = param.get_int('version_minor_minor')
    if load_mask:
        hfile_mask = h5py.File(mask_loc,'r')
    else:
        selection1 = galmatcher.read_selections(yamlfile='galmatcher/yaml/vet_protoDC2.yaml')
        selection2 = galmatcher.read_selections(yamlfile='galmatcher/yaml/colors_protoDC2.yaml')
    stepz = dtk.StepZ(200,0,500)
    output_step_list = []
    step_size = steps.size
    for i in range(0,step_size-1):
        t0 = time.time()
        step = steps[i+1]
        step2 = steps[i]
        print("\n\n=================================")
        print(" STEP: ",step)
        gltcs_step_fname = gltcs_fname.replace("${step}",str(step)) 
        gltcs_slope_step_fname = gltcs_slope_fname.replace("${step}",str(step))
        lightcone_step_fname = lightcone_fname.replace("${step}",str(step))
        output_step_loc = output_fname.replace("${step}",str(step))
        output_step_list.append(output_step_loc)
        sod_step_loc = sod_fname.replace("${step}",str(step))
        halo_shape_step_loc = halo_shape_fname.replace("${step}",str(step))
        halo_shape_red_step_loc = halo_shape_red_fname.replace("${step}",str(step))
        if load_mask:
            mask1 = hfile_mask['{}'.format(step)].value
            mask2 = hfile_mask['{}'.format(step2)].value
        else:
            mask_a = galmatcher.mask_cat(h5py.File(gltcs_step_fname, 'r'), selections=selection1)
            mask_b = galmatcher.mask_cat(h5py.File(gltcs_step_fname, 'r'), selections=selection2)
            mask1 = mask_a & mask_b
            gltcs_step2_fname = gltcs_fname.replace("${step}",str(step2))
            mask_a = galmatcher.mask_cat(h5py.File(gltcs_step2_fname, 'r'), selections=selection1)
            mask_b = galmatcher.mask_cat(h5py.File(gltcs_step2_fname, 'r'), selections=selection2)
            mask2 = mask_a & mask_b
        #The index remap galaxies in step2 to the same order as they were in step1
        index_2to1 = h5py.File(index_loc.replace("${step}",str(step)), 'r')['match_2to1']
        verbose = True
        lc_data = construct_lc_data(lightcone_step_fname, verbose = verbose, recolor=recolor)
        if(use_slope):
            print("using slope", step)
            lc_a = 1.0/(1.0 +lc_data['redshift'])
            lc_a_cc = np.copy(lc_a) # galaxy scale factor for copy columns
            del_lc_a =  np.max(lc_a) - np.min(lc_a)
            step_a = np.min(lc_a)-0.01*del_lc_a #so no galaxy is exactly on the egdge of the bins
            step2_a = np.max(lc_a)+0.01*del_lc_a
            print("min a: {} max a {}".format(step_a,step2_a))
            print("raw min a: {} raw max a: {}".format(np.min(lc_a),np.max(lc_a)))
            abins = np.linspace(step_a, step2_a,substeps+1)
            abins_avg = dtk.bins_avg(abins)
            index = -1*np.ones(lc_data['redshift'].size,dtype='i8')
            match_dust_factors = -1*np.ones(lc_data['redshift'].size,dtype='i4')
            for k in range(0,abins_avg.size):
                print("\t{}/{} substeps".format(k,abins_avg.size))
                slct_lc_abins1 = (abins[k]<=lc_a) 
                slct_lc_abins2 = (lc_a<abins[k+1])
                print("\t\t {} -> {}".format(abins[k],abins[k+1]))
                slct_lc_abin = slct_lc_abins1 & slct_lc_abins2
                print("\t\t num gals: {}".format(np.sum(slct_lc_abin)))
                lc_data_a = dic_select(lc_data, slct_lc_abin)
                if lc_data_a['redshift'].size == 0:
                    print("\t\t\t no galaxies for this redshift bin")
                    continue #nothing to match for this redshift bin
                if use_dust_factor:
                    gal_prop_list = [] 
                    for dust_factor in np.concatenate(([1.0],dust_factors)):
                        print("dust_factor********=",dust_factor)
                        # gal_prop_tmp,_ = construct_gal_prop_redshift_dust(gltcs_step_fname, gltcs_slope_step_fname,
                        #                                                          step_a, abins_avg[k],
                        #                                                          verbose = verbose,
                        #                                                          mask = mask1,
                        #                                                          dust_factor=dust_factor)
                        gal_prop_tmp2 = construct_gal_prop_redshift_dust_raw(
                            gltcs_fname, index_2to1, step, step2, abins_avg[k],
                            mask1, mask2, dust_factor)
                        # plot_gal_prop_dist([gal_prop_tmp, gal_prop_tmp2], ["old method", "new method"])
                        # plt.show()
                        gal_prop_list.append(gal_prop_tmp2)
                    gal_prop_a = cat_dics(gal_prop_list)
                else:
                    gal_prop_a, mask = construct_gal_prop_redshift_dust(gltcs_step_fname, gltcs_slope_step_fname,
                                                                        step_a, abins_avg[k],
                                                                        verbose = verbose,
                                                                        mask = mask1)
                # Find the closest Galacticus galaxy
                index_abin = resample_index(lc_data_a, gal_prop_a, ignore_mstar = ignore_mstar, verbose = verbose)
                if use_dust_factor:
                    # Get the Galacticus galaxy index, the division is to correctly
                    # offset the index for the extra dust gal_prop 
                    print('index_abin: ', np.min(index_abin), np.max(index_abin))
                    print('gal_prop_a: ', np.min(gal_prop_a['index']), np.max(gal_prop_a['index']))
                    a = gal_prop_a['index'][index_abin]
                    index[slct_lc_abin] = a
                    # = index_abin%(index_abin.size//(1+len(dust_factors)))
                    # Record the dust factor for the matched galaxy so that it can be applied 
                    # to other columns in copy_columns()
                    match_dust_factors[slct_lc_abin] = gal_prop_a['dust_factor'][index_abin]
                    if use_substep_redshift:
                        lc_a_cc[slct_lc_abin] = abins_avg[k]
                else:
                    # record the Galacticus galaxy index, 
                    index[slct_lc_abin] = index_abin
                    # all matches have 1.0 dust factor since aren't applying extra dust factors
                    match_dust_factors[slct_lc_abin] = 1.0
                if plot_substep:
                    print("subplots")
                    plot_differences(lc_data_a, gal_prop_a, index_abin);print("subplots")
                    plot_differences_2d(lc_data_a, gal_prop_a, index_abin);print("subplots")
                    plot_side_by_side(lc_data_a, gal_prop_a, index_abin);print("subplots")
                    mag_bins = (-21,-20,-19);print("subplots")
                    plot_mag_r(lc_data_a, gal_prop_a, index_abin);print("subplots")
                    plot_clr_mag(lc_data_a, gal_prop_a, index_abin, mag_bins, 'clr_gr', 'g-r color');print("subplots")
                    #plot_clr_mag(lc_data, gal_prop_a, index_abin, mag_bins, 'clr_ri', 'r-i color')
                    plot_ri_gr_mag(lc_data_a, gal_prop_a, index_abin, mag_bins);print("subplots")
                    plt.show()
                    # h_in_gp = h5py.File(gltcs_step_fname,'r')['galaxyProperties']
                    # h_in_slope_gp = h5py.File(gltcs_slope_step_fname,'r')['galaxyProperties']
                    # del_a = lc_a-step_a
                    # lum_g = get_column_slope_dust("SDSS_filters/totalLuminositiesStellar:SDSS_g:rest:dustAtlas", h_in_gp, h_in_slope_gp, index, del_a, mask1, match_dust_factors)
                    # lum_r = get_column_slope_dust("SDSS_filters/totalLuminositiesStellar:SDSS_r:rest:dustAtlas", h_in_gp, h_in_slope_gp, index, del_a, mask1, match_dust_factors)
                    # lum_i = get_column_slope_dust("SDSS_filters/totalLuminositiesStellar:SDSS_i:rest:dustAtlas", h_in_gp, h_in_slope_gp, index, del_a, mask1, match_dust_factors)
                    # sm = get_column_slope_dust("totalMassStellar", h_in_gp, h_in_slope_gp, index, del_a, mask1, match_dust_factors)
                    # mag_g = -2.5*np.log10(lum_g)
                    # mag_r = -2.5*np.log10(lum_r)
                    # mag_i = -2.5*np.log10(lum_i)
                    # clr_gr = mag_g-mag_r
                    # clr_ri = mag_r-mag_i
                    
            slct_neg = index == -1
            print(match_dust_factors)
            print("not assigned: {}/{}: {:.2f}".format( np.sum(slct_neg), slct_neg.size, np.float(np.sum(slct_neg))/np.float(slct_neg.size)))
            assert(np.sum(slct_neg) == 0)
            # TODO Only allow 
            if use_dust_factor: 
                # copy_columns_slope_dust(gltcs_step_fname, gltcs_slope_step_fname, 
                #                         output_step_loc, index, 
                #                         step_a, lc_a_cc,
                #                         verbose=verbose, mask=mask1, short = short, supershort=supershort,
                #                         step = step, dust_factors = match_dust_factors)
                copy_columns_interpolation_dust_raw(gltcs_fname, output_step_loc, index, 
                                                step, step2, step_a, step2_a, mask1, mask2, 
                                                index_2to1, lc_a_cc, verbose = verbose,
                                                short = short, supershort = supershort,
                                                    dust_factors = match_dust_factors, step = step)
            else:
                copy_columns_slope(gltcs_step_fname, gltcs_slope_step_fname, 
                                   output_step_loc, index, 
                                   step_a, lc_a,
                                   verbose=verbose, mask=mask1, 
                                   short = short, supershort = supershort, step = step)
            
        else:
            print("using no slope", step)
            gal_prop_simple,mask_simple = construct_gal_prop(gltcs_step_fname, verbose = verbose,mask =mask1)
            if use_dust_factor:
                masks = []
                gal_props = []
                masks.append(mask_simple)
                gal_props.append(gal_prop_simple)
                for dust_factor in dust_factors:
                    gp, m = construct_gal_prop_dust_factor(gltcs_step_fname, dust_factor, verbose = verbose,mask =mask1)
                    gal_props.append(gp)
                    masks.append(m)
                gal_prop = cat_dics(gal_props)
                mask = mask_simple#//np.concatenate(masks)
            else:
                gal_prop = gal_prop_simple
                mask = mask_simple
            print("mstar before call: ", ignore_mstar)
            index = resample_index(lc_data, gal_prop, ignore_mstar = ignore_mstar, verbose = verbose)

            if plot:
                #plot_differences(lc_data, gal_prop, index)
                #plot_differences_2d(lc_data, gal_prop, index)
                #plot_differences_2d(lc_data, gal_prop, index, x='m_star')
                mag_bins = (-21,-20,-19)
                #plot_mag_r(lc_data, gal_prop, index)
                plot_m_star(lc_data, gal_prop, index)
                #plot_side_by_side(lc_data, gal_prop, index)
                plot_single_dist(lc_data, gal_prop, index, 'clr_gr', 'g-r',bins=np.linspace(0,1.2,100))
                plot_clr_mag(lc_data, gal_prop, index, mag_bins, 'clr_gr', 'g-r color')
                #plot_clr_mag(lc_data, gal_prop, index, mag_bins, 'clr_ri', 'r-i color')
                # plot_ri_gr_mag(lc_data, gal_prop, index, mag_bins)
                
            if use_dust_factor:
                gal_prop_size = gal_prop['Mag_r'].size//(1+len(dust_factors)) #The data is replicated 1+N where n is the number of dust factors used
                dust_mod_index = index //( gal_prop_size) # which dust factor was matched in the frankenstein
                raw_index = index % (gal_prop_size) #gal_prop index not looking at dust factor
                #Now we get the dust factor for each match
                matched_dust_factor = np.zeros(lc_data['redshift'].size,dtype='f4')
                print(0,"--",len(dust_factors)+1)
                print("dust_mod_index", dust_mod_index)
                print("raw_index", raw_index)
                for ii in range(0,len(dust_factors)+1):
                    slct = dust_mod_index == ii
                    print(ii,'dust_index {}/{}'.format(np.sum(slct),slct.size))
                    if ii == 0:
                        print("ii == 0")
                        matched_dust_factor[slct] = 1.0
                        print(matched_dust_factor[slct])
                    else:
                        print("ii != 0")
                        print(dust_factors)
                        print(dust_factors[ii-1])
                        matched_dust_factor[slct] = dust_factors[ii-1]
                        print(matched_dust_factor[slct])
                print("matched_dust_factor")
                print(matched_dust_factor)
                print(np.min(matched_dust_factor),"->",np.max(matched_dust_factor))
                copy_columns_dust(gltcs_step_fname, output_step_loc, 
                                  raw_index, matched_dust_factor,
                                  verbose = verbose,mask = mask, short = short, step = step)
            else:
                copy_columns(gltcs_step_fname, output_step_loc, index, verbose = verbose,mask = mask, short = short, step = step)
   
        overwrite_columns(lightcone_step_fname, output_step_loc, ignore_mstar = ignore_mstar, verbose = verbose)
        overwrite_host_halo(output_step_loc,sod_step_loc, halo_shape_step_loc, halo_shape_red_step_loc, verbose=verbose)
        add_native_umachine(output_step_loc, lightcone_step_fname)
        add_blackhole_quantities(output_step_loc, np.average(lc_data['redshift']), lc_data['sfr_percentile'])
        add_size_quantities(output_step_loc)
        add_ellipticity_quantities(output_step_loc)
        if plot:
            dummy_mask = np.ones(lc_data['redshift'].size,dtype=bool)
            new_gal_prop,new_mask = construct_gal_prop(output_step_loc, verbose=verbose,mask=dummy_mask)
            index = np.arange(lc_data['redshift'].size)
            plt.figure()
            plt.title("Post match")
            # plt.figure()
            # plt.plot(lc_data['clr_gr'], new_gal_prop['clr_gr'], '.', alpha=0.3)
            # plt.figure()
            # plt.plot(lc_data['Mag_r'], new_gal_prop['Mag_r'], '.', alpha=0.3)
            # plt.figure()
            # plt.plot(lc_data['m_star'], new_gal_prop['m_star'], '.', alpha=0.3)
            # plt.figure()
            # plt.title(" org gal prop vs new gal prop")
            # plt.plot(new_gal_prop['m_star'][index], new_gal_prop['m_star'],'.',alpha=0.3)
            mag_bins = (-21,-20,-19)
            plot_differences(lc_data, new_gal_prop, index)
#            plot_differences_2d(lc_data_a, new_gal_prop, index)
            plot_side_by_side(lc_data, new_gal_prop, index)
            plot_mag_r(lc_data, new_gal_prop, index)
            plot_side_by_side(lc_data, new_gal_prop, index)
            plot_clr_mag(lc_data, new_gal_prop, index, mag_bins, 'clr_gr', 'g-r color')
            #plot_clr_mag(lc_data, new_gal_prop, index, mag_bins, 'clr_ri', 'r-i color')
            plot_ri_gr_mag(lc_data, new_gal_prop, index, mag_bins)
        if plot or plot_substep:
            plt.show()
        print("\n=====\ndone. {}".format(time.time()-t0))
    output_all = output_fname.replace("${step}","all")
    combine_step_lc_into_one(output_step_list, output_all)
    add_metadata(gltcs_metadata_ref, output_all, version_major, version_minor, version_minor_minor)
    print("\n\n========\nALL DONE. Answer correct. \ntime: {:.2f}".format(time.time()-t00))
