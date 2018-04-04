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
        #TODO check
        mask = (mag_r < -10) & mask
    gal_prop['m_star'] = m_star[mask]
    gal_prop['Mag_r']  = mag_r[mask]
    gal_prop['clr_gr'] = mag_g[mask]-mag_r[mask]
    gal_prop['clr_ri'] = mag_r[mask]-mag_i[mask]
    # gal_prop['Mag_g']  = mag_g[mask]
    # gal_prop['Mag_i']  = mag_i[mask]
    print(mag_g)
    print(mag_r)
    print(mag_i)
    print(gal_prop['clr_gr'])
    print(gal_prop['clr_ri'])
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
        #TODO check
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
        #TODO check
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


def construct_gal_prop_redshift_dust(fname,slope_fname,snap_a,target_a,verbose=False,mask=None,mag_r_cut = False,index = None, dust_factor = 1.0 ):
    t1 = time.time()
    del_a = target_a - snap_a
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
        mask = np.ones(mag_r.size,dtype=bool)
    if mag_r_cut:
        #TODO check
        mask = (mag_r < -10) & mask
    if index is None:
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
    m_star = np.log10(m_star[mask][index] + m_star_slp[mask][index]*del_a)
    gal_prop['m_star'] = m_star
    gal_prop['Mag_r'] = mag_r
    gal_prop['clr_gr'] = mag_g - mag_r
    gal_prop['clr_ri'] = mag_r - mag_i
    gal_prop['dust_factor'] = np.ones(m_star.size,dtype='f4')*dust_factor
    if verbose:
        print('done loading slope gal prop. {}'.format(time.time()-t1))
    return gal_prop,mask


def construct_lc_data(fname,verbose = False):
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
    m_star = lc_data['m_star']
    mag_r  = lc_data['Mag_r']
    clr_gr = lc_data['clr_gr']
    clr_ri = lc_data['clr_ri']
    if not ignore_mstar:
        lc_mat = np.stack((m_star,mag_r,clr_gr,clr_ri),axis=1)
        gal_mat = np.stack((gal_prop['m_star'],
                            gal_prop['Mag_r'],
                            gal_prop['clr_gr'],
                            gal_prop['clr_ri']),axis=1)
    else:
        lc_mat = np.stack((mag_r,clr_gr,clr_ri),axis=1)
        gal_mat = np.stack((gal_prop['Mag_r'],
                            gal_prop['clr_gr'],
                            gal_prop['clr_ri']),axis=1)

    if verbose:
        t2 = time.time()
        print('\tdone formating data. {}'.format(t2-t1))
    ckdtree = cKDTree(gal_mat, balanced_tree = False, compact_nodes = False)
    if verbose:
        t3 = time.time()
        print('\tdone making tree. {}'.format(t3-t2))
    dist, index = ckdtree.query(lc_mat, nnk,n_jobs=16)
    if verbose:
        t4= time.time()
        print('\tdone querying. {}'.format(t4-t3))
    if nnk > 1:
        rand = np.random.randint(nnk,size=dist.shape[0])
        aa = np.arange(dist.shape[0])
        #dist = dist[aa,rand]
        index = index[aa,rand]
    return index
                                

def get_keys(hgroup):
    keys = []
    def _collect_keys(name, obj):
        if isinstance(obj, h5py.Dataset): 
            keys.append(name)
    hgroup.visititems(_collect_keys)
    return keys


copy_avoids = ('x','y','z','vx','vy','vz', 'peculiarVelocity','galaxyID','redshift','redshiftHubble','placementType','isCentral','hostIndex')
copy_avoids_ptrn = ('hostHalo','magnitude')
no_slope_var = ('x','y','z','vx','vy','vz', 'peculiarVelocity','galaxyID','redshift','redshiftHubble','inclination','positionAngle')
no_slope_ptrn  =('morphology','hostHalo','infall')

def copy_columns(input_fname, output_fname, index, verbose = False,mask = None, short = False, step = -1):
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
        if "SDSS" not in key and "total" not in key and "totalMassStellar" != key:
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
        #TODO add units
        #a = h_in_gp[key].attrs['units'].value
        #h_out_gp[key].attrs['units'] = a
    return


def copy_columns_dust(input_fname, output_fname, raw_index, dust_factor, verbose = False,mask = None, short = False, step = -1, dust_factors = 1.0):
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
        if "SDSS" not in key and "total" not in key and ":rest" not in key and "totalMassStellar" != key:
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
        #TODO add units
        #a = h_in_gp[key].attrs['units'].value
        #h_out_gp[key].attrs['units'] = a
    return
    

def copy_columns_slope(input_fname, input_slope_fname, 
                       output_fname, index,  
                       input_a, lc_a, 
                       verbose = False, mask = None, short = False, step = -1):
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
        if "SDSS" not in key and "total" not in key and ":rest" not in key and "totalMassStellar" != key:
            if supershort:
                continue
        if any([ ca == key for ca in copy_avoids]) or any([ cap in key for cap in copy_avoids_ptrn ]):
            print("{} isn't copied".format(key))
            continue
   
        if verbose:
            print('{}/{},{} {} {} {}'.format(i,len(keys),step,float(i)/float(len(keys)), key,output_fname))
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


def copy_columns_slope_dust(input_fname, input_slope_fname, 
                       output_fname, index,  
                       input_a, lc_a, 
                       verbose = False, mask = None, short = False, step = -1, dust_factors = 1.0):
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
        if "LSST" in key or "SED" in key or "other" in key or "Lines" in key or "morphology" in key:
            if short:
                continue
        if "SDSS" not in key and "total" not in key and ":rest" not in key and "totalMassStellar" != key:
            if supershort:
                continue
        if any([ ca == key for ca in copy_avoids]) or any([ cap in key for cap in copy_avoids_ptrn ]):
            print("{} isn't copied".format(key))
            continue
   
        if verbose:
            print('{}/{} [{}] {}'.format(i,len(keys),step, key))
        if ":dustAtlas" in key:
            print("\thas dust")
            key_nd = key.replace(":dustAtlas","")
            data = h_in_gp[key_nd].value
            slope = h_in_slope_gp[key_nd].value
            # multiplicitive effect of dust
            data_dust = h_in_gp[key].value
            dust_effect = data_dust/data 
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
                new_data = (data[index] + slope[index]*del_a)*dust_effect[index]*dust_factors
            else:
                print("\t\t does not have dustAtlas")
                new_data = data[index] + slope[index]*del_a
        else:
            new_data = data[index]
        slct_finite = np.isfinite(new_data)
        #If the data is a double, record it as a float to save on disk space
        if(new_data.dtype == np.float64 and np.sum(new_data[slct_finite]>max_float) == 0):
            h_out_gp[key]= new_data.astype(np.float32)
        else:
            h_out_gp[key] = new_data
    return


def overwrite_columns(input_fname, output_fname, verbose=False):
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
    #TODO
    #metadata
    #galaxyID
    h_out_gp['galaxyID']=h_in['lightcone_id'].value

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
    print(keys)
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
    return 


def add_metadata(gal_ref_fname, out_fname):
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
    hfile_out.copy(hfile_gf['metaData'],'metaData')
    del hfile_out['/metaData/catalogCreationDate']
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
    reduced_dist_list =['Reduced','EigenVector'];reduced_dist_unit = 'unitless'
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
    print("assigning units")
    keys = get_keys(hfile)
    print( keys)
    for key in keys:
        print(key)
        #add magnitude units
        if(any(l in key for l in mag_list)):
            hfile[key].attrs['units']=mag_unit
            print("\t mag")
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
            #step
        elif(any(l in key for l in step_list)):
            hfile[key].attrs['units']=step_unit
            #Everything should have a unit!
        else:
            print("column", key, "was not assigned a unit :(")
            raise;


def plot_differences(lc_data, gal_prop, index):
    keys = gal_prop.keys()
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
        h,xbins = np.histogram(dist[key],bins=100)
        plt.plot(dtk.bins_avg(xbins),h,label=key)
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('original value - matched value')
    plt.ylabel('count')
    plt.figure()
    h,xbins = np.histogram(dist_all,bins=100)
    plt.plot(dtk.bins_avg(xbins),h,label='all',lw=2.0)
    for key in keys:
        h,xbins = np.histogram(dist[key],bins=xbins)
        plt.plot(dtk.bins_avg(xbins),h,label=key)
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('distance in match')
    plt.ylabel('count')
    return


def plot_differences_2d(lc_data, gal_prop,index, x='Mag_r'):
    keys = gal_prop.keys()
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
    return


def plot_mag_r(lc_data,gal_prop,index):
    plt.figure()
    max_r = np.max((np.max(lc_data['Mag_r']),np.max(gal_prop['Mag_r'][index])))
    min_r = np.min((np.min(lc_data['Mag_r']),np.min(gal_prop['Mag_r'][index])))
    bins = np.linspace(min_r,max_r,100)
    h_lc,_ = np.histogram(lc_data['Mag_r'],bins=bins)
    h_mg,_ = np.histogram(gal_prop['Mag_r'][index],bins=bins)
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
    axs[0].set_ylabel(key)
    h,xbins,ybins = np.histogram2d(gal_prop['redshift'][index],gal_prop[clr][index],bins=(xbins,ybins))
    axs[1].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    axs[1].grid()
    axs[1].set_title('Matched Galacticus')
    axs[1].set_xlabel('redshift')
    axs[1].set_ylabel(key)
    h,xbins,ybins = np.histogram2d(gal_prop['redshift'],gal_prop[clr],bins=(xbins,ybins))
    axs[2].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
    axs[2].grid()
    axs[2].set_title('Galacticus ')
    axs[2].set_xlabel('redshift')
    axs[2].set_ylabel(key)
    

if __name__ == "__main__":
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
    plot = param.get_bool('plot')
    use_dust_factor = param.get_bool('use_dust_factor')
    dust_factors = param.get_float_list('dust_factors')
    ignore_mstar = param.get_bool('ignore_mstar')
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
        mask1 = galmatcher.mask_cat(h5py.File(gltcs_step_fname, 'r'), selections=selection1)
        mask2 = galmatcher.mask_cat(h5py.File(gltcs_step_fname, 'r'), selections=selection2)
        mask = mask1 & mask2
        verbose = True
        lc_data = construct_lc_data(lightcone_step_fname, verbose = verbose)
        if(use_slope):
            print("using slope", step)
            lc_a = 1.0/(1.0 +lc_data['redshift'])
            step_a = np.min(lc_a)
            step2_a = np.max(lc_a)
            abins = np.linspace(step_a, step2_a,substeps+1)
            abins_avg = dtk.bins_avg(abins)
            index = -1*np.ones(lc_data['redshift'].size,dtype='i4')
            match_dust_factors = -1*np.ones(lc_data['redshift'].size,dtype='i4')
            for k in range(0,abins_avg.size):
                print("\t{}/{} substeps".format(k,abins_avg.size))
                h,xbins = np.histogram(lc_a,bins=100)
                slct_lc_abins1 = (abins[k]<lc_a) 
                slct_lc_abins2 = (lc_a<abins[k+1])
                slct_lc_abin = slct_lc_abins1 & slct_lc_abins2
                lc_data_a = dic_select(lc_data, slct_lc_abin)
                if lc_data_a['redshift'].size == 0:
                    print("\t no galaxies for this redshift bin")
                    continue #nothing to match for this redshift bin
                if use_dust_factor:
                    gal_prop_list = [] 
                    for dust_factor in np.concatenate(([1.0],dust_factors)):
                        print("dust_factor********=",dust_factor)
                        gal_prop_tmp,_ = construct_gal_prop_redshift_dust(gltcs_step_fname, gltcs_slope_step_fname,
                                                                                 step_a, abins_avg[k],
                                                                                 verbose = verbose,
                                                                                 mask = mask,
                                                                                 dust_factor=dust_factor)
                        gal_prop_list.append(gal_prop_tmp)
                    gal_prop_a = cat_dics(gal_prop_list)
                else:
                    gal_prop_a, mask = construct_gal_prop_redshift_dust(gltcs_step_fname, gltcs_slope_step_fname,
                                                                        step_a, abins_avg[k],
                                                                        verbose = verbose,
                                                                        mask = mask)
                # Find the closest Galacticus galaxy
                index_abin = resample_index(lc_data_a, gal_prop_a, verbose = verbose)
                if use_dust_factor:
                    # Get the Galacticus galaxy index, the division is to correctly
                    # offset the index for the extra dust gal_prop 
                    index[slct_lc_abin] = index_abin//(1+len(dust_factors))
                    # Record the dust factor for the matched galaxy so that it can be applied 
                    # to other columns in copy_columns()
                    match_dust_factors[slct_lc_abin] = gal_prop_a['dust_factor'][index_abin]
                else:
                    # record the Galacticus galaxy index, 
                    index[slct_lc_abin] = index_abin
                    # all matches have 1.0 dust factor since aren't applying extra dust factors
                    match_dust_factors[slct_lc_abin] = 1.0
                if plot:
                    plot_differences(lc_data_a, gal_prop_a, index_abin)
                    plot_differences_2d(lc_data_a, gal_prop_a, index_abin)
                    plot_side_by_side(lc_data_a, gal_prop_a, index_abin)
                    mag_bins = (-21,-20,-19)
                    plot_mag_r(lc_data, gal_prop_a, index_abin)
                    plot_clr_mag(lc_data, gal_prop_a, index_abin, mag_bins, 'clr_gr', 'g-r color')
                    plot_clr_mag(lc_data, gal_prop_a, index_abin, mag_bins, 'clr_ri', 'r-i color')
                    plot_ri_gr_mag(lc_data, gal_prop_a, index_abin, mag_bins)
                    plt.show()
            slct_neg = index == -1
            print(match_dust_factors)
            print("not assigned: {}/{}: {:.2f}".format( np.sum(slct_neg), slct_neg.size, np.float(np.sum(slct_neg))/np.float(slct_neg.size)))
            if use_dust_factor: 
                copy_columns_slope_dust(gltcs_step_fname, gltcs_slope_step_fname, 
                                        output_step_loc, index, 
                                        step_a, lc_a,
                                        verbose=verbose, mask=mask, short = short, step = step,
                                        dust_factors = match_dust_factors)
            else:
                copy_columns_slope(gltcs_step_fname, gltcs_slope_step_fname, 
                                   output_step_loc, index, 
                                   step_a, lc_a,
                                   verbose=verbose, mask=mask, short = short, step = step)
            
        else:
            print("using no slope", step)
            gal_prop_simple,mask_simple = construct_gal_prop(gltcs_step_fname, verbose = verbose,mask =mask)
            if use_dust_factor:
                masks = []
                gal_props = []
                masks.append(mask_simple)
                gal_props.append(gal_prop_simple)
                for dust_factor in dust_factors:
                    gp, m = construct_gal_prop_dust_factor(gltcs_step_fname, dust_factor, verbose = verbose,mask =mask)
                    gal_props.append(gp)
                    masks.append(m)
                gal_prop = cat_dics(gal_props)
                mask = mask_simple#//np.concatenate(masks)
            else:
                gal_prop = gal_prop_simple
                mask = mask_simple
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
                #plot_clr_mag(lc_data, gal_prop, index, mag_bins, 'clr_gr', 'g-r color')
                # plot_clr_mag(lc_data, gal_prop, index, mag_bins, 'clr_ri', 'r-i color')
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
   
        overwrite_columns(lightcone_step_fname, output_step_loc, verbose = verbose)
        overwrite_host_halo(output_step_loc,sod_step_loc, halo_shape_step_loc, halo_shape_red_step_loc, verbose=verbose)
        if plot:
            dummy_mask = np.ones(lc_data['redshift'].size,dtype=bool)
            new_gal_prop,new_mask = construct_gal_prop(output_step_loc, verbose=verbose,mask=dummy_mask)
            index = np.arange(lc_data['redshift'].size)
            # plt.figure()
            # plt.plot(lc_data['clr_gr'], new_gal_prop['clr_gr'], '.', alpha=0.3)
            # plt.figure()
            # plt.plot(lc_data['Mag_r'], new_gal_prop['Mag_r'], '.', alpha=0.3)
            # plt.figure()
            # plt.plot(lc_data['m_star'], new_gal_prop['m_star'], '.', alpha=0.3)
            
            plt.figure()
            plt.title(" org gal prop vs new gal prop")
            plt.plot(new_gal_prop['m_star'][index], new_gal_prop['m_star'],'.',alpha=0.3)
            plot_mag_r(lc_data, new_gal_prop, index)
            plot_side_by_side(lc_data, new_gal_prop, index)
            plot_clr_mag(lc_data, new_gal_prop, index, mag_bins, 'clr_gr', 'g-r color')
            plot_clr_mag(lc_data, new_gal_prop, index, mag_bins, 'clr_ri', 'r-i color')
            plot_ri_gr_mag(lc_data, new_gal_prop, index, mag_bins)

            plt.show()
        print("\n=====\ndone. {}".format(time.time()-t0))
    output_all = output_fname.replace("${step}","all")
    combine_step_lc_into_one(output_step_list, output_all)
    add_metadata(gltcs_metadata_ref, output_all)
    
