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

def construct_gal_prop(fname,verbose=False,mask = None,mag_r_cut = False):
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
    if verbose:
        print('done loading gal prop. {}'.format(time.time()-t1))
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


def select_by_index(data,index):
    new_data = {}
    for key in data.keys():
        new_data[key] = data[key][index]
    return new_data


def resample_index(lc_data, gal_prop, nnk = 10, verbose = False):
    if verbose:
        t1 = time.time()
        print("Starting kdtree resampling")
    m_star = lc_data['m_star']
    mag_r  = lc_data['Mag_r']
    clr_gr = lc_data['clr_gr']
    clr_ri = lc_data['clr_ri']
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


def copy_columns(input_fname, output_fname, index, verbose = False,mask = None, short = False):
    h_in = h5py.File(input_fname,'r')
    h_out = h5py.File(output_fname,'w')
    h_in_gp = h_in['galaxyProperties']
    h_out_gp = h_out.create_group('galaxyProperties')
    keys = get_keys(h_in_gp)
    for i in range(0,len(keys)):
        key = keys[i]
        if "dust" in key or "LSST" in key or "SED" in key or "other" in key or "Lines" in key:
            if short:
                continue
        print('{}/{}, {} {}'.format(i,len(keys),float(i)/float(len(keys)), key))
        data = h_in_gp[key].value
        if mask is not None:
            data = data[mask]
        h_out_gp[key]=data[index]
        #TODO add units
        #a = h_in_gp[key].attrs['units'].value
        #h_out_gp[key].attrs['units'] = a
    return
    

def copy_columns_slope(input_fname, input_slope_fname, 
                       output_fname, index,  
                       input_redshift, lc_redshift, 
                       verbose = False, mask = None, short = False):
    lc_a = 1.0/(1.0+lc_redshift)
    input_a = 1.0/(1.0 + input_redshift)
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
        if "dust" in key or "LSST" in key or "SED" in key or "other" in key or "Lines" in key:
            if short:
                continue
        print('{}/{}, {} {}'.format(i,len(keys),float(i)/float(len(keys)), key))
        data = h_in_gp[key].value
        slope = h_in_slope_gp[key].value
        if mask is not None:
            data = data[mask]
            slope = slope[mask]
        print(data.dtype)
        if data.dtype == np.float64 or data.dtype == np.float32:
            print("\tfloaty type")
            new_data = data[index] + slope[index]*del_a
        else:
            print("\tinty type")
            new_data = data[index]
        #TODO Does anything need to stored as double?
        slct_finite = np.isfinite(new_data)
        if(new_data.dtype == np.float64 and np.sum(new_data[slct_finite]>max_float) == 0):
            h_out_gp[key]= new_data.astype(np.float32)
        else:
            h_out_gp[key] = new_data
        #TODO add units
        #a = h_in_gp[key].attrs['units'].value
        #h_out_gp[key].attrs['units'] = a
    return


def overwrite_columns(input_fname, output_fname, verbose=False):
    t1 = time.time()
    if verbose:
        print("Overwriting columns.")
        #sdss = Table.read(input_fname,path='data')
    h_in = h5py.File(input_fname,'r')
    #redshift = np.ones(sdss['x'].quantity.size)*0.1
    h_out = h5py.File(output_fname, 'a')
    print(output_fname)
    print(h_out.keys())
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
    h_out_gp['x'][:]=x
    h_out_gp['y'][:]=y
    h_out_gp['z'][:]=z
    h_out_gp['vx'][:]=vx
    h_out_gp['vy'][:]=vy
    h_out_gp['vz'][:]=vz
    h_out_gp['lightcone_rotation'] = h_in['lightcone_rotation'].value
    h_out_gp['lightcone_replication'] = h_in['lightcone_replication'].value
    t3 = time.time()
    if verbose:
        print("\t done overwriting xyz, v_(xyz)",t3-t2)
    #peculiar velocity
    _,z_obs,v_pec,_,_,_,_ = pecZ(x,y,z,vx,vy,vz,redshift)
    h_out_gp['peculiarVelocity'][:] = v_pec
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
            h_out_gp[new_key][:]=adjust_mag -2.5*np.log10(h_out_gp[key].value)
        # Calculating new rest frame magnitudes
        if("totalLuminositiesStellar" in key and  ":rest" in key and ("SDSS" in key or "LSST" in key)):
            new_key = key.replace("totalLuminositiesStellar","magnitude",1)
            print("making: "+new_key+" from "+key)
            h_out_gp[new_key][:]=-2.5*np.log10(h_out_gp[key].value)

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
    h_out_gp['redshift'][:] = z_obs
    h_out_gp['redshiftHubble'][:] = redshift
    #TODO
    #metadata
    #galaxyID
    h_out_gp['galaxyID'][:]=h_in['lightcone_id'].value

    h_out_gp['ra'] = h_in['ra'].value
    h_out_gp['dec'] = h_in['dec'].value
    central = (h_in['host_centric_x'].value ==0) & (h_in['host_centric_y'].value ==0) & (h_in['host_centric_z'].value == 0)
    h_out_gp['isCentral'][:] = central
    h_out_gp['hostHaloTag'][:] = h_in['target_halo_id'].value
    h_out_gp['hostHaloMass'][:] = h_in['target_halo_mass'].value
    #No longer used
    del h_out_gp['placementType']
    tf = time.time()
    if verbose:
        print("\tDone overwrite columns", tf-t1)


def swap(slct, x1, x2):
    print( slct.size)
    print( x1.size)
    print( x2.size)
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


def overwrite_host_halo(output_fname, sod_loc, halo_shape_loc, halo_shape_red_loc, verbose=False):
    hgroup = h5py.File(output_fname,'r+')['galaxyProperties']
    del hgroup['hostHaloSODTag']
    del hgroup['hostIndex']

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
    hgroup['hostHaloSODMass'][:]=sod_mass

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

    srt = np.argsort(eg_cat_htag)
    indx = dtk.search_sorted(eg_cat_htag,htag_real[slct_step],sorter=srt)
    slct_indx = indx != -1
    slct = slct_step
    slct[slct_step]=slct_indx



    
    
def combine_step_lc_into_one(step_fname_list, out_fname):
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
    for key in keys:
        data_list = []
        #units = None
        for h_gp in hfile_steps_gp:
            data_list.append(h_gp[key].value)
            #units = h_gp[key].attrs['units']
        data = np.concatenate(data_list)
        hfile_gp_out[key]=date
        #hfile_gp_out[key].attrs['units']=units
    return 


def add_metadata_units(gal_prop_fname, out_fname):
    """
    Takes the metadata group and copies it over the final output product. 
    Also for each data column, copies the units attribute. 
    """
    hfile_gp = h5py.File(gal_prop_fname,'r')
    hfile_out = h5py.File(out_fname,'a')
    keys_a = get_keys(hfile_gp['galaxyProperties'])
    keys_b = get_keys(hfile_out['galaxyProperties'])
    assert(len(keys_a) == len(keys_b))
    for key in keys_a:
        assert(key in keys_b)
    for key in keys_a:
        hfile_out['galaxyProperties'][key].attrs['units'] = hfile_gp['galaxyProperties'][key].attrs['units']
    #copy over metadata
    hfile_out.copy(hfile_gp['metaData'],'metaData')
    del hfile_out['/metaData/catalogCreationDate']
    hfile_out['/metaData/catalogCreationDate']=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")


def plot_differences(lc_data,gal_prop,index):
    keys = lc_data.keys()
    dist = {}
    dist_all = None
    for key in keys:
        d = lc_data[key]-gal_prop[key][index]
        dist[key] = d
        if(dist_all == None):
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
    plt.plot(dtk.bins_avg(xbins),h,label='all')
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('distance in match')
    plt.ylabel('count')


def plot_differences_2d(lc_data,gal_prop,index,x='Mag_r'):
    keys = lc_data.keys()
    for key in keys:
        if key == x:
            continue
        plt.figure()
        h,xbins,ybins = np.histogram2d(lc_data[x],lc_data[key]-gal_prop[key][index],bins=(100,100))
        plt.pcolor(xbins,ybins,h.T,cmap='PuBu',norm = clr.LogNorm())
        plt.ylabel("diff {} (orginal-new)".format(key))
        plt.xlabel(x)
        plt.grid()
   
 
def plot_side_by_side(lc_data,gal_prop,index,x='Mag_r'):
    keys = lc_data.keys()
    for key in keys:
        if key == x:
            continue
        fig,axs = plt.subplots(1,3,sharey=True,sharex=True,figsize=(15,5))
        h,xbins,ybins = np.histogram2d(lc_data[x],lc_data[key],bins=(100,100))
        axs[0].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
        axs[0].grid()
        axs[0].set_title('Original Catalog')
        axs[0].set_xlabel(x)
        axs[0].set_ylabel(key)

        h,xbins,ybins = np.histogram2d(gal_prop[x][index],gal_prop[key][index],bins=(xbins,ybins))
        axs[1].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
        axs[1].grid()
        axs[1].set_title('Original Sampled from Target')
        axs[1].set_xlabel(x)
        axs[1].set_ylabel(key)

        h,xbins,ybins = np.histogram2d(gal_prop[x],gal_prop[key],bins=(xbins,ybins))
        axs[2].pcolor(xbins,ybins,h.T,cmap='PuBu',norm=clr.LogNorm())
        axs[2].grid()
        axs[2].set_title('Target Catalog ')
        axs[2].set_xlabel(x)
        axs[2].set_ylabel(key)
    return




if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    lightcone_fname = param.get_string('lightcone_fname')
    gltcs_fname = param.get_string('gltcs_fname')
    gltcs_slope_fname = param.get_string('gltcs_slope_fname')
    sod_fname = param.get_string("sod_fname")
    halo_shape_fname = param.get_string("halo_shape_fname")
    halo_shape_red_fname = param.get_string("halo_shape_red_fname")
    output_fname = param.get_string('output_fname')
    steps = param.get_int_list('steps')
    use_slope = param.get_bool('use_slope')
    selection = galmatcher.read_selections()
    stepz = dtk.StepZ(200,0,500)
    for step in steps:
        t0 = time.time()
        print("\n\n=================================")
        print(" STEP: ",step)
        gltcs_step_fname = gltcs_fname.replace("${step}",str(step))
        gltcs_slope_step_fname = gltcs_slope_fname.replace("${step}",str(step))
        lightcone_step_fname = lightcone_fname.replace("${step}",str(step))
        output_step_loc = output_fname.replace("${step}",str(step))
        sod_step_loc = sod_fname.replace("${step}",str(step))
        halo_shape_step_loc = halo_shape_fname.replace("${step}",str(step))
        halo_shape_red_step_loc = halo_shape_red_fname.replace("${step}",str(step))
        mask = galmatcher.mask_cat(h5py.File(gltcs_step_fname, 'r'), selection)
        verbose = True
        lc_data = construct_lc_data(lightcone_step_fname, verbose = verbose)
        gal_prop,mask = construct_gal_prop(gltcs_step_fname, verbose = verbose,mask =mask)
        index = resample_index(lc_data, gal_prop, verbose = verbose)
        # plot_differences(lc_data,gal_prop,index)
        # plot_differences_2d(lc_data,gal_prop,index)
        # plot_differences_2d(lc_data,gal_prop,index,x='m_star')
        # plot_side_by_side(lc_data,gal_prop,index)
        # dtk.save_figs('figs/{}/{}'.format(sys.argv[1],__file__))
        # plt.show()
        if(use_slope):
            print("using slope", step)
            step_redshift = stepz.get_z(step)
            copy_columns_slope(gltcs_step_fname, gltcs_slope_step_fname, 
                               output_step_loc, index, 
                               step_redshift, lc_data['redshift'],
                               verbose=verbose, mask=mask, short = True)
        else:
            print("using no slope", step)
            copy_columns(gltcs_step_fname, output_step_loc, index, verbose = verbose,mask = mask, short = True)
        overwrite_columns(lightcone_step_fname, output_step_loc, verbose = verbose)
        overwrite_host_halo(output_step_loc,sod_step_loc, halo_shape_step_loc, halo_shape_red_step_loc, verbose=verbose)
        
        print("\n=====\ndone. {}".format(time.time()-t0))
