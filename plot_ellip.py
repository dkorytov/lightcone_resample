#!/usr/bin/env python2.7

from __future__ import division,print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk 
import h5py
import sys
import time
from numpy.random import normal



if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    output_fname = param.get_string("output_fname")
    output_fname = output_fname.replace("${step}","all")
    hgroup = h5py.File(output_fname,'r')['galaxyProperties']
    minor = hgroup['morphology/spheroidMinorAxisArcsec'].value
    major = hgroup['morphology/spheroidMajorAxisArcsec'].value
    q = hgroup['morphology/spheroidAxisRatio'].value
    e = hgroup['morphology/spheroidEllipticity'].value
    a  = (minor/major)
    b = (1.0-e)/(1.0+e)
    equal =  a==b 
    equal2 = np.isclose(a,b)
    print( equal)
    print(np.sum(equal),equal.size)
    print(np.sum(equal2),equal2.size)
    for i in range(0,10):
        print(equal[i], "  {} == {} \t\t diff: {}".format(a[i],b[i],a[i]-b[i] ))
    print("max: ", np.nanmax(np.abs(a-b)))
    print(np.sum(np.isfinite(a)))
    print(np.sum(np.isfinite(b)))
    # plt.figure()
    # plt.plot(minor,major,'.',alpha=0.3)
    # plt.yscale('log');plt.xscale('log')

    # plt.figure()
    # plt.plot(minor/major,q,'.',alpha=0.3)


    # plt.figure()
    # plt.plot((1.0 - q)/(1.0 +q), e,'.',alpha=0.3)
    

    
    

    # plt.show()
    
