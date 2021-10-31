import os, sys, math, datetime, dateutil

import copy

import xml.dom as dom
from xml.dom import minidom

import numpy as np
from astropy.io import fits as pf
import scipy.interpolate as intp
import scipy.optimize    as opt
from .config import *
from .obslog import *
from .overscan import *


#-----------------------------------------------------------------
# classes and functions for substracting bias
#-----------------------------------------------------------------
def comb_bias(bias_lst):
    return np.median(np.array([pf.getdata(fp) for fp in bias_lst], dtype=np.float32), axis=0)

def biassubtract(logfile=None, configfile=None):

    ''' overscan corrections'''
    if configfile is None:
       configfile = 'config/instrument.conf'
    conf = load_conf(configfile)
    ccd_map = conf['ccd_map']

    readout_lst = get_readout_lst(ccd_map)
    large_rect = get_large_rect(readout_lst)

    # weight list
    w_lst = np.array([1,3,5,3,1])
    scan_p = len(w_lst)
    if logfile is None:
       logfile = find_logfile()
    log = read_log(logfile)
    print('Obs Log File =',logfile)
    ## list bias file
    bias_lst = []
    for item in log.item_list:
        if item.object.lower() == 'bias':
           filename = item.get_filename(log.filename_composition)
           bias_lst.append(filename)
    master_bias = comb_bias(bias_lst)
    bias = large_rect.cut_data(master_bias)
    for item in log.item_list:
        if (not item.skip) and (item.object.lower() !='bias'):
            filename = item.get_filename(log.filename_composition)

            print(filename,)

            data,head = pf.getdata(filename, header=True)
            data = large_rect.cut_data(data)
            frame = data - bias

            #pf.writeto('bias.fits',bias)
            ovs_filename = '.'.join(filename.split('.')[0:-1])+'_ovr.fits'
            pf.writeto(ovs_filename,frame,head, overwrite=True)
            print('-> %s'%(ovs_filename))
