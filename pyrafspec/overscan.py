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


#-----------------------------------------------------------------
# classes and functions for overscan
#-----------------------------------------------------------------

class Rectangle(object):
    def __init__(self, zero, shape):
        self.zero  = zero
        self.shape = shape
        self.area  = self.shape[0]*self.shape[1]
        self.x0 = self.zero[0]
        self.x1 = self.zero[0] + self.shape[0]
        self.y0 = self.zero[1]
        self.y1 = self.zero[1] + self.shape[1]

    def cut_data(self,data):
        return data[self.y0:self.y1, self.x0:self.x1]

class ReadOut(object):
    def __init__(self, data_rect, overscan_rect, align_axis):
        self.data_rect     = data_rect
        self.overscan_rect = overscan_rect
        self.align_axis    = align_axis.strip().lower()

def get_tuple_by_node(node):
    text_lst = []
    for i in node.childNodes:
        if i.nodeType == i.TEXT_NODE:
             text_lst.append(i.data)

    text = ''.join(text_lst).strip()
    g = text.split(',')
    return int(g[0]),int(g[1])

def get_readout_lst(filename):
    xmldoc = minidom.parse(filename)
    image_node = xmldoc.childNodes[0]
    readout_lst = []
    for readout_node in image_node.childNodes:
        if readout_node.nodeName == 'readout':

            for data_node in readout_node.childNodes:


                if data_node.nodeName=='data':
                    for shape_node in data_node.childNodes:
                        if shape_node.nodeName=='zero':
                            zero = get_tuple_by_node(shape_node)
                        elif shape_node.nodeName=='shape':
                            shape = get_tuple_by_node(shape_node)
                    data_rect = Rectangle(zero,shape)


                elif data_node.nodeName=='overscan':

                    align_axis = data_node.getAttribute('align_axis')

                    for shape_node in data_node.childNodes:
                        if shape_node.nodeName=='zero':
                            zero = get_tuple_by_node(shape_node)
                        elif shape_node.nodeName=='shape':
                            shape = get_tuple_by_node(shape_node)
                    overscan_rect = Rectangle(zero,shape)

            readout = ReadOut(data_rect,overscan_rect,align_axis)
            readout_lst.append(readout)
    return readout_lst

def get_large_rect(readout_lst):
    srect = readout_lst[0].data_rect
    xedge = [srect.x0, srect.x1]
    yedge = [srect.y0, srect.y1]

    for readout in readout_lst[1:]:

        # rect of science
        srect = readout.data_rect

        if srect.x0 < xedge[0]: xedge[0] = srect.x0
        if srect.x1 > xedge[1]: xedge[1] = srect.x1
        if srect.y0 < yedge[0]: yedge[0] = srect.y0
        if srect.y1 > yedge[1]: yedge[1] = srect.y1

    zero = (xedge[0],yedge[0])
    shape = (xedge[1]-xedge[0], yedge[1]-yedge[0])

    return Rectangle(zero,shape)



def get_rainbow_color(wv):
    '''
    see http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    if 3800 <= wv <= 4400:
        r = (4400.-wv)/(4400.-3800.)
        g = 0
        b = 1
    elif 4400 <= wv <= 4900:
        r = 0
        g = (wv-4400.)/(4900.-4400.)
        b = 1
    elif 4900 <= wv <= 5100:
        r = 0
        g = 1
        b = (5100.-wv)/(5100.-4900.)
    elif 5100 <= wv <= 5800:
        r = (wv-5100.)/(5800.-5100.)
        g = 1
        b = 0
    elif 5800 <= wv <= 6450:
        r = 1
        g = (6450.-wv)/(6450.-5800.)
        b = 0
    elif 6450 <= wv <= 7800:
        r = 1
        g = 0
        b = 0
    elif wv > 7800:
        r = 1
        g = 0
        b = 0
    if wv < 3800:
        r = (4400.-wv)/(4400.-3800.)
        if r < 0:
            r = 0
        g = 0
        b = 1

    if wv > 7000:
        s = 0.3 + 0.7*(7800.-wv)/(7800.-7000.)
    elif wv < 4200:
        s = 0.3 + 0.7*(wv-3800.)/(4200.-3800.)
    else:
        s = 1.0

    if s < 0:
        s = 0.

    gamma = 0.8
    red   = int((s*r)**gamma*255)
    green = int((s*g)**gamma*255)
    blue  = int((s*b)**gamma*255)
    return (red,green,blue)


def overscan(logfile=None, configfile=None, rot90=False):

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

    for item in log.item_list:
        if not item.skip:
            filename = item.get_filename(log.filename_composition)

            print(filename,)

            data,head = pf.getdata(filename, header=True)
            bias = np.float32(np.zeros_like(data))

            for readout in readout_lst:
                # rect of overscan region
                orect = readout.overscan_rect
                # rect of science region
                srect = readout.data_rect
                overscan = data[orect.y0:orect.y1,orect.x0:orect.x1]

                if readout.align_axis=='x':
                    over_mean = overscan.mean(axis=1)
                    if srect.shape[1]!=orect.shape[1]:
                        print('Error: Shape does not match')
                        exit()
                    b0,b1 = orect.y0, orect.y1
                elif readout.align_axis=='y':
                    over_mean = overscan.mean(axis=0)
                    if srect.shape[0]!=orect.shape[0]:
                        print('Error: Shape does not match')
                        exit()
                    b0,b1 = orect.x0, orect.x1
                else:
                    print('Align axis error')
                    exit()

                # weighted over mean
                b = np.float32(np.zeros_like(over_mean))
                for i in np.arange(b0,b1):
                    i1 = i-np.int((scan_p-1)//2)
                    i2 = i1 + scan_p
                    if i1 < b0:
                        print(f'i1={i1}, b0={b0}')
                        w = w_lst[b0-i1:]
                        i1 = b0
                    elif i2 > b1:
                        w = w_lst[:b1-i2]
                        i2 = b1
                    else:
                        w = w_lst

                    b[i-b0] = (over_mean[i1-b0:i2-b0]*w).sum()/w.sum()

                if readout.align_axis=='x':
                    for j in np.arange(srect.x0,srect.x1):
                        bias[srect.y0:srect.y1,j] = b
                elif readout.align_axis=='y':
                    for k in np.arange(srect.y0,srect.y1):
                        bias[k,srect.x0:srect.x1] = b
                else:
                    print('Error align axis')
                    exit()

            data = large_rect.cut_data(data)
            bias = large_rect.cut_data(bias)
            frame = data - bias
            if rot90:
               frame = np.rot90(frame)
            #pf.writeto('bias.fits',bias)
            ovs_filename = '.'.join(filename.split('.')[0:-1])+'_ovr.fits'
            pf.writeto(ovs_filename,frame,head, overwrite=True)
            print('-> %s'%(ovs_filename))
