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


#-----------------------------------------------------------------
# classes and functions for spectrum
#-----------------------------------------------------------------

class Spec(object):
    def __init__(self):
        self.orders = []
        self.order  = {}

    def add_order(self, order, specdata):
        if order not in self.orders:
            self.order[order] = copy.deepcopy(specdata)
            self.orders.append(order)
        else:
            print('Error: Order',order,'already in Spec')
            exit()
        self.orders = self.get_sorted_orders()

    def get_sorted_orders(self):
        """ Sort the order by wv. blue -> red"""
        blue_red_wv_lst = np.array([])
        blue_red_orders = []

        for order in self.order.keys():
            wv = self[order].wv_from
            # find index
            i = np.searchsorted(blue_red_wv_lst, wv)
            # insert wv and order number
            blue_red_wv_lst = np.insert(blue_red_wv_lst,i,wv)
            blue_red_orders.insert(i,order)
        return blue_red_orders

    def __getitem__(self,key):
        '''spec[ord], return specdata instance'''
        order = int(key)
        if order not in self.orders:
            raise IndexError
        return self.order[order]

    def find_best_order(self, wv):

        in_ords = []
        pos_ords = []
        for ord in self.orders:
            if self[ord].wv_from < wv < self[ord].wv_to:
                in_ords.append(ord)
                pos = abs(wv - (self[ord].wv_from+self[ord].wv_to)/2.0)/(self[ord].wv_to-self[ord].wv_from)
                pos_ords.append(pos)

        if len(in_ords)==0:
            return None

        i = np.array(pos_ords).argmin()
        return in_ords[i]


class SpecData(object):
    def __init__(self,wv,flux):
        self.wv      = np.array(wv)
        self.flux    = np.array(flux)
        if self.wv[0] > self.wv[-1]:
            self.wv   = self.wv[::-1]
            self.flux = self.flux[::-1]
        self.wv_from = self.wv[0]
        self.wv_to   = self.wv[-1]
        self.pts = self.wv.size
        self.wv_range = self.wv_to - self.wv_from

    def save_txt(self,filename):
        file = open(filename,'w')
        for i in np.arange(self.pts):
            file.write('%10.4f %10.4f%s'%(self.wv[i],self.flux[i],os.linesep))
        file.close()

    def save_fits(self,filename,samp=1.0):
        '''save fits file
        '''
        pts = int(self.pts*samp)
        dwv = self.wv_range/(pts-1)
        wvnew = np.arange(pts)*dwv + self.wv_from
        tck = intp.splrep(self.wv, self.flux, s=0)
        fluxnew = intp.splev(wvnew,tck,der=0)

        head = pf.Header()
        head['SIMPLE']=(True)
        #head.update('SIMPLE',True)

        #if fluxnew.dtype == 'float64':
        #    head.update('BITPIX', -64)
        #elif fluxnew.dtype == 'float32':
        #    head.update('BITPIX', -32)
        #else:
        #    print 'head Error'
        wvnew   = np.float32(wvnew)
        fluxnew = np.float32(fluxnew)
        head['BITPIX']=(-32)
        #head.update('BITPIX',-32)

        wv_offset = wvnew[0]

        head['NAXIS']=(1)
        head['NAXIS1']=(pts)
        head['CTYPE1']=('WAVELENGTH')
        head['CRPIX1']=(1.0)
        head['CRVAL1']=(np.float32(wv_offset))
        head['CDELT1']=(np.float32(dwv))
        #head['NAXIS']=(1)
        #head.update('NAXIS1',pts)
        #head.update('CTYPE1','WAVELENGTH')
        #head.update('CRPIX1',1.0)
        #head.update('CRVAL1',np.float32(wv_offset))
        #head.update('CDELT1',np.float32(dwv))

        if os.path.exists(filename):
            os.remove(filename)
        pf.writeto(filename,fluxnew,head)

    def measure_snr(self, argmax=20,gain=1.0):
        """Measure SNR from flux. Only works when flux is photon count
        call signature::
            specdata.measureSNR(argmax=50)
        """
        count = np.sort(self.flux)[self.pts-argmax]
        return int(math.sqrt(count*gain))

def load_multispec(filename):
    wv_group, flux_group = specwcs.calib(filename)

    spec = Spec()
    for ord in wv_group.keys():
        specdata = SpecData(wv   = wv_group[ord],
                            flux = flux_group[ord])
        spec.add_order(ord,specdata)

    return spec


class specwcs(object):

    @staticmethod
    def calib(filename):

        data,head = pf.getdata(filename,header=True)
        multispec_items = MultiSpecItem.get_wat2(head)

        wv_group   = {}
        flux_group = {}

        for N in multispec_items.keys():
            ord = multispec_items[N].beam
            wv = multispec_items[N].get_wv()
            if multispec_items[N].dtype==0:
                flux = data[N-1,:]
            elif multispec_items[N].dtype==2:
                pmin = multispec_items[N].pmin[0]
                pmax = multispec_items[N].pmax[0]
                flux = data[N-1,int(pmin)-1:int(pmax)-int(pmin)+1]

            wv_group[ord]   = wv
            flux_group[ord] = flux

        return wv_group,flux_group

class MultiSpecItem(object):
    """Object for MultiSpec Format"""

    def __init__(self,string):
        self.string = string
        self.parse_string()

    def parse_string(self):
        # fix "xxxE-40.xxx problem"

        g = self.string.split(" ")
        for i in np.arange(len(g)):
            if g[i].count('.')==2:
                ## find first dot
                #pos1 = e.find(".")
                ## find second dot
                #pos2 = e[pos1+1:].find(".") + pos1 + 1
                ## replace the second dot
                #g[i] = g[i].replace('0.',' 0.')

                t = g[i].split('.')
                t[2] = t[1][-1]+'.'+t[2]
                t[1] = t[1][:-1]
                g[i] = t[0]+'.'+t[1]+' '+t[2]

        string = ' '.join(g)

        g = string.split()
        self.ap    =   int(g[0])
        self.beam  =   int(g[1])
        self.dtype =   int(g[2])
        self.w1    = float(g[3])
        self.dw    = float(g[4])
        self.nw    =   int(g[5])
        self.z     = float(g[6])
        self.aplow = float(g[7])
        self.aphigh= float(g[8])

        if self.dtype == 0:
            # linear coordinate
            pass

        elif self.dtype == 2:
            # nonlinear coordiante
            func = 0
            self.wt    = {}
            self.wt0   = {}
            self.ftype = {}
            self.order = {}
            self.pmin  = {}
            self.pmax  = {}
            self.c     = {}
            pos       = 0
            while(pos+9<len(g)):
                self.wt[func]    = float(g[ 9+pos])
                self.wt0[func]   = float(g[10+pos])
                self.ftype[func] = round(float(g[11+pos]))
                if self.ftype[func] in [1,2]:
                    # Chebyshev (ftype=1) or Legendre (ftype=2) Polynomial
                    self.order[func] = int(g[12+pos])
                    self.pmin[func]  = float(g[13+pos])
                    self.pmax[func]  = float(g[14+pos])
                    self.c[func]     = []
                    for j in np.arange(self.order[func]):
                        self.c[func].append(float(g[15+pos+j]))
                    pos += 3+3+self.order[func]
                func += 1

    @staticmethod
    def get_wat2(head):
        ''' get WAT2 string, return a turple, consist of
        MultiSpecItem instances '''
        string = ""
        for item in head.items():
            if item[0][0:4]=="WAT2":
                string += item[1].ljust(68)
        ap_num = head['NAXIS2']
        g = string.split('spec')
        pos_list = []
        for ap in np.arange(1,ap_num+1):
            label = "spec"+str(ap)+" = "
            pos_list.append(string.find(label))
        pos_list.append(-1)

        multispec_items = {}
        for ap in np.arange(1,ap_num+1):
            pos1 = pos_list[ap-1]
            pos2 = pos_list[ap]
            tmp  = string[pos1:pos2].split("=")[1].replace("\"","").strip()
            multispec_items[ap]=MultiSpecItem(tmp)
        return multispec_items

    def get_wv(self):
        ''' get wavelength for a record in multispec'''
        if self.dtype == 0:
            # linear coordinate
            wv = self.w1 + np.arange(self.nw)*self.dw

        elif self.dtype == 2:
            # nonlinear coordiante
            p   = np.arange(self.pmin[0], self.pmax[0]+1e-6)
            n   = (p - (self.pmax[0] + self.pmin[0])/2.0)/((
                        self.pmax[0] - self.pmin[0])/2.0)
            wv  = np.zeros_like(n)
            for func in self.ftype.keys():
                if int(self.ftype[func])==1:
                    # ftype = 1: Chebyshev Polynomial
                    wvi = chebyshev_poly(n,self.c[func])
                elif int(self.ftype[func])==2:
                    # ftype = 2: Legendre Polynomial
                    wvi = legendre_poly(n,self.c[func])
                elif int(self.ftype[func])==3:
                    # ftype = 3: Cubic Spline
                    pass
                wv += self.wt[func]*wvi
        return wv
    @staticmethod
    def calib(filename):

        data,head = pf.getdata(filename,header=True)
        multispec_items = MultiSpecItem.get_wat2(head)

        wv_group   = {}
        flux_group = {}

        for N in multispec_items.keys():
            ord = multispec_items[N].beam
            wv = multispec_items[N].get_wv()
            if multispec_items[N].dtype==0:
                flux = data[N-1,:]
            elif multispec_items[N].dtype==2:
                pmin = multispec_items[N].pmin[0]
                pmax = multispec_items[N].pmax[0]
                flux = data[N-1,int(pmin)-1:int(pmax)-int(pmin)+1]

            wv_group[ord]   = wv
            flux_group[ord] = flux

        return wv_group,flux_group

class MultiSpecItem(object):
    """Object for MultiSpec Format"""

    def __init__(self,string):
        self.string = string
        self.parse_string()

    def parse_string(self):
        # fix "xxxE-40.xxx problem"

        g = self.string.split(" ")
        for i in np.arange(len(g)):
            if g[i].count('.')==2:
                ## find first dot
                #pos1 = e.find(".")
                ## find second dot
                #pos2 = e[pos1+1:].find(".") + pos1 + 1
                ## replace the second dot
                #g[i] = g[i].replace('0.',' 0.')

                t = g[i].split('.')
                t[2] = t[1][-1]+'.'+t[2]
                t[1] = t[1][:-1]
                g[i] = t[0]+'.'+t[1]+' '+t[2]

        string = ' '.join(g)

        g = string.split()
        self.ap    =   int(g[0])
        self.beam  =   int(g[1])
        self.dtype =   int(g[2])
        self.w1    = float(g[3])
        self.dw    = float(g[4])
        self.nw    =   int(g[5])
        self.z     = float(g[6])
        self.aplow = float(g[7])
        self.aphigh= float(g[8])

        if self.dtype == 0:
            # linear coordinate
            pass

        elif self.dtype == 2:
            # nonlinear coordiante
            func = 0
            self.wt    = {}
            self.wt0   = {}
            self.ftype = {}
            self.order = {}
            self.pmin  = {}
            self.pmax  = {}
            self.c     = {}
            pos       = 0
            while(pos+9<len(g)):
                self.wt[func]    = float(g[ 9+pos])
                self.wt0[func]   = float(g[10+pos])
                self.ftype[func] = round(float(g[11+pos]))
                if self.ftype[func] in [1,2]:
                    # Chebyshev (ftype=1) or Legendre (ftype=2) Polynomial
                    self.order[func] = int(g[12+pos])
                    self.pmin[func]  = float(g[13+pos])
                    self.pmax[func]  = float(g[14+pos])
                    self.c[func]     = []
                    for j in np.arange(self.order[func]):
                        self.c[func].append(float(g[15+pos+j]))
                    pos += 3+3+self.order[func]
                func += 1

    @staticmethod
    def get_wat2(head):
        ''' get WAT2 string, return a turple, consist of
        MultiSpecItem instances '''
        string = ""
        for item in head.items():
            if item[0][0:4]=="WAT2":
                string += item[1].ljust(68)
        ap_num = head['NAXIS2']
        g = string.split('spec')
        pos_list = []
        for ap in np.arange(1,ap_num+1):
            label = "spec"+str(ap)+" = "
            pos_list.append(string.find(label))
        pos_list.append(-1)

        multispec_items = {}
        for ap in np.arange(1,ap_num+1):
            pos1 = pos_list[ap-1]
            pos2 = pos_list[ap]
            tmp  = string[pos1:pos2].split("=")[1].replace("\"","").strip()
            multispec_items[ap]=MultiSpecItem(tmp)
        return multispec_items

    def get_wv(self):
        ''' get wavelength for a record in multispec'''
        if self.dtype == 0:
            # linear coordinate
            wv = self.w1 + np.arange(self.nw)*self.dw

        elif self.dtype == 2:
            # nonlinear coordiante
            p   = np.arange(self.pmin[0], self.pmax[0]+1e-6)
            n   = (p - (self.pmax[0] + self.pmin[0])/2.0)/((
                        self.pmax[0] - self.pmin[0])/2.0)
            wv  = np.zeros_like(n)
            for func in self.ftype.keys():
                if int(self.ftype[func])==1:
                    # ftype = 1: Chebyshev Polynomial
                    wvi = chebyshev_poly(n,self.c[func])
                elif int(self.ftype[func])==2:
                    # ftype = 2: Legendre Polynomial
                    wvi = legendre_poly(n,self.c[func])
                elif int(self.ftype[func])==3:
                    # ftype = 3: Cubic Spline
                    pass
                wv += self.wt[func]*wvi
        return wv

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


def overscan(logfile=None, configfile=None):

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

            #pf.writeto('bias.fits',bias)
            ovs_filename = '.'.join(filename.split('.')[0:-1])+'_ovr.fits'
            pf.writeto(ovs_filename,frame,head, overwrite=True)
            print('-> %s'%(ovs_filename))
