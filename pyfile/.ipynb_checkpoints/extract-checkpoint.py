#!/usr/bin/env python
import os, sys, math, datetime, dateutil

import copy

import xml.dom as dom
from xml.dom import minidom

import numpy as np
from astropy.io import fits as pf
import scipy.interpolate as intp
import scipy.optimize    as opt
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from sys import exit
import matplotlib 
matplotlib.use("TkAgg") 

import pylab


#-----------------------------------------------------------------
# classes and functions for obslog
#-----------------------------------------------------------------


OBJECT_LST = ['Bias','Flat','Comp','Dark']
CALIB_LST  = ['ThAr','Iodn','Mo','FeAr']

class LogItem(object):

    '''
        .date : = log.date, a datetime.date instance
        .time : a string with standard time description
        .datetime: a datetime.datetime instance

    '''
    def __init__(self,**kwargs):
        for name in kwargs:
            value = kwargs[name]
            object.__setattr__(self,name,value)

    def get_filename(self,filename_composition=None):
        '''get filename for given item id
        %Y Year with century as a decimal number.
        %m Month as a decimal number [01,12].
        %d Day of the month as a decimal number [01,31].
        %iN Frame ID as a decimal number, with length of N filled by zero.
        '''
        if filename_composition == None:
            filename_composition = self.parent.filename_composition
        fn_lst = []
        i=0
        while(True):
            if i==len(filename_composition):
                break
            elif filename_composition[i]=='%':
                key = filename_composition[i+1]
                if key == 'Y':
                    fn_lst.append(str(self.date.year))
                    i += 2
                elif key == 'm':
                    fn_lst.append(str(self.date.month).rjust(2,'0'))
                    i += 2
                elif key == 'd':
                    fn_lst.append(str(self.date.day).rjust(2,'0'))
                    i += 2
                elif key == 'i':
                    id_len = int(filename_composition[i+2])
                    fn_lst.append(str(self.id).rjust(id_len,'0'))
                    i += 3
            else:
                fn_lst.append(filename_composition[i])
                i += 1
        return ''.join(fn_lst)

    def utc_datetime(self):
        g = self.timezone.split(":")
        h = int(g[0])
        m = int(g[1])
        if h<0:
            m = -m
        tz = datetime.timedelta(hours=h,minutes=m)

        if type(self.time)==type('a') and len(self.time)>0:
            return combine_datetime(self.date, self.time)-tz
        else:
            return combine_datetime(self.date, '00:00:00')

def combine_datetime(date,timestring):
    '''get actually datetime'''
    g = timestring.split(':')
    days = int(int(g[0])/24.)
    rest_hours = int(g[0])%24
    g[0] = str(rest_hours).rjust(2,'0')
    tstring = ':'.join(g)
    t = dateutil.parser.parse(tstring).time()
    dt = datetime.datetime.combine(date+datetime.timedelta(days=days),t)
    return dt

class Log(object):
    def __init__(self,**kwargs):
        self.item_list = kwargs.pop('list',[])

    def add_item(self,item):
        self.item_list.append(item)

    def add_itemlist(self,list):
        self.item_list = list

    def check_files(self):
        #check file exists
        for item in self.item_list:
           if not item.skip:
                filename = item.get_filename(self.filename_composition)
                if not os.path.exists(filename):
                    print('Error', filename, 'does not exist')
                    return False
        return True

    def show_items(self):
        for item in self.item_list:
            if not item.skip:
                print(item.id, item.date, item.object,)
                print(item.exptime, item.note, item.i2,item.skip)

    def __str__(self):
        string_lst = []
        for item in self.item_list:
            lst = [str(item.id),
                   item.object.rjust(6),
                   str(item.date),
                   str(item.time),
                   str(item.datetime),
                   str(item.exptime).rjust(8),
                  ]
            string_lst.append(' '.join(lst))
        return os.linesep.join(string_lst)

    def save_file(self,filename,object=None,exptime=None):
        file = open(filename,'w')
        for item in self.item_list:
            if not item.skip:
                if object == 'Star':
                    if ((item.object in OBJECT_LST) or
                        (item.object in  CALIB_LST)):
                        continue
                elif object != None and item.object!=object:
                    continue
                if exptime != None and abs(item.exptime-exptime)>1e-6:
                    continue
                filename = item.get_filename(self.filename_composition)
                file.write(filename+os.linesep)
        file.close()

    def get_itemlist(self,**kwargs):
        lst = []
        for item in self.item_list:
            match = True
            for key in kwargs:
                value = kwargs[key]

                # judge string
                if type(value) == type('a'):
                    this_match = getattr(item,key).lower() == value.lower()

                # judge integer or long integer
                #elif type(value) == type(1) or type(value) == type(1L):
                elif type(value) == type(1):
                    # print('a')
                    # exception: if type do not match each other
                    # item attr is float but compare value is integer
                    if type(getattr(item,key)) == type(1.0):
                        # convert value to float
                        value = float(value)
                        this_match = abs(getattr(item,key)-value) < 1e-6
                    else:
                        this_match = getattr(item,key) == value

                # judge float
                elif type(value) == type(1.0):
                    this_match = abs(getattr(item,key)-value) < 1e-6

                else:
                    this_match = False

                match = match and this_match

            if match and (not item.skip):
                lst.append(item)

        return lst

    def get_filenamelist(self,**kwargs):

        itemlst = self.get_itemlist(**kwargs)

        lst = []

        for item in itemlst:
            lst.append(item.get_filename(self.filename_composition))

        return lst

    def set_filename_composition(self,list):
        self.filename_composition=list

    def check_file_exist(self):
        label = True
        for item in self.item_list:
            if not item.skip:
                filename = item.get_filename(self.filename_composition)
                if not os.path.exists(filename):
                    print('Error:', filename, 'doe not exist')
                    label = False
        return label

    def has_object(self, object):
        find = False
        for item in self.item_list:
            if item.skip:
                continue
            if item.object == object:
                find = True
                break
        return find

def read_log(filename):

    '''read log file, return a Log instance '''

    object_dict = {
    'bias':'Bias',
    'flat':'Flat',
    'dark':'Dark',
    'i2':  'Iodn',
    'iodn':'Iodn',
    'iodine':'Iodn',
    'comp':'Comp',
    'thar':'ThAr',
    }

    logfile = open(filename)

    log = Log()

    log.info = {}

    for line in logfile:
        line = line.strip()

        # blank line
        if len(line)==0:
            continue

        # read header
        if line[0]=='%':
            line = line[1:].strip()
            g = line.split('=')

            # read columns information
            if g[0].strip().lower()=='cols':
                col_lst = []
                tmp_g = g[1].split(',')
                for e in tmp_g:
                    col_lst.append(e.strip().lower())

            # read obs date
            elif g[0].strip().lower()=='date':
                date = g[1].strip()
                log.date = dateutil.parser.parse(date).date()

            # read filename composition
            elif g[0].strip().lower()=='filename composition':
                log.filename_composition = g[1].strip()

            # read program name
            elif g[0].strip().lower()=='program':
                log.program = g[1].strip()

            # read observatory
            elif g[0].strip().lower()=='observatory':
                log.observatory = g[1].strip()

            # read telescope
            elif g[0].strip().lower()=='telescope':
                log.telescope = g[1].strip()

            # read instrument
            elif g[0].strip().lower()=='instrument':
                log.instrument = g[1].strip()

            # read detector
            elif g[0].strip().lower()=='detector':
                log.detector = g[1].strip()

            # read observer
            elif g[0].strip().lower()=='observer':
                log.observer = g[1].strip()

            # read operator
            elif g[0].strip().lower()=='operator':
                log.operator = g[1].strip()

            # read timezone
            elif g[0].strip().lower()=='time zone':
                log.timezone = g[1].strip()

            else:
                log.info[g[0].strip()]=g[1].strip()

        #elif line[0]!='#':
        else:
            g = line.split('|')
            if g[0].strip()=='#':
                skip = True
            else:
                skip = False

            pos = col_lst.index('id')
            idg = g[pos+1].split(',')
            id_lst = []
            for idge in idg:
                ids = idge.split('-')
                if len(ids)==2:
                    for i in range(int(ids[0]),int(ids[1])+1):
                        id_lst.append(i)
                else:
                    id_lst.append(int(idge))

            # generate tmp item list
            item_lst = []
            for id in id_lst:
                item = LogItem(id=id,date=log.date,skip=skip,parent=log)
                item_lst.append(item)

            # find object
            pos = col_lst.index('object')
            obj_str = g[pos+1]
            if len(id_lst)>1:
                #print(x)
                obj_g = obj_str.split('x')
                print(obj_g)
                object = obj_g[0].strip()
                count = int(obj_g[1])
                if count != len(id_lst):
                    print("Error: count and id don't match",obj_str,id_lst)
            else:
                object = obj_str.strip()

            if object.strip().lower() in object_dict:
                object = object_dict[object.strip().lower()]

            for item in item_lst:
                item.object = object

            # find time
            if 'time' in col_lst:
                pos = col_lst.index('time')
                time = g[pos+1].strip()
                for item in item_lst:
                    item.time = time
                    # item.time is a string
                    if time.strip()!='':
                        item.datetime = combine_datetime(log.date,item.time)
                    else:
                        item.datetime = combine_datetime(log.date,'00:00:00')
                    item.timezone = log.timezone

            # find exptime
            pos = col_lst.index('exptime')
            if object == 'Bias':
                exptime = 0.0
            elif g[pos+1].strip()=='':
                exptime = None
            else:
                exptime = float(g[pos+1])
            for item in item_lst:
                item.exptime = exptime

            # find note
            if 'note' in col_lst:
                pos = col_lst.index('note')
                note = g[pos+1].strip()
                if 'with i2' in note.lower():
                    i2 = True
                    p = note.lower().index('with i2')
                    l = len('with i2')
                    note = note[:p]+note[p+l:]
                elif 'without i2' in note.lower():
                    i2 = False
                    p = note.lower().index('without i2')
                    l = len('without i2')
                    note = note[:p]+note[p+l:]
                else:
                    i2 = None
                note = note.strip()
                if len(note)==0:
                    note = None

                for item in item_lst:
                    item.note = note
                    item.i2 = i2

            for item in item_lst:
                log.add_item(item)
                #print item.date, item.id, item.object, item.time, item.exptime, item.note, item.i2
    logfile.close()

    return log

def find_logfile():

    '''find log file with surfix '.obslog' in current directory'''

    log_lst = []
    surfix  = '.obslog'

    for filename in os.listdir(os.curdir):
        if filename[-len(surfix):] == surfix:
            log_lst.append(filename)

    if len(log_lst)==1:
        return log_lst[0]
    else:
        print(log_lst)
        return None

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
#-----------------------------------------------------------------
# classes and functions for math
#-----------------------------------------------------------------

def poly(x,coeffs):
    '''return values for polynomial
    x is either a value or numpy array, coeffs=[c0,c1,c2...cn],
    return c0 + c1*x + c2*x^2 + ... + cn*x^n,
    '''
    # c0 + c1*x + c2*x^2 + ... + cn*x^n =
    # ((C[n]*x + C[n-1])*x + C[n-2])*x ... +C[0]
    # much faster than normal method
    l = len(coeffs)
    res = coeffs[l-1]
    for i in np.arange(l-1):
        res = res * x + coeffs[l-i-2]
    return res
    # res = 0.0
    # for i in np.arange(len(coeffs)):
    #     res += coeffs[i]*np.power(x,i)
    # return res

def chebyshev_poly(x,coeffs):
    """ return values for Fist Kind Chebyshev Polynomial
    x is either a value or numpy array, coeffs=[c0,c1,c2...cn],
    return c0*P[0] + c1*P[1] + c2*P[2] + ... + cn*P[n],
    where ord[i] is defined as:
        P[0] = 1
        P[1] = x
        P[2] = 2*x^2 - 1
        ......
        P[i] = 2*x*P[i-1] - P[i-2]
    """
    # append 0 for order 1 if only order 0
    if len(coeffs)==1:
        coeffs = np.array(coeffs)
        coeffs = np.append(coeffs,0)

    P   = [1.0,x]
    res = 0.0

    while(len(P)<len(coeffs)):
        P.append(2*x*P[-1]-P[-2])

    for i in np.arange(len(coeffs)):
        res += coeffs[i]*P[i]
    return res

def chebyshev2_poly(x,coeffs):
    """ return values for Second Kind Chebyshev Polynomial
    x is either a value or numpy array, coeffs=[c0,c1,c2...cn],
    return c0*P[0] + c1*P[1] + c2*P[2] + ... + cn*P[n],
    where ord[i] is defined as:
        P[0] = 1
        P[1] = 2x
        P[2] = 2*x^2 - 1
        ......
        P[i] = 2*x*P[i-1] - P[i-2]
    """
    # append 0 for order 1 if only order 0
    if len(coeffs)==1:
        coeffs = np.array(coeffs)
        coeffs = np.append(coeffs,0)

    P   = [1.0,2.0*x]
    res = 0.0

    while(len(P)<len(coeffs)):
        P.append(2*x*P[-1]-P[-2])

    for i in np.arange(len(coeffs)):
        res += coeffs[i]*P[i]
    return res


def legendre_poly(x,coeffs):
    """ return values for Legendre Polynomial
    x is either a value or numpy array, coeffs=[c0,c1,c2...cn],
    return c0*P[0] + c1*P[1] + c2*P[2] + ... + cn*P[n],
    where P[i] is defined as:
        P[0] = 1
        P[1] = x
        P[2] = 1/2*(3*x^2 - 1)
        ......
        P[i] = ( (2*i-1)*x*P[i-1] - (i-1)*P[i-2] )/i
    """
    # append 0 for order 1 if only order 0
    if len(coeffs)==1:
        coeffs = np.array(coeffs)
        coeffs = np.append(coeffs,0)

    P   = [1.0,x]
    res = 0.0

    i = 2
    while(len(P)<len(coeffs)):
        P.append(((2*i-1)*x*P[-1]-(i-1)*P[-2])/i)
        i += 1

    for i in np.arange(len(coeffs)):
        res += coeffs[i]*P[i]
    return res


#-----------------------------------------------------------------


def delete_fits(filename):
    if filename[-5:]=='.fits' and os.path.exists(filename):
        os.remove(filename)
    elif filename[-4:]=='.lst':
        file = open(filename)
        for row in file:
            row = row.strip()
            if row[-5:]=='.fits' and os.path.exists(row):
                os.remove(row)
        file.close()

def prepare_lst(reflst, surfix, filename):
    file = open(filename,'w')
    fileref = open(reflst)
    for row in fileref:
        row = row.strip()
        if row[-5:]=='.fits':
            file.write('%s_%s.fits%s'%(row[0:-5],surfix,os.linesep))
    fileref.close()
    file.close()

def check_data():
    ''' check data -
    the existences of frames in the log file
    the exposure times
    the offset between log file and control computer
    '''

    logfile = find_logfile()
    log = read_log(logfile)
    print('Obs Log File =',logfile)

    if not log.check_file_exist():
        exit()

    prev_time_end = None
    prev_fn = None

    time_offset_lst = []

    for item in log.item_list:
        if not item.skip:
            fn = item.get_filename(log.filename_composition)

            # get fits header
            head = pf.getheader(fn)
            exptime  = head['EXPTIME']
            time_sta = dateutil.parser.parse(head['DATE-STA'])
            time_end = dateutil.parser.parse(head['DATE-END'])

            # check exptime
            caltime = (time_end-time_sta).seconds
            if abs(caltime - exptime)>1.1:
                print(' Warning: exptime for',fn,':',)
                print('t2-t1 =',caltime,'sec, exptime =',exptime,'sec')

            # check start time
            if prev_time_end != None:
                dtime = time_sta - prev_time_end
                if time_sta < prev_time_end:
                    print(' Error: Start time of',fn,'earlier than',prev_fn)

            # determine time offset
            if item.time != '':
                dtime = item.datetime - time_sta
                time_offset_lst.append(dtime)

            prev_time_end = time_end
            prev_fn       = fn

    # find median dtime value
    time_offset_lst.sort()
    time_offset = time_offset_lst[int(len(time_offset_lst)/2)]
    if time_offset > datetime.timedelta(0,0):
        log_early = True
        print('Mean Time Offset: log = fits +', abs(time_offset))
    else:
        log_early = False
        print('Mean Time Offset: log = fits -', abs(time_offset))

    for item in log.item_list:
        if not item.skip:
            fn = item.get_filename(log.filename_composition)

            # get fits header
            head = pf.getheader(fn)
            exptime  = head['EXPTIME']
            time_sta = dateutil.parser.parse(head['DATE-STA'])

            if item.time != '':
                #corrected_ctime = item.datetime - time_offset
                #abstd = abs(corrected_ctime-time_sta)
                #if abstd > datetime.timedelta(0,1):
                #    print 'Check log time of',item.id,':',
                #    if corrected_ctime > time_sta:
                #        print abstd,'faster'
                #    else:
                #        print abstd,'slower'
                print(item.id,)
                if log_early:
                    print('LOG - FITS =',item.datetime - time_sta)
                else:
                    print('FITS - LOG =',time_sta - item.datetime)

def plot_flat():
    '''
    plot a portion of flat frame to check the brightness
    variations
    '''

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import colorConverter

    r1 = 1100
    r2 = 1200

    logfile = find_logfile()
    log = read_log(logfile)
    print('Obs Log File =',logfile)

    xs = np.arange(0,r2-r1)
    zs = []
    verts = []
    fcolors = []

    cc = lambda arg: colorConverter.to_rgba(arg,alpha=0.4)

    for item in log.item_list:
        if not item.skip and item.object=='Flat':
            filename = item.get_filename(log.filename_composition)
            datai = pf.getdata(filename)
            f = datai[r1:r2,1650]
            f[0],f[-1]=0,0
            zs.append(item.id)
            verts.append(zip(xs,f))
            fcolors.append(cc('b'))

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    poly = PolyCollection(verts,facecolors = fcolors)
    ax.add_collection3d(poly,zs=zs,zdir='y')

    ax.set_xlim3d(0,r2-r1)
    ax.set_ylim3d(zs[0]-1,zs[-1]+1)
    ax.set_zlim3d(0,50000)

    ax.set_xlabel('Rows - '+str(r1))
    ax.set_ylabel('ID')
    ax.set_zlabel('Flux')

    plt.show()
    plt.close(fig)


def overscan():

    ''' overscan corrections'''

    conf = load_conf('HRS.conf')
    ccd_map = conf['ccd_map']

    readout_lst = get_readout_lst(ccd_map)
    large_rect = get_large_rect(readout_lst)

    # weight list
    w_lst = np.array([1,3,5,3,1])
    scan_p = len(w_lst)

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
            pf.writeto(ovs_filename,frame,head)
            print('-> %s'%(ovs_filename))

def extract1d(calib=False):
    import matplotlib 
    matplotlib.use("TkAgg") 
    from pyraf import iraf

    # remove the old iraf.log file
    if not calib and os.path.exists('iraf.log'):
        os.remove('iraf.log')

    if not calib:
        os.system('mkiraf')

    logfile = find_logfile()
    log = read_log(logfile)
    print('Obs Log File =',logfile)

    # initialize IRAF
    print( )
    print('Initialize IRAF')
    print( )
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.echelle(_doprint=0)

    has_iodn = log.has_object('Iodn')

    # make list of different types
    log.save_file('thar.lst', object='ThAr')
    log.save_file('star.lst', object='Star')
    if has_iodn:
        log.save_file('iodn.lst', object='Iodn')

    # make flat list of different exposure times
    conf = load_conf('HRS.conf')
    flat_exptime = conf['flat_exptime']
    apsize = conf['apsize']

    for row in flat_exptime:
        color   = row[0]
        exptime = row[1]
        flat_lst = log.get_filenamelist(
                   object  = 'Flat',
                   exptime = exptime
                   )

        # write filenames to list
        flat_lst_file = open('flat_%s.lst'%color,'w')
        flat_lst_file.write(os.linesep.join(flat_lst))
        flat_lst_file.close()

        prepare_lst('flat_%s.lst'%color,'ovr','flat_%s_ovr.lst'%color)

        if calib and os.path.exists('flat_%s.lst'%color):
            continue

        # Combine flat of different colors
        print( )
        print('Combine %d %s Flats (EXPTIME = %s seconds)'%(len(flat_lst),color.title(),str(float(exptime))))
        print( )
        iraf.imcomb.unlearn()
        iraf.imcomb.combine = 'average'
        iraf.imcomb.logfile = 'iraf.log'
        delete_fits('flat_%s.fits'%color)
        iraf.imcomb('@flat_%s_ovr.lst'%color, 'flat_%s.fits'%color)

    if calib and os.path.exists('flat.fits'):
        pass
    else:
        # mosaic flats
        print( )
        print('Mosaic Flats')
        print( )
        mosaic_flat('flat.fits')

    # locating order
    print( )
    print('Order Location')
    print( )
    iraf.echelle.unlearn()
    iraf.echelle.dispaxis = 1

    iraf.apdefault.unlearn()
    iraf.apdefault.lower = -apsize
    iraf.apdefault.upper =  apsize

    iraf.apedit.unlearn()
    iraf.apedit.interac = 'yes'
    iraf.apedit.find    = 'no'
    iraf.apedit.recente = 'no'
    iraf.apedit.resize  = 'no'
    iraf.apedit.edit    = 'yes'
    iraf.apedit.nsum    = 10
    iraf.apedit.width   = 20.0
    iraf.apedit.radius  = 10.0
    if not calib:
        iraf.apedit('flat.fits')

    iraf.aptrace.unlearn()
    iraf.aptrace.interac = 'yes'
    iraf.aptrace.find    = 'yes'
    iraf.aptrace.recente = 'no'
    iraf.aptrace.resize  = 'no'
    iraf.aptrace.edit    = 'no'
    iraf.aptrace.trace   = 'yes'
    iraf.aptrace.fittrac = 'yes'
    iraf.aptrace.functio = 'legendre'
    iraf.aptrace.order   = 5
    if not calib:
        iraf.aptrace('flat.fits')

    print( )
    print('Flat Fielding Correction')
    print( )
    iraf.apresize.unlearn()
    iraf.apresize.llimit = -apsize
    iraf.apresize.ulimit =  apsize
    iraf.apresize.ylevel = 'INDEF'
    iraf.apresize.bkg    = 'no'

    iraf.apflatten.unlearn()
    iraf.apflatten.find    = 'no'
    iraf.apflatten.trace   = 'no'
    iraf.apflatten.fittrac = 'no'
    iraf.apflatten.functio = "legendre"
    iraf.apflatten.order   = 8
    iraf.apflatten.niterat = 5
    if not calib:
        delete_fits('flat_n.fits')
        iraf.apflatten('flat.fits','flat_n.fits')

    iraf.imarith.unlearn()

    # parse thar
    prepare_lst('thar.lst', 'ovr', 'thar_ovr.lst')
    prepare_lst('thar.lst', 'flt', 'thar_flt.lst')
    if not calib:
        delete_fits('thar_flt.lst')
        iraf.imarith('@thar_ovr.lst','/','flat_n.fits','@thar_flt.lst')

    # parse star
    prepare_lst('star.lst', 'ovr', 'star_ovr.lst')
    prepare_lst('star.lst', 'flt', 'star_flt.lst')
    if not calib:
        delete_fits('star_flt.lst')
        iraf.imarith('@star_ovr.lst','/','flat_n.fits','@star_flt.lst')

    # parse iodine
    if has_iodn:
        prepare_lst('iodn.lst', 'ovr', 'iodn_ovr.lst')
        prepare_lst('iodn.lst', 'flt', 'iodn_flt.lst')
        if not calib:
            delete_fits('iodn_flt.lst')
            iraf.imarith('@iodn_ovr.lst','/','flat_n.fits','@iodn_flt.lst')

    print( )
    print('Background Correction')
    print( )
    iraf.apresize.unlearn()
    iraf.apresize.llimit = -apsize
    iraf.apresize.ulimit =  apsize
    iraf.apresize.ylevel = 'INDEF'
    iraf.apresize.bkg    = 'no'

    iraf.apscatter.unlearn()
    iraf.apscatter.referen = 'flat.fits'
    iraf.apscatter.interac = 'no'  # Run task interactively?
    iraf.apscatter.find    = 'no'
    iraf.apscatter.recente = 'yes'
    iraf.apscatter.resize  = 'yes'
    iraf.apscatter.edit    = 'yes'
    iraf.apscatter.trace   = 'no'
    iraf.apscatter.fittrac = 'no'
    iraf.apscatter.subtrac = 'yes'
    iraf.apscatter.smooth  = 'yes'
    iraf.apscatter.fitscat = 'no'  # Fit scattered light interactively?
    iraf.apscatter.fitsmoo = 'no'  # Smooth the scattered light interactively?
    iraf.apscatter.nsum    = 100
    iraf.apscatter.buffer  = 2

    iraf.apscat1.unlearn()
    iraf.apscat1.function    = 'spline3'
    iraf.apscat1.order       = 4
    iraf.apscat1.low_reject  = 5.0
    iraf.apscat1.high_reject = 2.0
    iraf.apscat1.niterate    = 5

    iraf.apscat2.unlearn()
    iraf.apscat2.function    = 'spline3'
    iraf.apscat2.order       = 5
    iraf.apscat2.low_reject  = 3.0
    iraf.apscat2.high_reject = 3.0
    iraf.apscat2.niterate    = 5

    prepare_lst('star.lst', 'bkg', 'star_bkg.lst')
    if not calib:
        delete_fits('star_bkg.lst')
        iraf.apscatter('@star_flt.lst', '@star_bkg.lst')

    if has_iodn:
        prepare_lst('iodn.lst', 'bkg', 'iodn_bkg.lst')
        if not calib:
            delete_fits('iodn_bkg.lst')
            iraf.apscatter('@iodn_flt.lst', '@iodn_bkg.lst')

    print( )
    print('Extrcat 1d Spectra')
    print( )

    # optimal extraction or notinput
    prompt = 'Do you use optimal extraction or not? [y/N]: '
    while(True):
        ans = input(prompt)
        ans = ans.strip().lower()
        if len(ans)==0 or ans in ['n','no']:
            weights = 'none'
            break
        elif ans in ['y','yes']:
            weights = 'variance'
            break
        else:
            continue

    iraf.apsum.unlearn()
    iraf.apsum.output  = '@star_sum.lst'
    iraf.apsum.format  = 'echelle'
    iraf.apsum.referen = 'flat.fits'
    iraf.apsum.interac = 'no'    # Run task interactively?
    iraf.apsum.find    = 'no'    # Find apertures?
    iraf.apsum.recente = 'no'    # Recenter apertures?
    iraf.apsum.resize  = 'no'    # Resize apertures?
    iraf.apsum.edit    = 'no'    # Edit apertures?
    iraf.apsum.trace   = 'no'    # Trace apertures?
    iraf.apsum.fittrac = 'no'    # Fit the traced points interactively?
    iraf.apsum.extract = 'yes'   # Extract apertures?
    iraf.apsum.extras  = 'no'    # Extract sky, sigma, etc.?
    iraf.apsum.review  = 'no'    # Review extractions?
    iraf.apsum.weights = weights # Extraction weights (none|variance)
    prepare_lst('star.lst','sum','star_sum.lst')
    if not calib:
        delete_fits('star_sum.lst')
        iraf.apsum('@star_bkg.lst')

    iraf.apsum.output  = '@thar_sum.lst'
    prepare_lst('thar.lst','sum','thar_sum.lst')
    if not calib:
        delete_fits('thar_sum.lst')
        iraf.apsum('@thar_flt.lst')

    if has_iodn:
        iraf.apsum.output  = '@iodn_sum.lst'
        prepare_lst('iodn.lst','sum','iodn_sum.lst')
        if not calib:
            delete_fits('iodn_sum.lst')
            iraf.apsum('@iodn_bkg.lst')

    print( )
    print('Wavelength Calibration')
    print( )

    prompt = 'Do you want to use an [E]xisting reference or identify a [N]ew ThAr image? [E/N]: '
    while(True):
        ans = input(prompt)
        ans = ans.strip().lower()
        if len(ans)==0 or ans not in ['e','n']:
            continue
        elif ans == 'e':
            new_ecid = False
            break
        elif ans == 'n':
            new_ecid = True
            break
        else:
            continue


    if new_ecid:
        # delete the existing ec files corresponding to the filenames
        # in thar_sum.lst
        file = open('thar_sum.lst')
        for row in file:
            fn = row.strip()
            ecfn = 'database/ec'+fn.replace('.fits','')
            if os.path.exists(ecfn):
                os.remove(ecfn)
        file.close()

        # list files and their exptimes in thar_sum.lst
        ref_lst = {}
        no = 1
        file = open('thar_sum.lst')
        for row in file:
            row = row.strip()
            if len(row)==0 or row[0] in ['#'] or row[-5:]!='.fits':
                continue
            head = pf.getheader(row)
            exptime = head['EXPTIME']
            ref_lst[no] = row
            print(' [%02d] %s %6.2f s'%(no,row,exptime))
            no += 1
        file.close()

        # Now the user have to choose a thar file to identify wavelengths
        while (True):
            n = input('Select a ThAr frame for wavelength identification [1]:')
            try:
                n = int(n)
                if 0 < n <= len(ref_lst):
                    break
            except:
                pass

        ref_thar = ref_lst[int(n)]
        ref_name = ref_thar.replace('.fits','')

        iraf.ecid.unlearn()
        iraf.ecid.maxfeat = 100
        iraf.ecid(ref_thar)
    else:
        prompt = 'Cannot find reference in database.\n'
        prompt += 'Put an `ec` file in database/ and press [Enter] to refresh: '
        while(True):
            ref_lst = {}
            no = 1
            lst = os.listdir('./database')
            lst.sort()
            for fname in lst:
                if fname[0:2]=='ec':
                    ref_lst[no] = fname

                    # get information of this ec file
                    res = get_ec_info('./database/%s'%fname)
                    wvstd, rvstd    = res[0], res[1]
                    nused, nlines   = res[2], res[3]
                    norders         = res[4]
                    wv_min, wv_max  = res[5], res[6]
                    date = '%s-%s-%s'%(fname[2:6],fname[6:8],fname[8:10])

                    print(' [%d] %18s (%10s) %7.5f %6.3f %4d/%4d %3d (%7.1f - %7.1f)'%(
                            no, ref_lst[no][2:], date,
                            wvstd, rvstd,
                            nused, nlines,
                            norders,
                            wv_min, wv_max
                            ))
                    no += 1
            if len(ref_lst)==0:
                _ = input(prompt)
            else:
                break
        prompt = 'Which frame do you want to choose as the reference? [1]:'
        while(True):
            n = input(prompt)
            try:
                n = int(n)
                if 0 < n <= len(ref_lst):
                    ref_name = ref_lst[n][2:]
                    break
                else:
                    continue
            except:
                pass

    iraf.ecreid.unlearn()
    iraf.ecreid.logfile = 'STDOUT, iraf.log'
    iraf.ecreid('@thar_sum.lst', ref_name)

    _ = input('Press [Enter] to continue: ')

    iraf.refsp.unlearn()
    iraf.refsp.referen = '@thar_sum.lst'
    iraf.refsp.sort    = 'DATE-STA'
    iraf.refsp.group   = ''
    iraf.refsp.time    = 'yes'
    iraf.refsp.logfile = 'STDOUT, iraf.log'
    iraf.refsp('@star_sum.lst')
    iraf.refsp('@thar_sum.lst')
    if has_iodn:
        iraf.refsp('@iodn_sum.lst')

    iraf.dispcor.unlearn()
    iraf.dispcor.lineari = 'no'
    iraf.dispcor.flux    = 'no'
    iraf.dispcor.logfile = 'iraf.log'
    prepare_lst('star.lst', '1ds', 'star_1ds.lst')
    delete_fits('star_1ds.lst')
    iraf.dispcor('@star_sum.lst', '@star_1ds.lst')

    prepare_lst('thar.lst', '1ds', 'thar_1ds.lst')
    delete_fits('thar_1ds.lst')
    iraf.dispcor('@thar_sum.lst', '@thar_1ds.lst')

    if has_iodn:
        prepare_lst('iodn.lst', '1ds', 'iodn_1ds.lst')
        delete_fits('iodn_1ds.lst')
        iraf.dispcor('@iodn_sum.lst', '@iodn_1ds.lst')


    if calib:
        # if in calibration mode, force deleting blz files
        for item in log.item_list:
            if item.skip:
                continue
            filename = item.get_filename(log.filename_composition)
            fname = '%s_%s.fits'%('.'.join(filename.split('.')[0:-1]),'blz')
            if os.path.exists(fname):
                os.remove(fname)

    prompt = 'Correct the blaze function? [Y/N]: '
    while(True):
        ans = input(prompt)
        ans = ans.strip().lower()
        if ans in ['y','yes']:
            do_blaze = True
            break
        elif ans in ['n','no']:
            do_blaze = False
            break
        else:
            continue

    if do_blaze:

        print( )
        print('Correct Blaze Functions')
        print( )
        delete_fits('flat_f.fits')
        delete_fits('flat_a.fits')
        delete_fits('flat_as.fits')

        iraf.imarith('flat.fits','/','flat_n.fits','flat_f.fits')
        iraf.apsum.output  = 'flat_a.fits'
        iraf.apsum.interac = 'no'
        iraf.apsum.review  = 'no'
        iraf.apsum('flat_f.fits')

        iraf.continuum.unlearn()
        iraf.continuum.ask         = 'YES'
        iraf.continuum.type        = 'fit'
        iraf.continuum.wavescale   = 'no'
        iraf.continuum.logfiles    = 'iraf.log'
        iraf.continuum.interactive = 'no'
        iraf.continuum.order       = 6
        iraf.continuum.low_reject  = 3.0
        iraf.continuum.high_reject = 3.0
        iraf.continuum('flat_a.fits', 'flat_as.fits')

        prepare_lst('star.lst', 'blz', 'star_blz.lst')
        delete_fits('star_blz.lst')
        iraf.imarith('@star_1ds.lst','/','flat_as.fits','@star_blz.lst')

        if has_iodn:
            prepare_lst('iodn.lst', 'blz', 'iodn_blz.lst')
            delete_fits('iodn_blz.lst')
            iraf.imarith('@iodn_1ds.lst','/','flat_as.fits','@iodn_blz.lst')

    # delete deprecated files
    for row in flat_exptime:
        color = row[0]
        os.remove('flat_%s_ovr.lst'%color)

    os.remove('star_ovr.lst')
    os.remove('star_flt.lst')
    os.remove('star_bkg.lst')
    os.remove('star_sum.lst')
    os.remove('star_1ds.lst')
    if os.path.exists('star_blz.lst'):
        os.remove('star_blz.lst')

    os.remove('thar_ovr.lst')
    os.remove('thar_flt.lst')
    os.remove('thar_sum.lst')
    os.remove('thar_1ds.lst')

    if has_iodn:
        os.remove('iodn_ovr.lst')
        os.remove('iodn_flt.lst')
        os.remove('iodn_bkg.lst')
        os.remove('iodn_sum.lst')
        os.remove('iodn_1ds.lst')
        if os.path.exists('iodn_blz.lst'):
            os.remove('iodn_blz.lst')

    # update fits headers
    update_header()

    print('Finished')


def update_header():
    logfile = find_logfile()
    log = read_log(logfile)

    for item in log.get_itemlist():
        if item.skip:
            continue
        filename = item.get_filename(log.filename_composition)
        for key in ['1ds','blz']:
            fname = '%s_%s.fits'%('.'.join(filename.split('.')[0:-1]),key)
            if os.path.exists(fname):
                f = pf.open(fname,mode='update')
                print('fname:',fname)
                hd = f[0].header
                #print('hd:',hd)
                print('OBJECT',item.object,'Name of target')
                # parse object
                hd['OBJECT']=(item.object,'Name of target')
                #hd.update('OBJECT',item.object,'Name of target')

                # update observation infomation
                hd['TELESCOP']=(log.observatory+' '+log.telescope,'Telescope')
                hd['INSTRUME']=(log.instrument,'Instrument')
                hd['DETECTOR']=(log.detector,'Detector')
                hd['OBSERVER']=(log.observer,'Observers')
                hd['OPERATOR']=(log.operator,'Telescope Operator')

                #hd.update('TELESCOP',log.observatory+' '+log.telescope,'Telescope')
                #hd.update('INSTRUME',log.instrument,'Instrument')
                #hd.update('DETECTOR',log.detector,'Detector')
                #hd.update('OBSERVER',log.observer,'Observers')
                #hd.update('OPERATOR',log.operator,'Telescope Operator')
                if item.exptime != None:
                    hd['EXPTIME']=(item.exptime,'Exposure time (second)')
                    #hd.update('EXPTIME',item.exptime,'Exposure time (second)')
                fileid = '%04d%02d%02d%03d'%(item.date.year, item.date.month, item.date.day, item.id)
                hd['FILEID']=(fileid,'File ID')
                #hd.update('FILEID',fileid,'File ID')


                # parse i2 cell
                try:
                    hd["I2CELL"]=( item.i2, "Iodine cell")
                    #hd.update("I2CELL", item.i2, "Iodine cell")
                except:
                    pass

                # parse readout
                #print(log.info)
                if "Readout" in log.info:
                #if log.info.has_key("Readout"):
                    hd["READOUT"]=(log.info["Readout"], "Readout Mode")
                    #hd.update("READOUT", log.info["Readout"], "Readout Mode")

                # parse gain
                if "Gain" in log.info:
                #if log.info.has_key("Gain"):
                    hd["GAIN"]=(float(log.info["Gain"].split()[0]), "CCD Gain (e-/ADU)")
                    #hd.update("GAIN",float(log.info["Gain"].split()[0]), "CCD Gain (e-/ADU)")

                # parse fiber
                if "Fiber Diameter" in log.info:
                #if log.info.has_key("Fiber Diameter"):
                    hd["FIBERDIA"]=(log.info["Fiber Diameter"], "Fiber diameter")
                    #hd.update("FIBERDIA", log.info["Fiber Diameter"], "Fiber diameter")

                # parese slit width
                if "Slit Width" in log.info:
                #if log.info.has_key("Slit Width"):
                    hd["SLITWID"]=(log.info["Slit Width"],"Slit width")
                    #hd.update("SLITWID", log.info["Slit Width"],"Slit width")

                # parse time
                if item.time!='':
                    lt_sta = item.datetime
                    lt_mid = item.datetime + datetime.timedelta(seconds=item.exptime/2.0)
                    lt_end = item.datetime + datetime.timedelta(seconds=item.exptime)

                    hd["OBS-STA"]=(lt_sta.isoformat(),"Start time of exposure")
                    hd["OBS-MID"]=(lt_mid.isoformat(),"Mid time of exposure")
                    hd["OBS-END"]=(lt_end.isoformat(),"End time of exposure")
                    #hd.update("OBS-STA",lt_sta.isoformat(),"Start time of exposure")
                    #hd.update("OBS-MID",lt_mid.isoformat(),"Mid time of exposure")
                    #hd.update("OBS-END",lt_end.isoformat(),"End time of exposure")

                    utc_sta = item.utc_datetime()
                    utc_mid = utc_sta + datetime.timedelta(seconds=item.exptime/2.0)
                    utc_end = utc_sta + datetime.timedelta(seconds=item.exptime)

                    hd["UTC-STA"]=(utc_sta.isoformat(),"Start time of exposure")
                    hd["UTC-MID"]=(utc_mid.isoformat(),"Mid time of exposure")
                    hd["UTC-END"]=(utc_end.isoformat(),"End time of exposure")
                    #hd.update("UTC-STA",utc_sta.isoformat(),"Start time of exposure")
                    #hd.update("UTC-MID",utc_mid.isoformat(),"Mid time of exposure")
                    #hd.update("UTC-END",utc_end.isoformat(),"End time of exposure")

                    hd["TIMESYS"]=("GPS","Time system")
                    hd["TIMEZONE"]=(log.timezone,"Time Zone")
                    #hd.update("TIMESYS","GPS","Time system")
                    #hd.update("TIMEZONE",log.timezone,"Time Zone")

                if True:
                    gain = float(log.info["Gain"].split()[0])
                    spec = load_multispec(fname.replace('_blz','_1ds'))
                    for wv in np.arange(3000,9500,500):
                        bo = spec.find_best_order(wv)
                        if bo == None:
                            continue
                        snr = spec[bo].measure_snr(gain=gain)
                        hd["SNR"+str(wv)]=(snr,"SNR around %d Angstrom"%wv)
                        #hd.update("SNR"+str(wv),snr,"SNR around %d Angstrom"%wv)

                f.flush()
                f.close()


def combine(input_lst,outputname):

    color_lst = ['r','k','b','m','c','y','g']

    spec_lst = {}
    orders = None
    for fname in input_lst:
        data = pf.getdata(fname)
        spec = {}
        for o in np.arange(data.shape[0]):
            spec[o] = data[o,:]

        #spec = load_multispec(fname)
        spec_lst[fname] = spec
        print(fname,orders)
        
        try:
            if orders == None:
                #orders = spec.orders
                orders = np.arange(data.shape[0])
                #print(orders.any())
        #if len(orders) != len(spec.orders):
        except:
            if len(orders) != data.shape[0]:
                print('Numbers of orders in %s does not math others.')
                exit()

    # default display order
    default_display_order = orders[0]
    display_order = default_display_order

    # initialize figure
    fig = plt.figure(figsize=(10,6))
    ax = fig.gca()

    def plot_order(order,xlim=None,ylim=None):
        ax.cla()
        ax.currentorder = order
        #for i,fname in enumerate(spec_lst.keys()):
        for i,fname in enumerate(input_lst):

            factor = 1./scale_lst[fname][order]

            #xdata = np.arange(spec_lst[fname][order].pts)
            xdata = np.arange(spec_lst[fname][order].size)
            #ydata = spec_lst[fname][order].flux
            ydata = spec_lst[fname][order]
            ax.plot(xdata,
                    ydata*factor,
                    color_lst[i]+'-',
                    alpha=0.4
                    )
            mask = mask_lst[fname][order]
            rmask = np.logical_not(mask)
            ax.plot(xdata[mask],
                    ydata[mask]*factor,
                    color_lst[i]+'-',
                    alpha=0.8,
                    label=u'%s \xd7 %.2f'%(fname,factor)
                    )
            ax.plot(xdata[rmask],
                    ydata[rmask]/scale_lst[fname][order],
                    color_lst[i]+'x',
                    alpha=0.6,markersize=6)
        if xlim == None:
            #ax.set_xlim(0,spec_lst[fname][order].pts-1)
            ax.set_xlim(0,spec_lst[fname][order].size-1)
        else:
            ax.set_xlim(xlim[0],xlim[1])
        if ylim == None:
            ax.set_ylim()
        else:
            ax.set_ylim(ylim[0],ylim[1])
        leg = ax.legend()
        for t in leg.get_texts():
            t.set_fontsize('small')
        #ax.set_title('Order %d'%(order))
        ax.set_title('Aperture %d'%(order))
        ax.set_xlabel('Pixel')
        ax.set_ylabel('Relative Flux')
        fig.canvas.draw()

    # find scale_lst. scales are factors that made all fluxes in the
    # same order are displayed in nearly the same level.
    scale_lst = {}
    for i, fname in enumerate(input_lst):
        scale_lst[fname] = {}
        if i == 0:
            ref_spec = spec_lst[fname]
            for o in orders:
                scale_lst[fname][o]=1.0
        else:
            for o in orders:
                #ifrom = int(ref_spec[o].pts*0.25)
                ifrom = int(ref_spec[o].size*0.25)
                #ito   = int(ref_spec[o].pts*0.75)
                ito   = int(ref_spec[o].size*0.75)
                #tmp = spec_lst[fname][o].flux / ref_spec[o].flux
                tmp = spec_lst[fname][o] / ref_spec[o]
                tmp = np.sort(tmp)
                scale_lst[fname][o]=tmp[int(tmp.size/2.0)]

    # load_mask (.msk files)
    mask_lst = load_mask(input_lst)

    # if filename not in mask_lst (no .msk), initialize mask
    for fname in input_lst:
        if fname not in mask_lst:
            mask_lst[fname] = {}
            #for o in spec_lst[fname].orders:
            for o in orders:
                #mask_lst[fname][o] = np.zeros(spec_lst[fname][o].pts)<1
                mask_lst[fname][o] = np.zeros(spec_lst[fname][o].size)<1

    plot_order(display_order)


    # key event handler
    def on_key(event):
        if event.key == 'd':
            min_dis = 1e40
            min_fname = None
            x = int(round(event.xdata))
            for fname in input_lst:
                if not mask_lst[fname][ax.currentorder][x]:
                    continue
                #dis = abs(event.ydata-spec_lst[fname][ax.currentorder].flux[x]/scale_lst[fname][ax.currentorder])
                dis = abs(event.ydata-spec_lst[fname][ax.currentorder][x]/scale_lst[fname][ax.currentorder])
                if dis<min_dis:
                    min_dis = dis
                    min_fname = fname
            if min_fname != None:
                count = 0
                for fname in input_lst:
                    count += mask_lst[fname][ax.currentorder][x]
                if count > 1:
                    mask_lst[min_fname][ax.currentorder][x]=False
            plot_order(ax.currentorder,ax.get_xlim(),ax.get_ylim())

        elif event.key == 'a':
            min_dis = 1e40
            min_fname = None
            x = int(round(event.xdata))
            for fname in input_lst:
                if mask_lst[fname][ax.currentorder][x]:
                    continue
                #dis = abs(event.ydata-spec_lst[fname][ax.currentorder].flux[x]/scale_lst[fname][ax.currentorder])
                dis = abs(event.ydata-spec_lst[fname][ax.currentorder][x]/scale_lst[fname][ax.currentorder])
                if dis<min_dis:
                    min_dis = dis
                    min_fname = fname
            if min_fname != None:
                mask_lst[min_fname][ax.currentorder][x]=True
            plot_order(ax.currentorder,ax.get_xlim(),ax.get_ylim())

        elif event.key == 'c':
            if len(input_lst)<3:
                pass
            print('calculating the cosmic ray hints...')

            ref_spec = {}
            ts_lst = {}

            for o in orders:
                #npt = spec_lst[input_lst[0]][o].pts
                npt = spec_lst[input_lst[0]][o].size
                nsp = len(input_lst)
                ifrom = int(npt*0.2)
                ito   = int(npt*0.8)

                for i, fname in enumerate(input_lst):
                    tmp = spec_lst[fname][o]
                data = np.empty((nsp,npt))
                for i,fname in enumerate(input_lst):
                    #data[i,:] = spec_lst[fname][o].flux/scale_lst[fname][o]
                    data[i,:] = spec_lst[fname][o]/scale_lst[fname][o]

                maxargs = data.argmax(axis=0)
                y,x = np.mgrid[:data.shape[0],:data.shape[1]]
                mask = (y!=maxargs)
                ref_spd = ((mask*data).sum(axis=0))/(mask.sum(axis=0))
                ref_spec[o] = ref_spd


                ts_lst[o] = {}
                for i, fname in enumerate(input_lst):
                    t = data[i,:]/ref_spd
                    tm = t.mean()
                    ts = t[ifrom:ito].std()
                    ts_lst[o][fname] = ts

                    m = data[i,:]>ref_spd*(1+3.*ts)
                    m1 = np.arange(npt)>ifrom
                    m2 = np.arange(npt)<ito
                    m3 = np.logical_and(m1,m2)
                    m = np.logical_and(m3,m)
                    mask_lst[fname][o] = np.logical_not(m)
                #    std_scale = np.sqrt(1./abs(ref_spec/(ref_spec.max())))
                #    m = data[i,:]/ref_spec > (t.mean() + 3.*std_scale*t.std())
                #    mask_lst[fname][o]=np.logical_not(m)



                #m = abs(min_spec)<min_spec.max()*1e-6
                #factor_array = np.ones_like(data)

                #factor_array = data/min_spec
                #dd = np.ma.masked_array(data,mask=mask)



                #maxargs = data.argmax(axis=0)
                #y,x = np.mgrid[:data.shape[0],:data.shape[1]]
                #mask = (y==maxargs)
                #dd = np.ma.masked_array(data,mask=mask)
                #means = dd.mean(axis=0)
                #stds  = dd.std(axis=0)
                #crs = data*mask > (means + 100.*stds)
                #crs = crs.data
                #for i,fname in enumerate(input_lst):
                #    mask_lst[fname][o] = np.logical_not(crs[i,:])

            plot_order(ax.currentorder)

        elif event.key == 'up':
            #if orders.index(ax.currentorder)!=0:
            if ax.currentorder!=orders.size-1:
                display_order = ax.currentorder + 1
                plot_order(display_order)
        elif event.key == 'down':
            #if orders.index(ax.currentorder)!=len(orders)-1:
            if ax.currentorder!=0:
                display_order = ax.currentorder - 1
                plot_order(display_order)
        else:
            pass

    fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()

    save_mask(mask_lst)

    # now combine the final spectra

    data0, head = pf.getdata(input_lst[0],header=True)
    data = np.zeros_like(data0)

    for o in orders:
        flux = np.zeros(data.shape[1])
        w    = np.zeros(data.shape[1])
        for fname in input_lst:
            #fluxi = spec_lst[fname][o].flux
            fluxi = spec_lst[fname][o]
            maski = mask_lst[fname][o]
            scalei = scale_lst[fname][o]
            flux += fluxi/math.sqrt(scalei)*maski
            w    += math.sqrt(scalei)*maski
        flux = flux/w
        scale_sum = 0.0
        for fname in input_lst:
            scale_sum += scale_lst[fname][o]
        flux = flux*scale_sum

        data[o,:]=flux
    if os.path.exists(outputname):
        os.remove(outputname)
    pf.writeto(outputname,data,head)

def load_mask(input_lst):
    '''load masks (*.msk) for every file in input_lst'''
    mask_lst = {}
    for fname in input_lst:
        maskname = fname.replace('.fits','.msk')
        if not os.path.exists(maskname):
            continue
        mask_lst[fname] = {}
        file = open(maskname)
        for row in file:
            if row[0] in ['%','#'] or len(row.strip())==0:
                continue
            g = row.split(':')
            g0 = g[0].split(',')
            order = int(g0[0])
            pts   = int(g0[1])
            mask_lst[fname][order] = np.zeros(pts)<1
            if len(g[1].strip())>0:
                g1 = g[1].split(',')
                for v in g1:
                    mask_lst[fname][order][int(v)] = False
                    # cosmic ray: False
                    # normal points: True
        file.close()
    return mask_lst

def save_mask(mask_lst):
    for fname in mask_lst:
        maskname = fname.replace('.fits','.msk')
        file = open(maskname,'w')
        for o in mask_lst[fname]:
            rmask = np.logical_not(mask_lst[fname][o])
            file.write('%d,%d:'%(o,rmask.size))
            if rmask.sum() > 0:
                lst = ['%d'%v for v in np.nonzero(rmask)[0]]
                file.write(','.join(lst))
            file.write(os.linesep)
        file.close()


def get_ec_info(filename):
    '''get information of the ec file'''

    file = open(filename)
    find = False
    for row in file:
        row = row.strip()
        if len(row)==0 or row[0]=='#':
            continue
        g = row.split()
        if g[0]=='begin':
            # initialize lists
            use_dif_lst = []
            use_pix_lst = []
            use_wv_lst  = []
            use_ord_lst = []

            rej_dif_lst = []
            rej_pix_lst = []
            rej_wv_lst  = []
            rej_ord_lst = []

            ap_lst = []
        elif g[0]=='features':
            find = True
            continue
        elif g[0]=='offset':
            find = False
            continue
        elif find:
            ap    = int(g[0])
            order = int(g[1])
            pix   = float(g[2])
            wv_fit= float(g[3])
            if g[4]=='INDEF':
                continue
            wv_map= float(g[4])
            _     = float(g[5])
            _     = int(g[6])
            use   = int(g[7])
            if use == 0:
                rej_dif_lst.append(wv_fit-wv_map)
                rej_pix_lst.append(pix)
                rej_ord_lst.append(order)
            elif use == 1:
                use_dif_lst.append(wv_fit-wv_map)
                use_pix_lst.append(pix)
                use_ord_lst.append(order)
                use_wv_lst.append(wv_map)
        else:
            continue
    file.close()

    use_dif_lst = np.array(use_dif_lst)
    use_pix_lst = np.array(use_pix_lst)
    use_ord_lst = np.array(use_ord_lst)
    use_wv_lst  = np.array(use_wv_lst)

    rej_dif_lst = np.array(rej_dif_lst)
    rej_pix_lst = np.array(rej_pix_lst)
    rej_ord_lst = np.array(rej_ord_lst)

    wvstd = use_dif_lst.std()

    rv_dif_lst = use_dif_lst/use_wv_lst*299792458.
    rvstd = rv_dif_lst.std()*1e-3

    nlines = use_dif_lst.size+rej_dif_lst.size
    nused = use_dif_lst.size
    print(file,use_ord_lst,rej_ord_lst)
    order_from = min(use_ord_lst.min(), rej_ord_lst.min())
    order_to   = max(use_ord_lst.max(), rej_ord_lst.max())
    norders = order_to - order_from + 1

    wv_min = use_wv_lst.min()
    wv_max = use_wv_lst.max()

    return wvstd,rvstd,nused,nlines,norders,wv_min,wv_max


def load_reg(filename):
    xnodes = np.array([])
    ynodes = np.array([])
    file = open(filename)
    for row in file:
        row = row.strip()
        if row[0:5]=='point':
            row = row.split('#')[0]
            row = row.replace('point','')
            row = row.replace('(','')
            row = row.replace(')','')
            g = row.split(',')
            x = float(g[0]) - 1
            y = float(g[1]) - 1
            i = np.searchsorted(xnodes,x)
            xnodes = np.insert(xnodes,i,x)
            ynodes = np.insert(ynodes,i,y)
    file.close()
    return (xnodes,ynodes)

def mosaic_flat(flatname):

    conf = load_conf('HRS.conf')
    flat_exptime = conf['flat_exptime']

    # prepare data list
    data_lst = {}
    for row in flat_exptime:
        color = row[0]
        data = pf.getdata('flat_%s.fits'%color)
        yrows,xrows = data.shape
        xcenter = int(xrows/2)
        flat1d = data[:,xcenter]
        data_lst[color] = flat1d

    # load boundaries from *.reg files
    bound_lst = np.zeros(xrows,dtype=np.int32)
    for fname in os.listdir('./'):
        if fname[-4:]=='.reg':
            xnodes,ynodes = load_reg(fname)
            print(fname,xnodes,ynodes)
            f = intp.InterpolatedUnivariateSpline(xnodes,ynodes,k=3)
            bound = np.int32(f(np.arange(xrows)))
            bound_lst = np.vstack((bound_lst,bound))
    bound_lst.sort(axis=0)

    # check number of flats and boundaries
    if len(bound_lst)==0:
        print('Cannot find any .reg file')
        raise ValueError
    if len(flat_exptime) != len(bound_lst):
        print('number of flats = %d, while number of boundaries = %d'%(
                len(flat_exptime), len(bound_lst)))
        raise ValueError

    # prepare central nodes
    nodes = bound_lst[:,xcenter]

    # prepare select_area
    select_area = {}
    for row in flat_exptime:
        color = row[0]
        select_area[color] = np.zeros(len(nodes))>0

    nflats  = len(flat_exptime)
    nbounds = nodes.size

    fig = plt.figure(figsize=(12,7))

    for i in range(nflats):
        ax = fig.add_subplot(nflats+1,1,i+1)
        ax.flat_color = flat_exptime[i][0]

    # white flat is the final flat
    ax = fig.add_subplot(nflats+1,1,nflats+1)
    ax.flat_color = 'white'

    def replot():
        # process the color axes
        for k in range(nflats):
            ax = fig.get_axes()[k]
            x1,x2 = ax.get_xlim()
            y1,y2 = ax.get_ylim()
            ax.cla()
            color = ax.flat_color
            flat1d = data_lst[color]
            ax.plot(flat1d,color=color,ls='-',alpha=0.3)
            for x in nodes:
                ax.plot([x,x],[0,65535],'k--')
            for i in range(nbounds):
                xfrom = nodes[i]
                if i == nbounds-1:
                    xto = yrows
                else:
                    xto = nodes[i+1]
                if select_area[color][i]:
                    ax.plot(np.arange(xfrom,xto),flat1d[xfrom:xto],color=color,ls='-')
            ax.set_xlim(x1,x2)
            ax.set_ylim(y1,y2)

        # process the white flat
        ax = fig.get_axes()[nflats]
        x1,x2 = ax.get_xlim()
        y1,y2 = ax.get_ylim()
        ax.cla()
        for row in flat_exptime:
            color = row[0]
            flat1d = data_lst[color]
            for x in nodes:
                ax.plot([x,x],[0,65535],'k--')
            for i in range(nbounds):
                xfrom = nodes[i]
                if i == nodes.size-1:
                    xto = yrows
                else:
                    xto = nodes[i+1]
                if select_area[color][i]:
                    ax.plot(np.arange(xfrom,xto),flat1d[xfrom:xto],color='k',ls='-')
            ax.set_xlim(x1,x2)
            ax.set_ylim(y1,y2)
        ax.set_xlabel('Spatial Direction',family='serif')

        # process the ticks, y-labels
        for ax in fig.get_axes():
            ax.set_ylabel('ADU',family='serif')
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_family('serif')
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_family('serif')

        fig.canvas.draw()

    def onclick(event):
        ax = event.inaxes
        if ax != None:
            color = ax.flat_color
            if color != 'white':
                i = np.searchsorted(nodes, event.xdata) - 1
                select_area[color][i] = np.logical_not(select_area[color][i])
                if select_area[color][i]:
                    for row in flat_exptime:
                        if row[0] != color:
                            select_area[row[0]][i] = False
                replot()

    # first drawing
    for ax in fig.get_axes():
        ax.set_xlim(0,xrows-1)
        ax.set_ylim(0,65535)
        replot()

    fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show(block=False)

    _ = input('Press [Enter] to continue ')

    # check white flat
    all_area = np.zeros(nbounds,dtype=np.int32)
    for color in select_area:
        all_area += select_area[color]
    if all_area.all() != np.ones(nbounds,dtype=np.int32).all():
        print('The white flat is not completed')
        raise ValueError

    # get white flat
    flat = np.zeros((yrows,xrows))
    for i in range(nbounds):
        yfrom = bound_lst[i]
        if i == nbounds-1:
            yto = np.zeros(xrows,dtype=np.int32) + yrows
        else:
            yto = bound_lst[i+1]
        for color in select_area:
            if select_area[color][i]:
                y,x = np.mgrid[:yrows,:xrows]
                m1 = y>=yfrom
                m2 = y<yto
                m = np.logical_and(m1,m2)
                colorflat = pf.getdata('flat_%s.fits'%color)
                flat += m*colorflat

    # save white fits
    if os.path.exists(flatname):
        os.remove(flatname)
    pf.writeto(flatname,flat)

def view(filename):
    spec = load_multispec(filename)
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_axes([0.1,0.15,0.85,0.75])

    def plot_order(order):
        ax.cla()
        ax.currentorder = order
        spd = spec[order]

        wvc = (spd.wv_from+spd.wv_to)/2.0
        color = get_rainbow_color(wvc)
        ax.plot(spd.wv, spd.flux,'-',color=np.array(color)/255.)
        ax.set_xlabel(u'Wavelength (\xc5)')
        ax.set_ylabel('Flux (ADU)')
        ax.set_title('%s - Order %d'%(filename,order))
        if (spd.flux < 0).sum()==0:
            ax.set_ylim(0,)
        ax.set_xlim(spd.wv_from, spd.wv_to)
        ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g'))
        fig.canvas.draw()

    def on_key(event):
        if event.key == 'up':
            if spec.orders.index(ax.currentorder)!=0:
                plot_order(ax.currentorder + 1)
        elif event.key == 'down':
            if spec.orders.index(ax.currentorder)!=len(spec.orders)-1:
                plot_order(ax.currentorder - 1)

    order = spec.orders[int(len(spec.orders)/2.)]
    plot_order(order)

    fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()

def convert(filetype):
    for fname in os.listdir('./'):
        if fname[-9:] in ['_1ds.fits','_blz.fits']:
            spec = load_multispec(fname)
            for o in spec.orders:
                newname = '%s_o%03d.%s'%('.'.join(fname.split('.')[0:-1]),o,filetype)
                if filetype == 'txt':
                    spec[o].save_txt(newname)
                elif filetype == 'fits':
                    spec[o].save_fits(newname)
                print('%s order %03d -> %s'%(fname, o, newname))


def convert2(filetype):
    '''deprecated convert function'''
    logfile = find_logfile()
    log = read_log(logfile)

    for item in log.item_list:
        if item.skip:
            continue
        filename = item.get_filename(log.filename_composition)
        for key in ['1ds','blz']:
            fname = '%s_%s.fits'%('.'.join(filename.split('.')[0:-1]),key)
            if os.path.exists(fname):
                spec = load_multispec(fname)
                for o in spec.orders:
                    newname = '%s_o%03d.%s'%('.'.join(fname.split('.')[0:-1]),o,filetype)
                    if filetype == 'txt':
                        spec[o].save_txt(newname)
                    elif filetype == 'fits':
                        spec[o].save_fits(newname)
                    print('%s order %03d -> %s'%(fname, o, newname))

def clean():
    logfile = find_logfile()
    log = read_log(logfile)

    for item in log.item_list:
        if item.skip:
            continue
        filename = item.get_filename(log.filename_composition)
        for key in ['ovr','flt','bkg','sum']:
            fname = '%s_%s.fits'%('.'.join(filename.split('.')[0:-1]),key)
            if os.path.exists(fname):
                os.remove(fname)


def load_conf(filename):
    conf = {}
    file = open(filename)
    for row in file:
        row = row.strip()
        if len(row)==0 or row[0]=='#':
            continue
        g = row.split('=')
        key   = g[0].strip()
        value = g[1].split('#')[0].strip()
        if key[-12:]=='flat exptime':
            color = key.split()[0]
            try:
                exptime = int(value)
            except:
                exptime = float(value)
            if 'flat_exptime' not in conf:
                conf['flat_exptime'] = []
            conf['flat_exptime'].append([color,exptime])
        elif key=='ccd map':
            conf['ccd_map'] = value
        elif key=='apsize':
            conf['apsize'] = int(value)
    file.close
    return conf

def help():
    #print help
    print('''Echelle Spectra Reduction Software for Xinglong 2.16m HRS
    author: Liang Wang
    email: wang.leon@gmail.com
    last modified: 2013-10-08
    ''')
    print('''%s check
    check the exptime, start time, missing files.'''%__file__)
    print( )
    print('''%s flat
    plot a section of flat frames in a matplotlib window'''%__file__)
    print( )
    print('''%s overscan
    overscan correction'''%__file__)
    print( )
    print('''%s reduce
    extract the 1-D spectra'''%__file__)
    print( )
    print('''%s calib
    re-calibrate the wavelength'''%__file__)
    print( )
    print('''%s update-header
    copy log information to headers of _1ds.fits and _blz.fits'''%__file__)
    print( )
    print('''%s clean
    clean the directory after 1-D extraction'''%__file__)
    print( )
    print('''%s view FILENAME
    review the 1-D spectra in a matplotlib window'''%__file__)
    print( )
    print('''%s convert [txt/fits]
    convert the 1-D spectra to ascii or fits files'''%__file__)
    print( )

if __name__=='__main__':

    try:
        argv = sys.argv[1]
    except:
        help()
        exit()

    if argv=='clean':
        clean()

    elif argv=='check':
        # check fits files
        check_data()

    elif argv=='flat':
        # check flats
        plot_flat()

    elif argv=='overscan':
        # overscan correction
        overscan()

    elif argv=='mosaic':
        # mosaic flats
        mosaic_flat('flat.fits')

    elif argv=='reduce':
        # reduce 1-D spectra
        extract1d()

    elif argv=='calib':
        extract1d(calib=True)

    elif argv=='view':
        # review the extracted 1-D spectra in a matolotlib window
        if len(sys.argv)>=3:
            filename = sys.argv[2]
            if os.path.exists(filename):
                if filename[-9:] in ['_1ds.fits','_blz.fits']:
                    view(filename)
                else:
                    print('we do not support this type of fits'%filename)
            else:
                print('%s does not exist'%filename)
        else:
            print('please give the filename')

    elif argv=='convert':
        # convert the MULTISPEC fits to txt of fits
        if len(sys.argv)>=3:
            filetype = sys.argv[2]
            if filetype in ['txt','fits']:
                convert(filetype)
            else:
                print('we do not support this type')
        else:
            print('please give the file type (txt/fits)')

    elif argv=='combine':
        # combine several spectra
        args = sys.argv[2:]
        try:
            i = args.index('-o')
        except:
            print('Please specify the output name with -o [OUTNAME]')
            exit()
        outname = args[i+1]
        if os.path.exists(outname):
            ans = input('%s will be overidden. Confirm ? [y/N]: '%outname)
            if ans.strip().lower() in ['no','n','']:
                exit()
        input_files = args[0:i]

        for fname in input_files:
            if not os.path.exists(fname):
                print('Input file %s does not exist'%fname)
                exit()

        combine(input_files,outname)

    elif argv=='update-header':
        # update fits header of _1ds.fits and _blz.fits
        update_header()


    else:
        help()
