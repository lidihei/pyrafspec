# list files and handle fits file
import os, sys
import numpy as np
from astropy.io import fits
from .obslog import *
import scipy.interpolate as intp
import scipy.optimize    as opt
from .mathfunc import *


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

        head = fits.Header()
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
        fits.writeto(filename,fluxnew,head)

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

        data,head = fits.getdata(filename,header=True)
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

        data,head = fits.getdata(filename,header=True)
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
        elif row[-4:]=='.fit':
            file.write('%s_%s.fits%s'%(row[0:-4],surfix,os.linesep))
    fileref.close()
    file.close()

def copy_lstfile(reflst, filename, direcp = './'):
    '''copy files of *.lst and modify the lst'''
    file = open(filename,'w')
    fileref = open(reflst)
    for row in fileref:
        row = row.strip()
        if row[-5:]=='.fits':
           basename = os.path.basename(row)
           file.write('%s%s'%(basename,os.linesep))
           os.system(f'cp {row} {direcp}')
    fileref.close()
    file.close()


def update_header(logfile=None):
    if logfile is None:
       logfile = find_logfile()
    log = read_log(logfile)

    for item in log.get_itemlist():
        if item.skip:
            continue
        filename = item.get_filename(log.filename_composition)
        for key in ['1ds','blz']:
            fname = '%s_%s.fits'%('.'.join(filename.split('.')[0:-1]),key)
            if os.path.exists(fname):
                f = fits.open(fname,mode='update')
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
                    print(f'fname = {fname}')
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

