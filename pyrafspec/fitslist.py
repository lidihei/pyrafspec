# list files and handle fits file
import os, sys
import numpy as np
from astropy.io import fits
from .obslog import *

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

