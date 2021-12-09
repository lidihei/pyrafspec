'''This file is to get lists of flat, lamp, bias and targets from 216 cm log file
'''

import numpy as np
import pandas as pd
import os
from pyrafspec import default_logheader

def logs2list_2021(fname, dire=None, equipment='G10_E9',flat_expt = 900, lamp_expt= 300 ):
    '''convert log file of BFOSC (observed after 2021) to the lists of flat, bias, arc lamp and objecti
    parameters:
    ----------------
    fname: [str] the name of the log file  
    dire: [str] the directory of data
    equipment: [str]
    flat_expt: [float] the exposure time of flat
    lamp_expt: [float] the exposure time of lamp
    returns:
    ----------------
    bias_list: [list]
    flat_list: [list]
    star_list: [list] 
    lamp_list: [list]
    '''
    logs = pd.read_table(fname, header=0, delim_whitespace=True)
    flist = logs['FileName']
    expt = logs['ExpT']
    flat_list = []
    bias_list = []
    star_list = []
    lamp_list = []
    for i, fname in enumerate(flist):
        if ('SPECSFLAT' in fname) and (equipment in fname) and (expt[i] == flat_expt):
           flat_list.append(os.path.join(dire, f'{fname}.fit'))
        if ('SPECSLAMP' in fname) and (equipment in fname) and (expt[i] == lamp_expt):
           lamp_list.append(os.path.join(dire, f'{fname}.fit'))
        if ('SPECSTARGET' in fname) and (equipment in fname):
           star_list.append(os.path.join(dire, f'{fname}.fit'))
        if ('BIAS' in fname) and (equipment in fname):
           bias_list.append(os.path.join(dire, f'{fname}.fit'))
    return bias_list, flat_list, star_list, lamp_list


def getfnamelist(fname, prefixi=8):
    '''
    fname: [stri] e.g. 20201216001-3
    returns:
    ------------
    fnamelst: [lst] e.g. ['20201216-0001', '20201216-0002', '20201216-0003']
    '''
    fnamelst = []
    _ids = fname[prefixi:].strip()
    ids = _ids.split('-')
    prefix = fname[0:prefixi]
    if len(ids) == 1:
        _id = int(_ids)
        fnamelst.append(prefix+f'-{_id:04d}')
    else:
        for _id in np.arange(int(ids[0]), int(ids[1])+1):
            _fname = prefix+f'-{_id:04d}'
            fnamelst.append(_fname)
    return fnamelst


def logs2list_2020(fname, dire, equipment='1.6+G10+E9',flat_expt = 900, lamp_expt= 300, lamp='Fe/Ar', prefixi=8, filetype='.fit', skiprows=11):
    '''convert log file of BFOSC (observed after 2020) to the lists of flat, bias, arc lamp and objecti
    parameters:
    ----------------
    fname: [str] the name of the log file  
    dire: [str] the directory of data
    equipment: [str]
    flat_expt: [float] the exposure time of flat
    lamp_expt: [float] the exposure time of lamp
    returns:
    ----------------
    bias_list: [list]
    flat_list: [list]
    star_list: [list] 
    lamp_list: [list]
    '''
    filename = fname
    header = ['FileName', 'obj', 'Btime', 'ExpT', 'Ra', 'Dec', 'epoch', 'notes']
    logs = pd.read_table(filename, names=header, skiprows=skiprows,delim_whitespace=True)
    flist = logs['FileName']
    expt = logs['ExpT']
    obj = logs['obj']
    notes = logs['notes']
    flat_list = []
    bias_list = []
    star_list = []
    lamp_list = []
    for i, fname in enumerate(flist):
        expti = np.float(expt[i])
        if (obj[i].lower() == 'bias') and (expti == 0):
           fnamelst = getfnamelist(fname, prefixi=prefixi)
           bias_list += fnamelst
        fnamelst = getfnamelist(fname, prefixi=prefixi)
        if obj[i].lower() == 'flat':
           if (expti == flat_expt):
              flat_list += fnamelst
           else: continue
        if notes[i] != equipment:
           continue
        elif obj[i] == lamp:
           if expti == lamp_expt:
              lamp_list += fnamelst
           else: continue
        else:
           star_list += fnamelst
    bias_list = [os.path.join(dire, _i+filetype) for _i in bias_list]
    flat_list = [os.path.join(dire, _i+filetype) for _i in flat_list]
    star_list = [os.path.join(dire, _i+filetype) for _i in star_list]
    lamp_list = [os.path.join(dire, _i+filetype) for _i in lamp_list]
    return bias_list, flat_list, star_list, lamp_list


def match_star2lamp(fname, equipment='G10_E9', lamp_expt= 300, fout=None, **kwargs):
    ''' match star with lamp for 216 logs produced after 2021
    parameters:
    ----------------
    fname: [str] the name of the log file
    lamp_expt: [float] the exposure time of lamp
    returns:
    -----------------
    star_list: [list]
    star_lamps: [list] the list of lamp cooresponding to the star
    '''
    logs = pd.read_table(fname, header=0, delim_whitespace=True, **kwargs)
    flist = logs['FileName']
    expt = logs['ExpT']
    ra = logs['Ra']
    dec = logs['Dec']
    btime = logs['UTCTime']
    star_list = []
    star_ra = []
    star_dec = []
    lamp_list = []
    lamp_ra = []
    lamp_dec = []
    lamp_time = []
    star_time = []
    times = []
    for _i, _btime in enumerate(btime):
        hh, mm, ss = _btime.split(':')
        times.append(np.float32(hh) +  np.float32(mm)/60 + np.float32(ss)/3600)
    for i, fname in enumerate(flist):
        if ('SPECSLAMP' in fname) and (equipment in fname) and (expt[i] == lamp_expt):
           lamp_list.append(fname)
           lamp_ra.append(ra[i])
           lamp_dec.append(dec[i])
           lamp_time.append(times[i])
        if ('SPECSTARGET' in fname) and (equipment in fname):
          star_list.append(fname)
          star_ra.append(ra[i])
          star_dec.append(dec[i])
          star_time.append(times[i])
    lamp_list = np.array(lamp_list)
    lamp_ra = np.array(lamp_ra)
    lamp_dec = np.array(lamp_dec)
    lamp_time = np.array(lamp_time)
    star_list = np.array(star_list)
    star_ra = np.array(star_ra)
    star_dec = np.array(star_dec)
    star_time = np.array(star_time)
    star_lamps = []
    strs = 'starlist,lamplist\n'
    for i, star in enumerate(star_list):
        try:
           rai = star_ra[i]
           deci = star_dec[i]
           absdtime = np.abs(lamp_time-star_time[i])
           ind = (lamp_ra ==star_ra[i]) & (lamp_dec == star_dec[i])
           _indtime = np.argmin(absdtime[ind])
           starlamp = lamp_list[ind][_indtime]
        except:
           print(f'witout corresponding lamp file: {star}')
           starlamp = 'witout_corresponding_lamp'
        star_lamps.append(starlamp)
        strs = strs+ f'{star},{starlamp}\n'
    if fout is not None:
       ifile = open(fout, 'w')
       ifile.writelines(strs)
       ifile.close()
    return star_list, np.array(star_lamps)


def match_star2lamp_2020(fname, equipment='1.6+G10+E9', lamp_expt= 300, lamp='Fe/Ar', prefixi=8, fout=None):
    ''' match star with lamp for 216 logs produced before 2021
    parameters:
    ----------------
    fname: [str] the name of the log file  
    dire: [str] the directory of data
    equipment: [str]
    flat_expt: [float] the exposure time of flat
    lamp_expt: [float] the exposure time of lamp
    returns:
    ----------------
    star_list: [list]
    star_lamps: [list] the list of lamp cooresponding to the star
    '''
    header = ['FileName', 'obj', 'Btime', 'ExpT', 'Ra', 'Dec', 'epoch', 'notes']
    logs = pd.read_table(fname, names=header, skiprows=11,delim_whitespace=True)
    flist = logs['FileName']
    expt = logs['ExpT']
    obj = logs['obj']
    notes = logs['notes']
    btime = logs['Btime']
    ra = logs['Ra']
    dec = logs['Dec']
    times = []
    for _i, _btime in enumerate(btime):
        hh, mm, ss = _btime.split(':')
        times.append(np.float32(hh) +  np.float32(mm)/60 + np.float32(ss)/3600)
    lamp_list = []
    star_list = []
    lamp_ra = []
    lamp_dec = []
    star_ra = []
    star_dec = []
    lamp_time = []
    star_time = []
    for i, fname in enumerate(flist):
        expti = np.float(expt[i])
        if (obj[i].lower() == 'bias'):
           continue
        if obj[i].lower() == 'flat':
           continue
        fnamelst = getfnamelist(fname, prefixi=prefixi)
        nn = len(fnamelst) 
        if obj[i] == lamp:
           if expti == lamp_expt:
              lamp_list += fnamelst
              lamp_ra += [ra[i]]*nn
              lamp_dec += [dec[i]]*nn
              lamp_time += [times[i]]*nn
           else: continue
        else:
           star_list += fnamelst
           star_ra += [ra[i]]*nn
           star_dec += [dec[i]]*nn
           star_time += [times[i]]*nn
    lamp_list = np.array(lamp_list)
    lamp_ra = np.array(lamp_ra)
    lamp_dec = np.array(lamp_dec)
    lamp_time = np.array(lamp_time)
    star_list = np.array(star_list)
    star_ra = np.array(star_ra)
    star_dec = np.array(star_dec)
    star_time = np.array(star_time)
    star_lamps = []
    strs = 'starlist,lamplist,ra,dec\n'
    for i, star in enumerate(star_list):
        try:
           rai = star_ra[i]
           deci = star_dec[i]
           absdtime = np.abs(lamp_time-star_time[i])
           ind = (lamp_ra ==star_ra[i]) & (lamp_dec == star_dec[i])
           _indtime = np.argmin(absdtime[ind])
           starlamp = lamp_list[ind][_indtime]
        except:
           print(f'witout corresponding lamp file: {star}')
           starlamp = 'witout_corresponding_lamp'
        star_lamps.append(starlamp)
        strs = strs+ f'{star},{starlamp},{rai},{deci}\n'
    if fout is not None:
       ifile = open(fout, 'w')
       ifile.writelines(strs)
       ifile.close()
    return star_list, np.array(star_lamps)


class convertobslog:

      def __init__(self, dire, logname):
          '''
          dire: [str] the directory stored fit and log files
          logname: [str] the name of log file
          '''
          self.dire = dire
          self.logname = logname

      def convert2obslog(self, date=None, program=None, observer=None, operator=None,
                               observatory=None, telescope=None, instrument= None, detector=None,
                               fiberdiameter =None, slitwidth=None, T_I2=None,
                               readout=None, gain=None, timesystem=None, timezone=None,
                               filenamecomposition =None, cols=None, fout=None, obslogheader=None, equipment='G10_E9'):
          '''
          rewrite the log file of Xinglong 216 cm to the *.obslog of pyrafspec
          paramters:
          -------------
          date: [str] e.g. '2021-02-07'
          program: [str] e.g. 'Star Search'
          observer: [str] e.g. 'J. Li'
          operator: [str] e.g. 'J. Li'
          observatory: [str] e.g. 'NAOC-Xinglong'
          telescope: [str] e.g. '2.16m'
          instrument: [str] e.g. 'BFOSC'
          detector: [str] e.g. 'E2V CCD203-82'
          fiberdiameter:[str] e.g. '2"'
          slitwidth: [str] e.g. '1.6"'
          T_I2: [str] I2 Temperature, e.g. '60'
          readout: [str] e.g. 'Left Top & Bottom'
          gain: [str] e.g. '1.0 e-/ADU'
          timesystem: [str] e.g. 'UTC'
          timezone: [str] e.g. '+08:00'
          filenamecomposition: [str] e.g 'dire/%fname.fit'
          cols: [str] e.g. 'id,object,time,exptime,note'
          fout: [str] the output file name
          obslogheader: [str]
          equipment: [str], e.g. 'G10_E9'
          '''
          logname = self.logname
          dire = self.dire
          if obslogheader is None: obslogheader =default_logheader.obslogheader
          if date is not None: obslogheader=obslogheader.replace('<date>', date)
          if program is not None: obslogheader=obslogheader.replace('<program>', program)
          if observer is not None: obslogheader=obslogheader.replace('<observer>', observer)
          if operator is not None: obslogheader=obslogheader.replace('<operator>', operator)
          if observatory is not None: obslogheader=obslogheader.replace('<observatory>', observatory)
          if telescope is not None: obslogheader=obslogheader.replace('<telescope>', telescope)
          if instrument is not None: obslogheader=obslogheader.replace('<instrument>', instrument)
          if detector is not None: obslogheader=obslogheader.replace('<detector>', detector)
          if fiberdiameter is not None: obslogheader=obslogheader.replace('<fiberdiameter>', fiberdiameter)
          if slitwidth is not None: obslogheader=obslogheader.replace('<slitwidth>', slitwidth)
          if T_I2 is not None: obslogheader=obslogheader.replace('<T_I2>', T_I2)
          if readout is not None: obslogheader=obslogheader.replace('<readout>', readout)
          if gain is not None: obslogheader=obslogheader.replace('<gain>', gain)
          if timesystem is not None: obslogheader=obslogheader.replace('<timesystem>', timesystem)
          if timezone is not None: obslogheader=obslogheader.replace('<timezone>', timezone)
          if filenamecomposition is None:
             filenamecomposition = f'{dire}/%fname.fit'
          if os.path.dirname(logname) == '':
             logname = os.path.join(dire, logname)
          if fout is None:
             fout = logname.replace('.log', '.obslog')
 
          obslogheader=obslogheader.replace('<filenamecomposition>', filenamecomposition)
          logs = pd.read_table(logname, header=0, delim_whitespace=True)
          if cols is not None:
             obslogheader=obslogheader.replace('<cols>', cols)
          else:
            cols = ','.join(['id', 'object']+list(logs.keys()))
            cols = cols.replace('UTCTime', 'time')
            cols = cols.replace('ExpT', 'exptime')
            cols = cols.replace('FileName', 'fname')
            obslogheader=obslogheader.replace('<cols>', cols)
          FileName = logs['FileName']
          UTCTime = logs['UTCTime']
          ObservorName = logs['ObservorName']
          Ra = logs['Ra']
          Dec = logs['Dec']
          Epoch = logs['Epoch']
          ExpT = logs['ExpT']
          DataType = logs['DataType']
          for i, fname in enumerate(FileName):
              obj = None
              if ('SPECSFLAT' in fname) and (equipment in fname):
                 obj = 'flat'
              if ('SPECSLAMP' in fname) and (equipment in fname):
                 if 'FeAr' in fname:
                    obj = 'fear'
                 elif 'ThAr' in fname:
                    obj = 'thar'
                 else: obj = lamp
              if ('SPECSTARGET' in fname) and (equipment in fname):
                 obj = fname.split('_')[2]
              if ('BIAS' in fname) and (equipment in fname):
                 obj = 'bias'
              if obj is not None:
                 stri = f'|{i}|{obj}|{FileName[i]}|{UTCTime[i]}|{ObservorName[i]}|{Ra[i]}|{Dec[i]}|{Epoch[i]}|{ExpT[i]}|{DataType[i]}|\n'
                 obslogheader += stri
          ofile = open(fout, 'w')
          ofile.writelines(obslogheader)
        
