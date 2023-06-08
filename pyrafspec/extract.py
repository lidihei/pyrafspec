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

from .obslog import *
from .overscan import *
from .mathfunc import *
from .fitslist import *
from .config import *
from .mosaicflat import *
from .tools import *
from glob import glob


def extract1d(calib=False, apedit=True, apsum=True, logfile=None, lamp='ThAr', lamp_exptime=None,
              starlstname_wv=None, lamplstname_wv=None, refspsort='DATE-STA',\
             lampreference=None):
    '''
    starlstname: [str] e.g. star.lst
                 dire/ljg2m401-yf01-20210104-0182-e00.fits
                 dire/ljg2m401-yf01-20210104-0253-e00.fits
    lampreference: [str] e.g. lamp-ref.fits, which should locate in the working directory, and eclamp-ref be in database
    '''
    from pyraf import iraf

    # remove the old iraf.log file
    if not calib and os.path.exists('iraf.log'):
        os.remove('iraf.log')

    if not calib:
        os.system('mkiraf')
    if (logfile is None) or (logfile == ''):
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
    direname = os.path.dirname(log.filename_composition)
    log.save_file(f'{direname}/lamp.lst', object=lamp, exptime=lamp_exptime)
    log.save_file(f'{direname}/star.lst', object='Star')
    if has_iodn:
        log.save_file(f'{direname}/iodn.lst', object='Iodn')

    # make flat list of different exposure times
    conf = load_conf('config/instrument.conf')
    flat_exptime = conf['flat_exptime']
    apsize = conf['apsize']
    print(f'lijiao flat_exptime = {flat_exptime}')
    flatfitsname = os.path.join(direname, 'flat.fits')
    for row in flat_exptime:
        print(f'lijiao raw = {row[0]}')
        color   = row[0]
        exptime = row[1]
        flat_lst = log.get_filenamelist(
                   object  = 'Flat',
                   exptime = exptime
                   )

        # write filenames to list
        flat_lst_file = open(f'{direname}/flat_%s.lst'%color,'w')
        flat_lst_file.write(os.linesep.join(flat_lst))
        flat_lst_file.close()

        prepare_lst(f'{direname}/flat_{color}.lst','ovr',f'{direname}/flat_{color}_ovr.lst')

        if calib and os.path.exists(f'{direname}/flat_{color}.lst'):
            continue

        # Combine flat of different colors
        print( )
        print('Combine %d %s Flats (EXPTIME = %s seconds)'%(len(flat_lst),color.title(),str(float(exptime))))
        print( )
        iraf.imcomb.unlearn()
        iraf.imcomb.combine = 'median'
        iraf.imcomb.logfile = 'iraf.log'
        delete_fits(f'{direname}/flat_{color}.fits')
        print(f'lijiao {color}')
        flatcoloris1 = (len(flat_exptime) == 1) and (color=='flat')
        if flatcoloris1:
           iraf.imcomb(f'@{direname}/flat_{color}_ovr.lst', flatfitsname)
        else:
           iraf.imcomb(f'@{direname}/flat_{color}_ovr.lst', f'{direname}/flat_{color}.fits')
    
    if calib and os.path.exists(flatfitsname):
        pass
    else:
        # mosaic flats
        print( )
        print('Mosaic Flats')
        print( )
        if not flatcoloris1:
           mosaic_flat(flatfitsname)

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
    if not calib and apedit:
        #os.system(f'cp {direname}/flat.fits ./')
        iraf.apedit(flatfitsname)
        #os.system('rm flat.fits')

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
    iraf.aptrace.niterate   = 7
    if not calib and apedit:
        #os.system(f'cp {direname}/flat.fits ./')
        iraf.aptrace(flatfitsname)
        #os.system('rm flat.fits')

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
    flat_nfitsname = os.path.join(direname, 'flat_n.fits')
    if not calib and apedit:
        delete_fits(flat_nfitsname)
        #os.system(f'cp {direname}/flat.fits ./')
        iraf.apflatten(flatfitsname,flat_nfitsname)
        #os.system('rm flat.fits')

    iraf.imarith.unlearn()

    # parse lamp
    print('----------------------------parse lamp-----------------------------')
    lamplstname = os.path.join(direname, 'lamp.lst')
    lamp_ovrlstname = os.path.join(direname, 'lamp_ovr.lst')
    lamp_fltlstname = os.path.join(direname,'lamp_flt.lst')
    prepare_lst(lamplstname, 'ovr', lamp_ovrlstname)
    prepare_lst(lamplstname, 'flt', lamp_fltlstname)
    if not calib and apedit:
        delete_fits(lamp_fltlstname)
        iraf.imarith(f'@{lamp_ovrlstname}','/',flat_nfitsname,f'@{lamp_fltlstname}')

    # parse star
    print('----------------------------parse star-----------------------------')
    starlstname = os.path.join(direname, 'star.lst')
    star_ovrlstname = os.path.join(direname, 'star_ovr.lst')
    star_fltlstname = os.path.join(direname,'star_flt.lst')
    prepare_lst(starlstname, 'ovr', star_ovrlstname)
    prepare_lst(starlstname, 'flt', star_fltlstname)
    if not calib and apedit:
        delete_fits(star_fltlstname)
        iraf.imarith(f'@{star_ovrlstname}','/',flat_nfitsname, f'@{star_fltlstname}')

    # parse iodine
    print('----------------------------parse iodine-----------------------------')
    if has_iodn:
        iodnlstname = os.path.join(direname, 'iodn.lst')
        iodn_ovrlstname = os.path.join(direname, 'iodn_ovr.lst')
        iodn_fltlstname = os.path.join(direname, 'iodn_flt.lst')
        prepare_lst(iodnlstname, 'ovr', iodn_ovrlstname)
        prepare_lst(iodnlstname, 'flt', iodn_fltlstname)
        if not calib and apedit:
            delete_fits(iodn_fltlstname)
            iraf.imarith(f'@{iodn_ovrlstname}','/',flat_nfitsname,f'@{iodn_fltlstname}')

    #os.system(f'cp {direname}/flat.fits ./')
    print( )
    print('Background Correction')
    print( )
    iraf.apresize.unlearn()
    iraf.apresize.llimit = -apsize
    iraf.apresize.ulimit =  apsize
    iraf.apresize.ylevel = 'INDEF'
    iraf.apresize.bkg    = 'no'

    iraf.apscatter.unlearn()
    iraf.apscatter.referen = flatfitsname
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
    star_bkglstname =  os.path.join(direname,'star_bkg.lst')
    prepare_lst(starlstname, 'bkg', star_bkglstname)
    if not calib and apedit:
        delete_fits(star_bkglstname)
        iraf.apscatter(f'@{star_fltlstname}', f'@{star_bkglstname}')

    if has_iodn:
        iodn_bkglstname =  os.path.join(direname,'iodn_bkg.lst')
        prepare_lst(iodnlstname, 'bkg', iodn_bkglstname)
        if not calib and apedit:
            delete_fits(iodn_bkglstname)
            iraf.apscatter(f'@{iodn_fltlstname}', f'@{iodn_bkglstname}')

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

    star_sumlstname =  os.path.join(direname,'star_sum.lst')
    prepare_lst(starlstname,'sum',star_sumlstname)
    iraf.apsum.unlearn()
    iraf.apsum.output  = f'@{star_sumlstname}'
    iraf.apsum.format  = 'echelle'
    iraf.apsum.referen = flatfitsname
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
    if not calib or apsum:
        delete_fits(star_sumlstname)
        iraf.apsum(f'@{star_bkglstname}')

    lamplstname = os.path.join(direname, 'lamp.lst')
    lamp_sumlstname = os.path.join(direname, 'lamp_sum.lst')
    prepare_lst(lamplstname,'sum',lamp_sumlstname)
    iraf.apsum.output  = f'@{lamp_sumlstname}'
    if not calib and apedit:
        delete_fits(lamp_sumlstname)
        iraf.apsum(f'@{lamp_fltlstname}')

    if has_iodn:
        prepare_lst(iodnlstname,'sum',iodn_sumlstname)
        iodn_sumlstname = os.path.join(direname, 'iodn_sum.lst')
        iraf.apsum.output  = f'@{iodn_sumlstname}'
        if not calib and epedit:
            delete_fits(iodn_sumlstname)
            iraf.apsum(f'@{iodn_bkglstname}')

    print( )
    print('Wavelength Calibration')
    print( )

    prompt = 'Do you want to use an [E]xisting reference or identify a [N]ew lamp image? [E/N]: '
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
        # in lamp_sum.lst
        file = open(lamp_sumlstname)
        for row in file:
            fn = row.strip()
            _direname = os.path.dirname(fn)
            _basename = os.path.basename(fn)
            ecfn = os.path.join('database', 'ec'+_basename.replace('.fits',''))
            if os.path.exists(ecfn):
                os.remove(ecfn)
        file.close()

        # list files and their exptimes in lamp_sum.lst
        ref_lst = {}
        no = 1
        file = open(lamp_sumlstname)
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
        if lampreference is not None:
           ref_lst[no] = lampreference
           print(' [%02d] %s'%(no,lampreference))

        # Now the user have to choose a lamp file to identify wavelengths
        while (True):
            n = input('Select a lamp frame for wavelength identification [1]:')
            try:
                n = int(n)
                if 0 < n <= len(ref_lst):
                    break
            except:
                pass

        ref_lamp = ref_lst[int(n)]
        ref_lamp = os.path.basename(ref_lamp)
        ref_name = ref_lamp.replace('.fits','')

        iraf.ecid.unlearn()
        iraf.ecid.maxfeat = 100
        iraf.ecid.coordli = 'linelists$'+f'{lamp.lower()}.dat'
        iraf.ecid(ref_lamp)
    else:
        prompt = 'Cannot find reference in database.\n'
        prompt += 'Put an `ec` file in database/ and press [Enter] to refresh: '
        while(True):
            ref_lst = {}
            no = 1
            #lst = os.listdir('./database')
            lst = glob(os.path.join('database', 'ec*'))
            lst.sort()
            for fname in lst:
                #if fname[0:2]=='ec':
                ref_lst[no] = os.path.basename(fname)

                # get information of this ec file
                res = get_ec_info(fname)
                wvstd, rvstd    = res[0], res[1]
                nused, nlines   = res[2], res[3]
                norders         = res[4]
                wv_min, wv_max  = res[5], res[6]
                #date = '%s-%s-%s'%(fname[2:6],fname[6:8],fname[8:10])

                print(' [%d] %18s %7.5f %6.3f %4d/%4d %3d (%7.1f - %7.1f)'%(
                        no, ref_lst[no][2:], 
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

    
    copy_lstfile(lamp_sumlstname, 'lamp_sum_tmp.lst', direcp = './')
    if lampreference is not None:
       ifile = open('lamp_sum_tmp.lst', 'a')
       ifile.writelines(lampreference)
       ifile.close()
    #copy_lstfile(star_sumlstname, 'star_sum_tmp.lst', direcp = './') 
    iraf.ecreid.unlearn()
    iraf.ecreid.logfile = 'STDOUT, iraf.log'
    print(f'ref_name = "{ref_name}"')
    #iraf.ecreid('f@{lampsumlstname}', ref_name)
    iraf.ecreid('@lamp_sum_tmp.lst', ref_name)

    _ = input('Press [Enter] to continue: ')

    if starlstname_wv is None:
       star_sumlstname_wv =star_sumlstname
    else:
       star_sumlstname_wv = os.path.join(direname, 'star_sumlstname_wv.lst')
       star_1dslstname = os.path.join(direname, 'star_1dslstname_wv.lst')
       prepare_lst(starlstname_wv, '1ds', star_1dslstname)
       prepare_lst(starlstname_wv, 'sum', star_sumlstname_wv)
    if lamplstname_wv is None:
       lamp_sumlstname_wv =lamp_sumlstname
       lamp_sum_wvtmpname = 'lamp_sum_tmp.lst'
    else:
       lamp_sumlstname_wv = os.path.join(direname, 'lamp_sumlstname_wv.lst')
       prepare_lst(lamplstname_wv, 'sum', lamp_sumlstname_wv)
       lamp_sum_wvtmpname ='lamp_sum_wv_tmp.lst'
       copy_lstfile(lamp_sumlstname_wv, lamp_sum_wvtmpname, direcp = './')
       if lampreference is not None:
          ifile = open(lamp_sum_wvtmpname, 'a')
          ifile.writelines(lampreference)
          ifile.close()

    if not new_ecid:
       ifile = open('lamp_sum_tmp.lst', 'a')
       ifile.writelines(f'{ref_name}.fits\n')
       ifile.close()

    iraf.refsp.unlearn()
    #iraf.refsp.referen = f'@lamp_sum_tmp.lst'
    iraf.refsp.referen = f'@{lamp_sum_wvtmpname}'
    if refspsort is None:
       iraf.refsp.sort    = ''
    else:
       iraf.refsp.sort = refspsort
    if lamplstname_wv is None:
       iraf.refsp.time    = 'yes'
    else: 
      iraf.refsp.time    = 'no'
    iraf.refsp.group   = ''
    iraf.refsp.logfile = 'STDOUT, iraf.log'
    iraf.refsp(f'@{star_sumlstname_wv}')
    iraf.refsp(f'@{lamp_sumlstname_wv}')
    if has_iodn:
        #copy_lstfile(iodn_sumlstname, 'iodn_sum1.lst', direcp = './') 
        iraf.refsp(f'@{iodn_sumlstname}')


    iraf.dispcor.unlearn()
    iraf.dispcor.lineari = 'no'
    iraf.dispcor.flux    = 'no'
    iraf.dispcor.logfile = 'iraf.log'
    star_1dslstname = os.path.join(direname, 'star_1ds.lst')
    prepare_lst(starlstname, '1ds', star_1dslstname)
    delete_fits(star_1dslstname)
    iraf.dispcor(f'@{star_sumlstname_wv}', f'@{star_1dslstname}')

    lamp_1dslstname = os.path.join(direname, 'lamp_1ds.lst')
    prepare_lst(lamplstname, '1ds', lamp_1dslstname)
    delete_fits(lamp_1dslstname)
    iraf.dispcor(f'@{lamp_sumlstname}',f'@{lamp_1dslstname}')

    if has_iodn:
        iodn_1dslstname = os.path.join(direname, 'iodn_1ds.lst')
        prepare_lst(iodnlstname, '1ds', iodn_1dslstname)
        delete_fits(iodn_1dslstname)
        iraf.dispcor(f'@{iodn_sumlstname}', f'@{iodn1dslstname}')
        delete_fits(iodn_sumlstname)


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
        delete_fits(f'{direname}/flat_f.fits')
        delete_fits(f'{direname}/flat_a.fits')
        delete_fits(f'{direname}/flat_as.fits')

        iraf.imarith(f'{direname}/flat.fits','/',f'{direname}/flat_n.fits',f'{direname}/flat_f.fits')
        iraf.apsum.output  = f'{direname}/flat_a.fits'
        iraf.apsum.interac = 'no'
        iraf.apsum.review  = 'no'
        iraf.apsum(f'{direname}/flat_f.fits')

        iraf.continuum.unlearn()
        iraf.continuum.ask         = 'YES'
        iraf.continuum.type        = 'fit'
        iraf.continuum.wavescale   = 'no'
        iraf.continuum.logfiles    = 'iraf.log'
        iraf.continuum.interactive = 'no'
        iraf.continuum.order       = 6
        iraf.continuum.low_reject  = 3.0
        iraf.continuum.high_reject = 3.0
        iraf.continuum(f'{direname}/flat_a.fits',f'{direname}/flat_as.fits')

        star_blzlstname = os.path.join(direname, 'star_blz.lst')
        prepare_lst(starlstname, 'blz', star_blzlstname)
        delete_fits(star_blzlstname)
        iraf.imarith(f'@{star_1dslstname}','/',f'{direname}/flat_as.fits',f'@{star_blzlstname}')

        if has_iodn:
           iodn_blzlstname = os.path.join(direname, 'iodn_blz.lst')
           prepare_lst(iodnlstname, 'blz', iodn_blzlstname)
           delete_fits(iodn_blzlstname)
           iraf.imarith(f'@{iodn_1dslstname}','/',f'{direname}/flat_as.fits',f'@{iodn_blzlstname}')

    # delete deprecated files
    for row in flat_exptime:
        color = row[0]
        os.remove(f'{direname}/flat_{color}_ovr.lst')

    #delete_fits('lamp_sum1.lst')
    #delete_fits('star_sum1.lst')
    #os.remove(f'{star_sumlstname}')
    #os.remove(f'{lamp_sumlstname}')
    #os.remove(f'{star_ovrlstname}')
    #os.remove(f'{star_fltlstname}')
    #os.remove(f'{star_bkglstname}')
    #os.remove(f'{star_sumlstname}')
    #os.remove(f'{star_1dslstname}')
    #if os.path.exists(star_blzlstname):
    #    os.remove(star_blzlstname)

    #os.remove(lamp_ovrlstname)
    #os.remove(lamp_fltlstname)
    #os.remove(lamp_sumlstname)
    #os.remove(lamp_1dslstname)

    if has_iodn:
        os.remove(iodn_ovrlstname)
        os.remove(iodn_fltlstname)
        os.remove(iodn_bkglstname)
        os.remove(iodn_sumlstname)
        os.remove(iodn_1dslstname)
        if os.path.exists(iodn_blzlstname):
            os.remove(iodn_blzlstname)

    ## update fits headers
    update_header(logfile=logfile)

    print('Finished')
    print(f'direname = "{direname}"')
