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


def extract1d(calib=False, logfile=None):
    from pyraf import iraf

    # remove the old iraf.log file
    if not calib and os.path.exists('iraf.log'):
        os.remove('iraf.log')

    if not calib:
        os.system('mkiraf')
    if (logfile is None) or (logfile is ''):
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
    log.save_file(f'{direname}/lamp.lst', object='lamp')
    log.save_file(f'{direname}/star.lst', object='Star')
    if has_iodn:
        log.save_file(f'{direname}/iodn.lst', object='Iodn')

    # make flat list of different exposure times
    conf = load_conf('config/instrument.conf')
    flat_exptime = conf['flat_exptime']
    apsize = conf['apsize']
    print(f'lijiao flat_exptime = {flat_exptime}')
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
           iraf.imcomb(f'@{direname}/flat_{color}_ovr.lst', f'{direname}/flat.fits')
        else:
           iraf.imcomb(f'@{direname}/flat_{color}_ovr.lst', f'{direname}/flat_{color}.fits')

    if calib and os.path.exists(f'{direname}/flat.fits'):
        pass
    else:
        # mosaic flats
        print( )
        print('Mosaic Flats')
        print( )
        if not flatcoloris1:
           mosaic_flat(f'{direname}/flat.fits')

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
        #os.system(f'cp {direname}/flat.fits ./')
        iraf.apedit(os.path.join(direname, 'flat.fits'))
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
    if not calib:
        #os.system(f'cp {direname}/flat.fits ./')
        iraf.aptrace(os.path.join(direname, 'flat.fits'))
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
    if not calib:
        delete_fits(f'flat_n.fits')
        #os.system(f'cp {direname}/flat.fits ./')
        iraf.apflatten(os.path.join(direname, 'flat.fits'),os.path.join(direname, 'flat_n.fits')
        #os.system('rm flat.fits')

    iraf.imarith.unlearn()

    # parse lamp
    print('----------------------------parse lamp-----------------------------')
    prepare_lst(f'{direname}/lamp.lst', 'ovr', 'lamp_ovr.lst')
    prepare_lst(f'{direname}/lamp.lst', 'flt', 'lamp_flt.lst')
    if not calib:
        delete_fits('lamp_flt.lst')
        iraf.imarith(f'@lamp_ovr.lst','/',f'{direname}/flat_n.fits','@lamp_flt.lst')

    # parse star
    print('----------------------------parse star-----------------------------')
    prepare_lst(f'{direname}/star.lst', 'ovr', 'star_ovr.lst')
    prepare_lst(f'{direname}/star.lst', 'flt', 'star_flt.lst')
    if not calib:
        delete_fits('star_flt.lst')
        iraf.imarith('@star_ovr.lst','/',f'{direname}/flat_n.fits','@star_flt.lst')

    # parse iodine
    print('----------------------------parse iodine-----------------------------')
    if has_iodn:
        prepare_lst(f'iodn.lst', 'ovr', f'iodn_ovr.lst')
        prepare_lst(f'iodn.lst', 'flt', f'iodn_flt.lst')
        if not calib:
            delete_fits('iodn_flt.lst')
            iraf.imarith('@iodn_ovr.lst','/',f'{direname}/flat_n.fits','@iodn_flt.lst')

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
    iraf.apscatter.referen = os.path.join(direname, 'flat.fits')
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

    prepare_lst(f'{direname}/star.lst', 'bkg', f'star_bkg.lst')
    if not calib:
        delete_fits(f'star_bkg.lst')
        iraf.apscatter(f'@star_flt.lst', f'@star_bkg.lst')

    if has_iodn:
        prepare_lst(f'iodn.lst', 'bkg', f'iodn_bkg.lst')
        if not calib:
            delete_fits(f'iodn_bkg.lst')
            iraf.apscatter(f'@iodn_flt.lst', f'@iodn_bkg.lst')

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
    iraf.apsum.referen = os.path.join(direname, 'flat.fits')
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
    prepare_lst(f'{direname}/star.lst','sum',f'star_sum.lst')
    if not calib:
        delete_fits(f'star_sum.lst')
        iraf.apsum(f'@star_bkg.lst')

    iraf.apsum.output  = f'@lamp_sum.lst'
    prepare_lst(f'{direname}/lamp.lst','sum',f'lamp_sum.lst')
    if not calib:
        delete_fits(f'lamp_sum.lst')
        iraf.apsum(f'@lamp_flt.lst')

    if has_iodn:
        iraf.apsum.output  = f'@iodn_sum.lst'
        prepare_lst(f'iodn.lst','sum',f'iodn_sum.lst')
        if not calib:
            delete_fits(f'iodn_sum.lst')
            iraf.apsum(f'@iodn_bkg.lst')

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
        # in lamp_sum.lst
        file = open('lamp_sum.lst')
        for row in file:
            fn = row.strip()
            _direname = os.path.dirname(fn)
            _basename = os.path.dirname(fn)
            ecfn = f'{_direname}/ec'+_basename.replace('.fits','')
            if os.path.exists(ecfn):
                os.remove(ecfn)
        file.close()

        # list files and their exptimes in lamp_sum.lst
        ref_lst = {}
        no = 1
        file = open(f'lamp_sum.lst')
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

        # Now the user have to choose a lamp file to identify wavelengths
        while (True):
            n = input('Select a ThAr frame for wavelength identification [1]:')
            try:
                n = int(n)
                if 0 < n <= len(ref_lst):
                    break
            except:
                pass

        ref_lamp = ref_lst[int(n)]
        ref_name = ref_lamp.replace('.fits','')

        iraf.ecid.unlearn()
        iraf.ecid.maxfeat = 100
        iraf.ecid(ref_lamp)
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

    copy_lstfile('lamp_sum.lst', 'lamp_sum1.lst', direcp = './') 
    copy_lstfile('star_sum.lst', 'star_sum1.lst', direcp = './') 
    iraf.ecreid.unlearn()
    iraf.ecreid.logfile = 'STDOUT, iraf.log'
    print(f'ref_name = "{ref_name}"')
    #iraf.ecreid('@lamp_sum.lst', ref_name)
    iraf.ecreid('@lamp_sum1.lst', ref_name)

    _ = input('Press [Enter] to continue: ')

    iraf.refsp.unlearn()
    iraf.refsp.referen = f'@lamp_sum1.lst'
    iraf.refsp.sort    = 'DATE-STA'
    iraf.refsp.group   = ''
    iraf.refsp.time    = 'yes'
    iraf.refsp.logfile = 'STDOUT, iraf.log'
    iraf.refsp(f'@star_sum1.lst')
    iraf.refsp(f'@lamp_sum1.lst')
    if has_iodn:
        copy_lstfile('iodn_sum.lst', 'iodn_sum1.lst', direcp = './') 
        iraf.refsp(f'@iodn_sum.lst')


    iraf.dispcor.unlearn()
    iraf.dispcor.lineari = 'no'
    iraf.dispcor.flux    = 'no'
    iraf.dispcor.logfile = 'iraf.log'
    prepare_lst(f'{direname}/star.lst', '1ds', f'star_1ds.lst')
    delete_fits(f'star_1ds.lst')
    iraf.dispcor(f'@star_sum1.lst', f'@star_1ds.lst')

    prepare_lst(f'{direname}/lamp.lst', '1ds', f'lamp_1ds.lst')
    delete_fits(f'lamp_1ds.lst')
    iraf.dispcor(f'@lamp_sum1.lst',f'@lamp_1ds.lst')

    if has_iodn:
        prepare_lst(f'iodn.lst', '1ds', f'iodn_1ds.lst')
        delete_fits(f'iodn_1ds.lst')
        iraf.dispcor(f'@iodn_sum1.lst', f'@iodn_1ds.lst')
        delete_fits('iodn_sum1.lst')


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

        prepare_lst(f'{direname}/star.lst', 'blz', f'star_blz.lst')
        delete_fits(f'star_blz.lst')
        iraf.imarith(f'@star_1ds.lst','/',f'{direname}/flat_as.fits',f'@star_blz.lst')

        if has_iodn:
            prepare_lst('iodn.lst', 'blz', 'iodn_blz.lst')
            delete_fits('iodn_blz.lst')
            iraf.imarith(f'@iodn_1ds.lst','/',f'{direname}/flat_as.fits',f'@iodn_blz.lst')

    # delete deprecated files
    for row in flat_exptime:
        color = row[0]
        os.remove(f'{direname}/flat_{color}_ovr.lst')

    delete_fits('lamp_sum1.lst')
    delete_fits('star_sum1.lst')
    os.remove(f'star_sum1.lst')
    os.remove(f'lamp_sum1.lst')
    os.remove(f'star_ovr.lst')
    os.remove(f'star_flt.lst')
    os.remove(f'star_bkg.lst')
    os.remove(f'star_sum.lst')
    os.remove(f'star_1ds.lst')
    if os.path.exists(f'star_blz.lst'):
        os.remove(f'star_blz.lst')

    os.remove(f'lamp_ovr.lst')
    os.remove(f'lamp_flt.lst')
    os.remove(f'lamp_sum.lst')
    os.remove(f'lamp_1ds.lst')

    if has_iodn:
        os.remove(f'iodn_ovr.lst')
        os.remove(f'iodn_flt.lst')
        os.remove(f'iodn_bkg.lst')
        os.remove(f'iodn_sum.lst')
        os.remove(f'iodn_1ds.lst')
        if os.path.exists('iodn_blz.lst'):
            os.remove('iodn_blz.lst')

    ## update fits headers
    update_header(logfile=logfile)

    print('Finished')
    print(f'direname = "{direname}"')
