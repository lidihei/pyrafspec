#!/usr/bin/env python
import os, sys, math, datetime, dateutil
import argparse
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
from .extract import *
from .mosaicflat import *
from .extract import *
from .subtractbias import *
from .tools import *
#from sys import exit
#import pylab
#import matplotlib 
#matplotlib.use("TkAgg") 

def check_data(logfile=None):
    ''' check data -
    the existences of frames in the log file
    the exposure times
    the offset between log file and control computer
    '''
    if logfile is None:
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

def plot_flat(logfile=None):
    '''
    plot a portion of flat frame to check the brightness
    variations
    '''

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import colorConverter

    r1 = 1100
    r2 = 1200
    if loggfile is None:
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

def _compare(filenames, colors=['k', 'b', 'r', 'c', 'y']):
    '''
    comapre the each order 1d spectra of two or more fits file
    len(colors) > len(filenames)
    parameters:
    -----------
    filenames [list]
    '''
    specs = [load_multispec(filename) for filename in filenames]
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_axes([0.1,0.15,0.85,0.75])

    def plot_order(order):
        ax.cla()
        ax.currentorder = order
        spds = [spec[order] for spec in specs]

        wvc = (spds[0].wv_from+spds[0].wv_to)/2.0
        try:
           [ax.plot(spd.wv, spd.flux,'-',color=colors[_]) for _, spd in enumerate(spds)]
        except:
           print(f'the number of colors (Ncolor = {len(colors)}) should larger than the number of filenmaes (Nfiles = {len(filenames)})')
        ax.set_xlabel(u'Wavelength (\xc5)')
        ax.set_ylabel('Flux (ADU)')
        #basename = os.path.basename(filename)
        ax.set_title('%s - Order %d'%('spec',order))
        if (spds[0].flux < 0).sum()==0:
            ax.set_ylim(0,)
        #ax.set_xlim(spd.wv_from, spd.wv_to)
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

    order = specs[0].orders[int(len(specs[0].orders)/2.)]
    spec = specs[0]
    plot_order(order)

    fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()

def convert(filetype, direname = None):
    if direname is None: direname = './'
    for fname in os.listdir(direname):
        #if fname[-9:] in ['_com.fits']:
        if fname[-9:] in ['_1ds.fits','_blz.fits']:
            spec = load_multispec(fname)
            for o in spec.orders:
                newname = '%s_o%03d.%s'%('.'.join(fname.split('.')[0:-1]),o,filetype)
                if filetype == 'txt':
                    spec[o].save_txt(newname)
                elif filetype == 'fits':
                    spec[o].save_fits(newname)
                print('%s order %03d -> %s'%(fname, o, newname))


def convert2(filetype, logfile=None):
    '''deprecated convert function'''
    if logfile is None:
       logfile = find_logfile()
    log = read_log(logfile)

    for item in log.item_list:
        if item.skip:
            continue
        filename = item.get_filename(log.filename_composition)
        for key in ['1ds','blz','com']:
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

def clean(logfile =None):
    if logfile is None:
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

def main():
    parser = argparse.ArgumentParser(description="pyrafspec: extract spectrum by using pyraf")

    subparsers = parser.add_subparsers(dest='action')

    # --clean--
    clean_parser = subparsers.add_parser('clean', help='clean the directory after 1-D extraction')
    clean_parser.add_argument('logfile', default=None,
                              help='log file name')
    clean_parser.set_defaults(func=clean)
if __name__=='__main__':

    try:
        argv = sys.argv[1]
    except:
        help()
        #import sys
        exit()

    if argv=='clean':
        if len(sys.argv)>=2:
           logfile = sys.argv[2]
        clean(logfile)

    elif argv=='check':
        # check fits files
        if len(sys.argv)>=2:
           logfile = sys.argv[2]
        check_data(logfile)

    elif argv=='flat':
        # check flats
        plot_flat()

    elif argv=='overscan':
        if len(sys.argv)>=2:
           logfile = sys.argv[2]
        else: logfile =None
        # overscan correction
        overscan(logfile=logfile)

    elif argv=='mosaic':
        # mosaic flats
        mosaic_flat('flat.fits')

    elif argv=='reduce':
        if len(sys.argv)>=2:
           logfile = sys.argv[2]
        else: logfile =None
        # reduce 1-D spectra
        extract1d(logfile=logfile)

    elif argv=='calib':
        if len(sys.argv)>=2:
           logfile = sys.argv[2]
        else: logfile =None
        extract1d(calib=True, logfile=logfile)

    elif argv=='view':
        # review the extracted 1-D spectra in a matolotlib window
        if len(sys.argv)>=3:
            filename = sys.argv[2]
            if os.path.exists(filename):
                #if filename[-9:] in ['_1ds.fits','_blz.fits']:
                if filename[-9:] in ['_1ds.fits','_blz.fits','_com.fits']:
                    view(filename)
                else:
                    print('we do not support this type of fits'%filename)
            else:
                print('%s does not exist'%filename)
        else:
            print('please give the filename')

    elif argv=='compare':
        # review the extracted 1-D spectra in a matolotlib window
        if len(sys.argv)>=3:
            filenames = sys.argv[2:]
            for filename in filenames:
                if os.path.exists(filename):
                    #if filename[-9:] in ['_1ds.fits','_blz.fits']:
                    if filename[-9:] in ['_1ds.fits','_blz.fits','_com.fits']:
                        print(f'we support this type of fits {filename}')
                    else:
                        print(f'we do not support this type of fits {filename}')
                else:
                    print('%s does not exist'%filename)
            _compare(filenames)
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
