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


def mosaic_flat(flatname, confname=None):
    if confname is None:
       confname = 'config/instrument.conf'
    conf = load_conf(confname)
    flat_exptime = conf['flat_exptime']

    # prepare data list
    direname = os.path.dirname(flatname)
    data_lst = {}
    for row in flat_exptime:
        color = row[0]
        data = pf.getdata(f'{direname}/flat_{color}.fits')
        yrows,xrows = data.shape
        xcenter = int(xrows/2)
        flat1d = data[:,xcenter]
        data_lst[color] = flat1d

    # load boundaries from *.reg files
    bound_lst = np.zeros(xrows,dtype=np.int32)
    cwd = os.getcwd()
    for fname in os.listdir(cwd):
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
        return
    if len(flat_exptime) != len(bound_lst):
        print('number of flats = %d, while number of boundaries = %d'%(
                len(flat_exptime), len(bound_lst)))
        return

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
                colorflat = pf.getdata(f'{direname}/flat_{color}.fits')
                flat += m*colorflat

    # save white fits
    if os.path.exists(flatname):
        os.remove(flatname)
    pf.writeto(flatname,flat)
