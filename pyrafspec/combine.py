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

def combine(input_lst,outputname):

    color_lst = ['r','g','b','m','c','y','k']

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
                ifrom = int(npt*0.05)#0.2
                ito   = int(npt*0.95)#0.8

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

