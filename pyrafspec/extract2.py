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
from twodspec.polynomial import Poly1DFitter
from twodspec.trace import trace_naive_max
from astropy import table
from astropy.io import fits
from twodspec import thar

class extract1d():
      
      def __init__(self, logfile, lamp='ThAr', lamp_exptime=None,
                         instrumentsfile =None, direname=None):
          log = read_log(logfile)
          print('Obs Log File =',logfile)
          # make list of different types
          direname = os.path.dirname(log.filename_composition)
          log.save_file(f'{direname}/lamp.lst', object=lamp, exptime=lamp_exptime)
          log.save_file(f'{direname}/star.lst', object='Star')
          has_iodn = log.has_object('Iodn')
          self.has_iodn = has_iodn
          if has_iodn:
              log.save_file('{ename}/iodn.lst', object='Iodn')
          if instrumentsfile is None: instrumentsfile = './conf/instrument.conf'
          conf = load_conf('config/instrument.conf')
          flat_exptime = conf['flat_exptime']
          flatfitsname = os.path.join(direname, 'flat.fits')
          self.flatfitsname = flatfitsname
          self.direname = direname
          self.flat_exptime = conf['flat_exptime']
          self.apsize = conf['apsize']
          self.log = log

      def initialize_iraf(self):
          from pyraf import iraf
          # initialize IRAF
          print( )
          print('Initialize IRAF')
          print( )
          iraf.noao(_doprint=0)
          iraf.imred(_doprint=0)
          iraf.echelle(_doprint=0)
          return iraf


      def process_flat(self, mosaicflats=True):
          iraf = self.initialize_iraf()
          #make list of different exposure times
          print(self.flat_exptime)
          direname = self.direname
          log = self.log
          flatfitsname = self.flatfitsname
          flat_exptime = self.flat_exptime
          for row in flat_exptime:
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
          if mosaicflats:
             # mosaic flats
             print( )
             print('Mosaic Flats')
             print( )
             mosaic_flat(flatfitsname)

      def extract(self, apedit=True, apsum=True, flatfitsname=None):
          iraf = self.initialize_iraf()
          has_iodn = self.has_iodn
          apsize = self.apsize
          if flatfitsname is None:flatfitsname = self.flatfitsname 
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
          if apedit:
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
          if  apedit:
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
          if  apedit:
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
          if apedit:
              delete_fits(lamp_fltlstname)
              iraf.imarith(f'@{lamp_ovrlstname}','/',flat_nfitsname,f'@{lamp_fltlstname}')
          # parse star
          print('----------------------------parse star-----------------------------')
          starlstname = os.path.join(direname, 'star.lst')
          self.starlstname = starlstname
          star_ovrlstname = os.path.join(direname, 'star_ovr.lst')
          star_fltlstname = os.path.join(direname,'star_flt.lst')
          prepare_lst(starlstname, 'ovr', star_ovrlstname)
          prepare_lst(starlstname, 'flt', star_fltlstname)
          if apedit:
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
              if apedit:
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
          if apsum:
              delete_fits(star_sumlstname)
              iraf.apsum(f'@{star_bkglstname}')

          lamplstname = os.path.join(direname, 'lamp.lst')
          lamp_sumlstname = os.path.join(direname, 'lamp_sum.lst')
          slef.lamplstname = lamplstname
          prepare_lst(lamplstname,'sum',lamp_sumlstname)
          iraf.apsum.output  = f'@{lamp_sumlstname}'
          if apedit:
              delete_fits(lamp_sumlstname)
              iraf.apsum(f'@{lamp_fltlstname}')

          if has_iodn:
              prepare_lst(iodnlstname,'sum',iodn_sumlstname)
              iodn_sumlstname = os.path.join(direname, 'iodn_sum.lst')
              iraf.apsum.output  = f'@{iodn_sumlstname}'
              if epedit:
                  delete_fits(iodn_sumlstname)
                  iraf.apsum(f'@{iodn_bkglstname}')

      def do_blaze(self, starlstname=None, direname = None):
          iraf = self.initialize_iraf()
          if starlstname is None: starlstname = self.starlstname
          if direname is None: direname = self.direname
          star_sumlstname =  os.path.join(direname,'star_sum.lst')
          prepare_lst(starlstname,'sum',star_sumlstname)
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
          iraf.imarith(f'@{star_sumlstname}','/',f'{direname}/flat_as.fits',f'@{star_blzlstname}')



      def calib_wave_manual(self, direname= None, lamplstname=None, lampreference=None):
          ''' manually calibrate wavelength by iraf
          '''
          print( )
          print('Wavelength Calibration')
          print( )
          if lamplstname is None: lamplstname = lamplstname
          if direname is None: direname = self.direname
          lamplstname = os.path.join(direname, 'lamp.lst')
          lamp_sumlstname = os.path.join(direname, 'lamp_sum.lst')
          prepare_lst(lamplstname,'sum',lamp_sumlstname)

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

      def lstname2list(self, lstname, suffix=None, split='.fit'):
          ''' obtain list from a list file
          parameters:
          -----------------
          lstname: [str]  e.g. star.lst
          suffix: [str], e.g. '.fits'
          returns:
          lists: [list]
          '''
          with open(lstname, 'r') as f:
               lists = f.readlines()
          if suffix is None: suffix = '.fits'
          lists = [f'{i.split(split)[0]}{suffix}' for i in lists]
          return lists

      def calib_wave_all(self, lampsumlst=None, suffix='_sum.fits', outfix='_1d.fits', **arguments):
          '''
          calibrate all lamp files
          lampsumlst: [list]
          '''
          if lampsumlst is None:
             lamplstname = slef.lamplstname
             lampsumlst = self.lstname2list(lamplstname, suffix=suffix)
          for lampsumfits in lampsumlst:
              self.calib_wave(lampsumlst, suffix=outfix, **arguments)

      def match_object2lamp(self, lampsumlst, starsumlst=None, suffix='_sum.fits', 
                           outfix='_1d.fits'):
          '''add wavelength to star spectra
          parameters:
          ----------------
          lampsumlst: [list]
          starsumlst: [list]   len(lampsumlst) == len(starsumlst)
          '''
          if starsumlst is None:
             starlstname = self.starlstname
             starsumlst = self.lstname2list(lamplstname, suffix=suffix)
          assert len(lampsumlst) == len(starsumlst)
          for i, fname in enumerate(starsumlst):
              lampname = lamplstname[i]
              waves  = fits.getdata(lampname)[1]
              basename = os,path.basename(lampname)
              flux,header= fits.getdata(fname, header=True)
              header['lamp'] = (basename, 'lamp name')
              hdu = fits.PrimaryHDU(np.array([flux, waves]))
              hdu.header = header
              fout = fname.replace(suffix, outfix)
              hdu.writeto(fout, overwrite=True)
              

      def calib_wave(self, lampsumfits, templatename= None, linelistname=None,
                   pw=1, deg=2, threshold=0.1, min_select=20,
                   nsigma=2.5, verbose=False, suffix='_wv', **arguments):
          '''calibrate the wavelength
          parameters:
          ----------------------
          lampsumfits [str] e.g.  20210102028_1ds.fits
          returns:
          --------------------
          calibration_dict [dictionary]
          '''
          if templatename is None:
             templatename = '../template/240HRStemplate/20210102028_1ds.fits'
          if 'dump' in templatename:
              import joblib
              dumpdata = joblib.load(templatename)
              wave_temp = dumpdata['wave']
              lamp_temp = dumpdata['flux']
          else:
              lamp_temp,head= fits.getdata(templatename,header=True)
              multispec_items = MultiSpecItem.get_wat2(head)
              orders = np.arange(len(multispec_items))+1
              wave_temp = np.zeros(lamp_temp.shape,dtype=np.float32)
              for _i, _order in enumerate(orders):
                  wave_temp[_i] = multispec_items[_order].get_wv()
          
          if linelistname is None:
             linelistname = '../linelist/thar.dat'
          linelist = np.loadtxt(linelistname)
          lamp1d, header = fits.getdata(lampsumfits,header=True)
          wave_init = thar.corr_thar(wave_temp, lamp_temp, lamp1d, maxshift=50)
          #find lamp lines
          tlines = thar.find_lines(wave_init, lamp1d, linelist, npix_chunk=20, ccf_kernel_width=1.5)
          ind_good = np.isfinite(tlines["line_x_ccf"]) & (np.abs(tlines["line_x_ccf"] - tlines["line_x_init"]) < 10) & (
                  (tlines["line_peakflux"] - tlines["line_base"]) > 100) & (
                             np.abs(tlines["line_wave_init_ccf"] - tlines["line"]) < 3)
          tlines.add_column(table.Column(ind_good, "ind_good"))
          # tlines.show_in_browser()
          #clean each order
          def clean(pw=1, deg=2, threshold=0.1, min_select=20):
              order = tlines["order"].data
              ind_good = tlines["ind_good"].data
              linex = tlines["line_x_ccf"].data
              z = tlines["line"].data

              u_order = np.unique(order)
              for _u_order in u_order:
                  ind = (order == _u_order) & ind_good
                  if np.sum(ind) > 10:
                      # in case some orders have only a few lines
                      p1f = Poly1DFitter(linex[ind], z[ind], deg=deg, pw=pw)
                      res = z[ind] - p1f.predict(linex[ind])
                      ind_good[ind] &= np.abs(res) < threshold
              tlines["ind_good"] = ind_good
              return
          
          print("  |- {} lines left".format(np.sum(tlines["ind_good"])))
          clean(pw=1, deg=2, threshold=0.8, min_select=20)
          clean(pw=1, deg=2, threshold=0.4, min_select=20)
          clean(pw=pw, deg=deg, threshold=threshold, min_select=min_select)
          print("  |- {} lines left".format(np.sum(tlines["ind_good"])))
          tlines = tlines[tlines["ind_good"]]
          # fitting grating equation
          x = tlines["line_x_ccf"]  # line_x_ccf/line_x_gf
          y = tlines["order"]
          z = tlines["line"]
          pf1, pf2, indselect = thar.grating_equation(
              x, y, z, deg=(3, 7), nsigma=nsigma, min_select=210, verbose=verbose)
          tlines.add_column(table.Column(indselect, "indselect"))
          rms = pf2.rms
          # reasonable
          nlines = np.sum(indselect)
          # mpflux
          mpflux = np.median(tlines["line_peakflux"][tlines["indselect"]])
          # rms
          rms = np.std((pf2.predict(x, y) - z)[indselect])
          print("  |- nlines={}  rms={:.4f}A  mpflux={:.1f}".format(nlines, rms, mpflux))
          # predict wavelength solution
          nx, norder = lamp1d.shape
          mx, morder = np.meshgrid(np.arange(norder), np.arange(nx))
          wave_solu = pf2.predict(mx, morder)  # polynomial fitter
          # result
          calibration_dict = collections.OrderedDict(
              wave_init=wave_init,
              wave_solu=wave_solu,
              tlines=tlines,
              nlines=nlines,
              pf1=pf1,
              pf2=pf2,
              rms = rms,
              mpflux=mpflux,
              # fear=fear,
              lamp1d=lamp1d
          )
          
          header['rms'] = rms
          hdu = fits.PrimaryHDU(np.array([lamp1d, wave_solu]))
          hdu.header = header
          fout = f'{lampsumfits.split(".fit")[0]}{suffix}.fits'
          hdu.writeto(fout, overwrite=True)
          return calibration_dict
