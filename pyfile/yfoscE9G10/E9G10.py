from pyrafspec.main import *
import os
from pyrafspec import bfosclog, default_logheader
dire = '/home/lcq/media/backup/240cm/yf01_reduced/20210104/'
logname = '20210104_E9G10.log'
logfile = os.path.join(dire, logname)

###############################
## get lst of la
#starlist, lamplist = bfosclog.match_star2lamp(logfile)
#starlstname_wv = f'{dire}/starlstname_wv.lst'
#starlist = [f'{dire}/{_i}.fits\n' for _i in starlist]
#fwv = open(starlstname_wv, 'w')
#fwv.writelines(starlist)
#fwv.close()
#lamplstname_wv = f'{dire}/lamplstname_wv.lst'
#lamplist = [f'{dire}/{_i}.fits\n' for _i in lamplist]
#fwv = open(lamplstname_wv, 'w')
#fwv.writelines(lamplist)
#fwv.close()



#########################################################################
# substract bias by overscan or bias files
logfile = logfile.replace('log', 'obslog')
#overscan(logfile = logfile) # by overscan

#biassubtract(logfile = logfile, rot90=True) #by the combine bias file
#------------------------------------------------------------------------

######################################################################
#list file
log = read_log(logfile)
#item = log.item_list[0]
direname = dire
#log.save_file(f'{direname}/lamp.lst', object='FeAr', exptime=300)
#log.save_file(f'{direname}/star.lst', object='Star')
#flat_lst = log.get_filenamelist(
#           object  = 'flat',
#           exptime = 900)

#for item in log.item_list:
#    print(item.object, item.skip)
#-----------------------------------------------------------------------

#######################################################################
## extract spec
#direname = dire

#calib = False
extract1d(calib=False, apsum=True, logfile=logfile, apedit=True, lamp='HeNe', lamp_exptime=90, 
         starlstname_wv=None, 
         lamplstname_wv = None, refspsort=None,
         lampreference=None)

#########################################################################


'''
from pyraf import iraf
lamp='HeNe'
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.echelle(_doprint=0)
ref_lamp = 'lamp-ref.fits'
iraf.ecid.unlearn()
iraf.ecid.maxfeat = 100 
iraf.ecid.coordli = 'linelists$hene.dat' 
iraf.ecid(ref_lamp)

ref_name = ref_lamp.replace('.fits','') 
iraf.ecreid.unlearn()
iraf.ecreid.logfile = 'STDOUT, iraf.log'
iraf.ecreid('@lamp_sum_tmp.lst', ref_name)


star_sumlstname_wv = os.path.join(dire, 'raw/star_sum.lst')
iraf.dispcor.unlearn()
iraf.refsp.referen = f'@lamp_sum_tmp.lst'
#iraf.refsp.group   = ''
iraf.refsp.sort    = ''
iraf.refsp.time    = 'no' 
iraf.refsp.logfile = 'STDOUT, iraf.log'
#iraf.refsp(f'@{star_sumlstname_wv}')
lamp_sumlstname_wv = 'lamp_sum_tmp.lst'
iraf.refsp(f'@{star_sumlstname_wv}')
'''
