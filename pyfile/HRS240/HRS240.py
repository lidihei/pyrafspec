from pyrafspec.main import *
import os
from pyrafspec import bfosclog, default_logheader

logfile = '20210207_test.obslog'
#########################################################################
# substract bias by overscan or bias files
#logfile = logfile.replace('log', 'obslog')
#overscan(logfile = logfile) # by overscan

#biassubtract(logfile = logfile) #subtract bias by the combined bias file
#------------------------------------------------------------------------

######################################################################
#list file
log = read_log(logfile)

#log.save_file(f'{direname}/lamp.lst', object='FeAr', exptime=300)
#log.save_file(f'{direname}/star.lst', object='Star')
#flat_lst = log.get_filenamelist(
#           object  = 'flat',
#           exptime = 900)

#for item in log.item_list:
#    print(item.object, item.skip)
#-----------------------------------------------------------------------

########################################################################
# extract spec
#direname = dire
#color = 'Blue'

#extract1d(logfile=logfile, lamp='thar')

from pyrafspec import extract2

extract = extract2.extract1d(logfile=logfile, lamp='ThAr')
#extract.process_flat(mosaicflats=True)
#extract.extract(apedit=True, apsum=True, flatfitsname=None)
## do blaze
#extract.do_blaze(starlstname=None, direname = None)
## manually calibrate wavelength
#extract.calib_wave_manual()

outfix = '_1d.fits'
suffix='_sum.fit'
#extract.calib_wave_all(lampsumlst=None, suffix=suffix, outfix=outfix)


lamplstname = os.path.join(extract.direname, 'lamp.lst')
lamplists = extract.lstname2list(lstname, suffix='_1d.fits')
starlstname = os.path.join(extract.direname, 'star.lst')
starlists = extract.lstname2list(lstname, suffix=suffix)
lamplists = [lamplists[0]]*len(starlists)
extract.match_object2lamp(lampsumlst, starsumlst=starsumlst, suffix=suffix,
                          outfix=outfix)


# check lamp
lampdata = fits.open(lamplists[0])
plt.plot(lampdata[1].T, lampdata[0].T)

# check star
starlists = extract.lstname2list(lstname, suffix=outfix)
lampdata = fits.open(starlists)
plt.plot(lampdata[1].T, lampdata[0].T)
