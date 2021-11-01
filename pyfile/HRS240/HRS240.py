from pyrafspec.main import *
import os
from pyrafspec import bfosclog, default_logheader

logfile = '20210207.obslog'
#########################################################################
# substract bias by overscan or bias files
logfile = logfile.replace('log', 'obslog')
#overscan(logfile = logfile) # by overscan

biassubtract(logfile = logfile) #subtract bias by the combined bias file
#------------------------------------------------------------------------

######################################################################
#list file
log = read_log(logfile)
direname = dire
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

extract1d(logfile=logfile, lamp='thar')
