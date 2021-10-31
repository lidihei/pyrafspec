from pyrafspec.main import *
import os
from pyrafspec import bfosclog, default_logheader

dire = '/home/lcq/media/backup/216BFOSC/20211022_bfosc_pyraf'
logname = 'liuchao_bfosc.log'
logfile = os.path.join(dire, logname)


####################################################################################
# convert the log file produced by 216 cm to the format which can be used by pyrafspec (*.obslog)
date = '2021-10-22'
program = 'sdB'
observer = 'Chao Liu'
operator = 'Junjun Jia'
observatory = 'NAOC-Xinglong'
telescope = '216cm'
instrument = 'BFOSC'
slitwidth = '16"'
readout = 'Left Top & Bottom'
gain = '1.0 e-/ADU'
timesystem = 'UTC'
timezone = '+08:00'

rewritelog = bfosclog.convertobslog(dire,logname)
rewritelog.convert2obslog(date=date, program=program, observer=observer,
                          operator=operator,observatory=observatory,
                          telescope=telescope,instrument=instrument,
                          slitwidth= slitwidth, readout=readout,
                          gain = gain,timesystem=timesystem,timezone=timezone,
                          )

#-----------------------------------------------------------------------------

#########################################################################
# substract bias by overscan or bias files
logfile = logfile.replace('log', 'obslog')
#overscan(logfile = logfile) # by overscan

biassubtract(logfile = logfile) #by the combine bias file
#------------------------------------------------------------------------

######################################################################
#list file
#log = read_log(logfile)
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

extract1d(logfile=logfile)
