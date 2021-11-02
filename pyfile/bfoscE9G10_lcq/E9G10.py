from pyrafspec.main import *
import os
from pyrafspec import bfosclog, default_logheader

dire = '/home/lcq/media/backup/216BFOSC/20211023_bfosc_lcq'
logname = 'liuchao_bfosc.log'
logfile = os.path.join(dire, logname)

##############################
# get lst of la
starlist, lamplist = bfosclog.match_star2lamp(logfile)
starlstname_wv = f'{dire}/starlstname_wv.lst'
starlist = [f'{dire}/{_i}.fits\n' for _i in starlist]
fwv = open(starlstname_wv, 'w')
fwv.writelines(starlist)
fwv.close()
lamplstname_wv = f'{dire}/lamplstname_wv.lst'
lamplist = [f'{dire}/{_i}.fits\n' for _i in lamplist]
fwv = open(lamplstname_wv, 'w')
fwv.writelines(lamplist)
fwv.close()








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

#######################################################################
# extract spec
direname = dire
color = 'Blue'

extract1d(logfile=logfile, lamp='FeAr', lamp_exptime=300)

#########################################################################
'''
#test pyraf
from pyraf import iraf
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.echelle(_doprint=0)

direname = dire
star_sumlstname =  os.path.join(direname,'star_sum.lst')
star_bkglstname =  os.path.join(direname,'star_bkg.lst')
starlstname =  os.path.join(direname,'star.lst')
flatfitsname = os.path.join(direname,'flat.fits')
lamp_sumlstname =  os.path.join(direname,'lamp_sum.lst')
lamplstname =  os.path.join(direname,'lamp.lst')
prepare_lst(starlstname,'sum',star_sumlstname)
new_ecid = True
if new_ecid:
    # delete the existing ec files corresponding to the filenames
    # in lamp_sum.lst
    file = open(lamp_sumlstname)
    for row in file:
        fn = row.strip()
        _direname = os.path.dirname(fn)
        _basename = os.path.basename(fn)
        ecfn = f'{_direname}/ec'+_basename.replace('.fits','')
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
    ref_lamp = os.path.basename(ref_lamp)
    ref_name = ref_lamp.replace('.fits','')

    iraf.ecid.unlearn()
    iraf.ecid.maxfeat = 100
    iraf.ecid.coordli = 'linelists${lamp.lower()}.dat'
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

starlstname_wv
star_sumlstname_wv = os.path.join(direname, 'star_sumlstname_wv.lst')
lamp_sumlstname_wv = os.path.join(direname, 'lamp_sumlstname_wv.lst')
#prepare_lst(starlstname, '1ds', star_1dslstname)
prepare_lst(starlstname_wv, 'sum', star_sumlstname_wv) 
lamplstname_wv
#prepare_lst(lamplstname, '1ds', lamp_1dslstname)
prepare_lst(lamplstname_wv, 'sum', lamp_sumlstname_wv)

copy_lstfile(lamp_sumlstname, 'lamp_sum_tmp.lst', direcp = './') 
#copy_lstfile(star_sumlstname, 'star_sum_tmp.lst', direcp = './') 
iraf.ecreid.unlearn()
iraf.ecreid.logfile = 'STDOUT, iraf.log'
print(f'ref_name = "{ref_name}"')
#iraf.ecreid(f'@{lamp_sumlstname}', ref_name)
iraf.ecreid('@lamp_sum_tmp.lst', ref_name)

_ = input('Press [Enter] to continue: ')

iraf.refsp.unlearn()
iraf.refsp.referen = f'@lamp_sum_tmp.lst'
#iraf.refsp.sort    = 'DATE-STA'
iraf.refsp.group   = ''
iraf.refsp.time    = 'yes'
iraf.refsp.logfile = 'STDOUT, iraf.log'
iraf.refsp(f'@{star_sumlstname_wv}')
iraf.refsp(f'@{lamp_sumlstname_wv}')
has_iodn = False
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
iraf.dispcor(f'@{star_sumlstname}', f'@{star_1dslstname}')

lamp_1dslstname = os.path.join(direname, 'lamp_1ds.lst')
prepare_lst(lamplstname, '1ds', lamp_1dslstname)
delete_fits(lamp_1dslstname)
iraf.dispcor(f'@{lamp_sumlstname}',f'@{lamp_1dslstname}')
umlstname
'''
