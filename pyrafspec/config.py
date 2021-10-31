import os

def load_conf(filename):
    conf = {}
    file = open(filename)
    for row in file:
        row = row.strip()
        if len(row)==0 or row[0]=='#':
            continue
        g = row.split('=')
        key   = g[0].strip()
        value = g[1].split('#')[0].strip()
        if key[-12:]=='flat exptime':
            color = key.split()[0]
            try:
                exptime = int(value)
            except:
                exptime = float(value)
            if 'flat_exptime' not in conf:
                conf['flat_exptime'] = []
            conf['flat_exptime'].append([color,exptime])
        elif key=='ccd map':
            conf['ccd_map'] = value
        elif key=='apsize':
            conf['apsize'] = int(value)
    file.close
    return conf
