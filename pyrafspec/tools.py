import numpy as np

def get_ec_info(filename):
    '''get information of the ec file'''

    file = open(filename)
    find = False
    for row in file:
        row = row.strip()
        if len(row)==0 or row[0]=='#':
            continue
        g = row.split()
        if g[0]=='begin':
            # initialize lists
            use_dif_lst = []
            use_pix_lst = []
            use_wv_lst  = []
            use_ord_lst = []

            rej_dif_lst = []
            rej_pix_lst = []
            rej_wv_lst  = []
            rej_ord_lst = []

            ap_lst = []
        elif g[0]=='features':
            find = True
            continue
        elif g[0]=='offset':
            find = False
            continue
        elif find:
            ap    = int(g[0])
            order = int(g[1])
            pix   = float(g[2])
            wv_fit= float(g[3])
            if g[4]=='INDEF':
                continue
            wv_map= float(g[4])
            _     = float(g[5])
            _     = int(g[6])
            use   = int(g[7])
            if use == 0:
                rej_dif_lst.append(wv_fit-wv_map)
                rej_pix_lst.append(pix)
                rej_ord_lst.append(order)
            elif use == 1:
                use_dif_lst.append(wv_fit-wv_map)
                use_pix_lst.append(pix)
                use_ord_lst.append(order)
                use_wv_lst.append(wv_map)
        else:
            continue
    file.close()
    #return use_ord_lst, rej_ord_lst
    use_dif_lst = np.array(use_dif_lst)
    use_pix_lst = np.array(use_pix_lst)
    use_ord_lst = np.array(use_ord_lst)
    use_wv_lst  = np.array(use_wv_lst)

    rej_dif_lst = np.array(rej_dif_lst)
    rej_pix_lst = np.array(rej_pix_lst)
    rej_ord_lst = np.array(rej_ord_lst)

    wvstd = use_dif_lst.std()

    rv_dif_lst = use_dif_lst/use_wv_lst*299792458.
    rvstd = rv_dif_lst.std()*1e-3

    nlines = use_dif_lst.size+rej_dif_lst.size
    nused = use_dif_lst.size
    print(file,use_ord_lst,rej_ord_lst)
    _n_use = len(use_ord_lst)
    _n_rej = len(rej_ord_lst)
    if  _n_use > 0 and _n_rej >0:
       order_from = np.min(use_ord_lst.min(), rej_ord_lst.min())
       order_to   = np.max(use_ord_lst.max(), rej_ord_lst.max())
    elif _n_use > 0 and _n_rej ==0:
       order_from = np.min(use_ord_lst)
       order_to   = np.max(use_ord_lst)
    elif _n_use == 0 and _n_rej > 0:
       order_from = np.min(rej_ord_lst)
       order_to   = np.max(rej_ord_lst)
    norders = order_to - order_from + 1

    wv_min = use_wv_lst.min()
    wv_max = use_wv_lst.max()

    return wvstd,rvstd,nused,nlines,norders,wv_min,wv_max
