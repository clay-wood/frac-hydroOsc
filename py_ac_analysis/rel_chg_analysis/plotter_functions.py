import numpy as np

def marker_shape(input):
    freqs = np.array([0.1,1,10,40])
    mrkr = np.array(['d','o','*','s'])
    return mrkr[np.where(input == freqs)][0]

def mrkr_run(run):
    sets = np.array([1,2,3,4,None])
    mrkr = np.array(['o','s','d','*','None'])
    mrkrsz = np.array([6,5,6,8,1])
    return mrkr[np.where(run == sets)][0], mrkrsz[np.where(run == sets)][0]

def marker_size(input):
    freqs = np.array([0.1,1,10,40])
    mrkrsz = np.array([6,6,9,5])
    return mrkrsz[np.where(input == freqs)][0]

def mrk_amp_sz(amp):
    mrkrsz = np.exp(amp)+3;
    return mrkrsz

def oneHz_filter(x): 
        freq_locs = np.where(x[:,1] != 1)[0]
        x[freq_locs,2] = np.nan
        return x

def tenHz_filter(x, y):
    freq_locs = np.where(x[:,1] > 1)[0]
    y[freq_locs] = np.nan
    return y

def oneMPa_filter(x): 
        amp_locs = np.where(np.logical_and(x[:,2] > 0.9, x[:,2] < 1.3))[0]
        x[amp_locs,1] = np.nan
        return x

def which_osc(osc_type, x):
    if osc_type == 'NS':
        if len(x.T) > 21:
            clr = '#fdb366' #'#ee7733' #'#fb9a29'
            perm = 8; amp = 2
            kdot = 14
        else:
            clr = '#B8221E'
            perm = 7; amp = 2
            kdot = 12
    else:
        if len(x.T) > 21:
            clr = '#92C5DE'
            perm = 8; amp = 2
            kdot = 14
        else:
            clr = '#2166AC'
            perm = 7; amp = 2
            kdot = 12
    return clr, kdot, amp, perm