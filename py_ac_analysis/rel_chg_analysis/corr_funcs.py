import numpy as np
import pandas as pd
import scipy.io as sio
from sklearn.metrics import r2_score
from scipy.signal import butter,filtfilt# Filter requirements.

from bokeh.plotting import figure, show, save
from bokeh.io import output_notebook, output_file, reset_output
from bokeh.layouts import gridplot, row, column
output_notebook()


def osc_beg_end(start, num_osc_sets):
    cycles = 20
    rec_freq = 100 #Hz
    hold_time = 90 * rec_freq # rec num
    osc_sets = np.ones(num_osc_sets)

    pp_osc_recs = cycles/osc_sets * rec_freq # num / Hz * Hz = rec num
    
    pp_start = np.ones(np.size(osc_sets))*start
    pp_end = np.ones(np.size(osc_sets))*start + pp_osc_recs[0]
    
    pp_start[1::] = pp_start[1::] + np.cumsum(pp_osc_recs[0:-1] + hold_time)
    pp_end[1::] = pp_end[0::-1] + np.cumsum(pp_osc_recs[1::] + hold_time)
    
    return pp_start.astype(int), pp_end.astype(int)


def relChgPct(param, time_before, time_after, pp_start, pp_end):
    rec_freq = 100 #Hz
    rec_before = time_before * rec_freq #rec num
    rec_after = time_after * rec_freq #rec num
    
    before = pp_start - rec_before
    after = pp_end + rec_after
    
    p0 = np.array(list(map(lambda x: np.nanmean(np.abs(param)[before[x]:pp_start[x]]), np.arange(np.size(pp_start)))))
    p1 = np.array(list(map(lambda x: np.nanmean(np.abs(param)[pp_end[x]:after[x]]), np.arange(np.size(pp_end)))))

    relchg = (p1-p0)/p0 * 100
    
    return relchg


def osc_amp(param, pp_start, pp_end):
    amp_max = np.array(list(map(lambda x: np.nanmax(param[pp_start[x]:pp_end[x]]), np.arange(np.size(pp_start)))))
    amp_min = np.array(list(map(lambda x: np.nanmin(param[pp_start[x]:pp_end[x]]), np.arange(np.size(pp_end)))))
    
    amps = np.abs(amp_max - amp_min)/2
    
    return amps


def recovFitter(pp_start, pp_end, Time, param, label):

    pp_start = np.append(pp_start[1::], pp_end[-1]+(90*100))
    
    p = np.zeros((len(pp_end), 3))
    q = np.empty((len(pp_end),90*100))
    
    for aa in range(len(pp_end)):
        fin_idx = np.isfinite(param[pp_end[aa]:pp_start[aa]])
        p[aa,0:2] = np.polyfit(np.log10(Time[pp_end[aa]:pp_start[aa]][fin_idx]-Time[pp_end[aa]:pp_start[aa]][fin_idx][0]+0.01), param[pp_end[aa]:pp_start[aa]][fin_idx], 1)
        q[aa,:] = np.polyval(p[aa,0:2], np.log10(Time[pp_end[aa]:pp_start[aa]]-Time[pp_end[aa]:pp_start[aa]][0]+0.01))
        p[aa,2] = r2_score(param[pp_end[aa]:pp_start[aa]][fin_idx], q[aa,:][fin_idx])

        fig = figure(title = 'r^2 = '+str(round(p[aa,2],3)), tools='pan,box_zoom,undo,hover,crosshair') 
        fig.circle(np.log10(Time[pp_end[aa]:pp_start[aa]][fin_idx]-Time[pp_end[aa]:pp_start[aa]][fin_idx][0]+0.01), param[pp_end[aa]:pp_start[aa]][fin_idx], size=5, fill_color='black', line_color="black")
        fig.line(np.log10(Time[pp_end[aa]:pp_start[aa]]-Time[pp_end[aa]:pp_start[aa]][0]+0.01), q[aa,:], line_color="red")
        fig.yaxis.axis_label = label+' recov'
        fig.xaxis.axis_label = 'Time (s)'

        fig = gridplot([fig], ncols=1, plot_width=600, plot_height=400)
        show(fig)
    
    return q, p


def detrend(x, y):
    fin_idx = np.isfinite(y)
    p = np.polyfit(x[fin_idx], y[fin_idx], 1)
    q = np.polyval(p, x)
    y_detrend = (y - q) + y[fin_idx][0]
    return y_detrend


def butter_filter(data, cutoff, fs, order, hilo):
#     print("Cutoff freq " + str(cutoff))
    nyq = 0.5 * fs # Nyquist Frequency
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype=hilo, analog=False)
    y = filtfilt(b, a,data)
    return y