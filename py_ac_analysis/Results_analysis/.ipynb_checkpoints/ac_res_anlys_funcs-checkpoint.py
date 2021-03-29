import numpy as np
import pandas as pd
import scipy.io as sio

from bokeh.plotting import figure, show, save
from bokeh.io import output_notebook, output_file, reset_output
from bokeh.layouts import gridplot, row, column
from bokeh.io import export_png
output_notebook()



def loadThickness(filename):
    InitTh,Th,ThChange,ThChangeLPF,ThLPF,idx_original_thickness = list(map(lambda x: sio.loadmat(filename)[x], ['InitTh','Th','ThChange','ThChangeLPF','ThLPF','idx_original_thickness']))
    
    return InitTh,Th,ThChange,ThChangeLPF,ThLPF,idx_original_thickness


def loadResFullrun(filename):
    acous_vars = ['Amp', 'AmpRef', 'Filter', 'FreqQ', 'LocalAcTime', 'MaxInter', 'NZPad', 'NtoStack', 'RmsAmp', 'RmsAmpRef', 'TimeRange', 'TimeShift', 'freqQAmp', 'fullWFref', 'idx2analyze', 'maxAmp', 'maxFreq', 'ref', 'threshold']
    Amp, AmpRef, Filter, FreqQ, LocalAcTime, MaxInter, NZPad, NtoStack, RmsAmp, RmsAmpRef, TimeRange, TimeShift, freqQAmp, fullWFref, idx2analyze, maxAmp, maxFreq, ref, threshold = list(map(lambda x: sio.loadmat(filename, squeeze_me=True)[x], acous_vars))

    return [Amp, AmpRef, Filter, FreqQ, LocalAcTime, MaxInter, NZPad, NtoStack, RmsAmp, RmsAmpRef, TimeRange, TimeShift, freqQAmp, fullWFref, idx2analyze, maxAmp, maxFreq, ref, threshold]


def loadAcSettings(AcSettingsfile):
#     AcSettingsfile = '../mat_ac_analysis/p5483_'+run+'.mat';

    acSettings = sio.loadmat(AcSettingsfile); # load acoustic settings
    numCHR = np.size(acSettings['channels2save'][0]); # number of channels
    numCHT = np.size(acSettings['channels2transmit'][0]); # number of channels
    Nsamples = acSettings['Nsamples'][0,0]; # waveform length
    fs = float(acSettings['samplingFreq'][0,0]); # acoustic sampling rate in MHz
#     del acSettings
    return acSettings, numCHR, numCHT, Nsamples, fs


def interpDat(Time, LocalAcTime, numCHR, numCHT, MaxInter, TimeShift, RmsAmp, freqQAmp, maxAmp, maxFreq):

    # Find sample number corresponding to the beginning and end of the acoustic run

    FirstIdxAc = np.where(Time > LocalAcTime[0,0])[0][0]
    LastIdxAc = np.where(Time < LocalAcTime[-1,0])[0][-1]
    idxAc = np.arange(FirstIdxAc-1,LastIdxAc)

    MaxInterI = np.zeros((len(Time[idxAc]),numCHR,numCHT))
    TimeShiftI = np.zeros((len(Time[idxAc]),numCHR,numCHT))

    RmsAmpI = np.zeros((len(Time[idxAc]),numCHR,numCHT))
    freqQAmpI = np.zeros((len(Time[idxAc]),numCHR,numCHT))
    maxAmpI = np.zeros((len(Time[idxAc]),numCHR,numCHT))
    maxFreqI = np.zeros((len(Time[idxAc]),numCHR,numCHT))

    # Interpolate

    for chnumr in range(numCHR):
        for chnumt in range(numCHT):
            MaxInterI[:,chnumr,chnumt] = np.interp(Time[idxAc], LocalAcTime[:,chnumt], MaxInter[:,chnumr,chnumt])
            TimeShiftI[:,chnumr,chnumt] = np.interp(Time[idxAc], LocalAcTime[:,chnumt], TimeShift[:,chnumr,chnumt])

            RmsAmpI[:,chnumr,chnumt] = np.interp(Time[idxAc], LocalAcTime[:,chnumt], RmsAmp[:,chnumr,chnumt])
            freqQAmpI[:,chnumr,chnumt] = np.interp(Time[idxAc], LocalAcTime[:,chnumt], freqQAmp[:,chnumr,chnumt])
            maxAmpI[:,chnumr,chnumt] = np.interp(Time[idxAc], LocalAcTime[:,chnumt], maxAmp[:,chnumr,chnumt])
            maxFreqI[:,chnumr,chnumt] = np.interp(Time[idxAc], LocalAcTime[:,chnumt], maxFreq[:,chnumr,chnumt])

    return idxAc, MaxInterI, TimeShiftI, RmsAmpI, freqQAmpI, maxAmpI, maxFreqI


def TSplot(idxAc, Time, effNS, ThChangeLPF, TimeShift, T, R, runname, run, SAVE):

    fig11 = figure(tools='pan,box_zoom,undo,hover') #, y_axis_type="log" , 
    fig11.line(Time[idxAc], effNS[idxAc], line_width = 1.5, line_color="crimson")
    fig11.yaxis.axis_label = 'effNS (MPa)'

    fig12 = figure(x_range=fig11.x_range, tools='pan,box_zoom,undo,save,hover')
    fig12.line(Time[idxAc], ThChangeLPF[idxAc], line_width=1.5, line_color="black")
    fig12.yaxis.axis_label = 'Norm. Disp. (um)'

    fig13 = figure(title='T'+str(T+1)+' --> R'+str(R+1),x_range=fig11.x_range, tools='pan,box_zoom,undo,save,hover')
    fig13.line(Time[idxAc], TimeShift[:,R,T], line_width=1.5, line_color="mediumblue")
    fig13.yaxis.axis_label = 'Time Shift (us)'
    fig13.xaxis.axis_label = 'Time (s)'

    fig1 = gridplot([fig11, fig12, fig13], ncols=1, plot_width=800, plot_height=333)
    
    if SAVE == 1:
        reset_output()
        filename = '../Results_wAmp/'+runname+'_'+run+'_timeshift_unFilt_T'+str(T+1)+'R'+str(R+1)
        output_file(filename+'.html')
        save(fig1)
        show(fig1, notebook_handle=True)
    else: 
        reset_output()
        output_notebook()
        show(fig1)
    
    
def pickTOF(idx2analyze, T, R, Nsamples, fs, fullWFref, runname, run, SAVE):
    
    winanalysis = idx2analyze[:,R,T];
    timeWF = np.arange(Nsamples)/fs;
    partWFref = np.full(Nsamples, np.nan)
    partWFref[winanalysis[0]:winanalysis[1]] = fullWFref[winanalysis[0]:winanalysis[1],R,T]

    # figure used to pick arrival by hand
    fig2 = figure(plot_width=800, plot_height=400, tools='pan,box_zoom,undo,hover,crosshair') 
    fig2.line(timeWF, fullWFref[:,R,T], line_width=2, line_color="indigo")
    fig2.line(timeWF, partWFref, line_width=2, line_color="fuchsia")
    fig2.yaxis.axis_label = 'Reference WF (.)'
    fig2.xaxis.axis_label = 'Time (us)'
    
    if SAVE == 1:
        reset_output()
        filename = '../Results_wAmp/'+runname+'_'+run+'_refWF_unFilt_T'+str(T+1)+'R'+str(R+1)
        output_file(filename+'.html')
        save(fig2)
        show(fig2, notebook_handle=True)
    else: 
        reset_output()
        output_notebook()
        show(fig2)
    
    return timeWF