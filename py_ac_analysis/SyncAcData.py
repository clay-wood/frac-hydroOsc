import numpy as np
import scipy.io as sio
import fnmatch
import os

from bokeh.plotting import figure, show
from bokeh.io import output_notebook
from bokeh.layouts import gridplot, row, column
output_notebook()

def SyncAcData(AcSettingsfile,Time,Sync,idxft,idxlt,idxref1,idxref2,runpath):
    # SyncAcData uses acoustic settings (mat file from Verasonics) and indexes
    # provided by user to sync mechanical and acoustical data. It accounts for
    # small variations in sampling frequency/period (of the order of
    # picoseconds). It uses the biax recorder (mechanical data) as a
    # reference and assume that the sampling frequency for acoustical data (and
    # therefore acoustic pulsing rate) is slightly off.

    # The main output is a time vector for acoustic data. It also outputs an
    # adjusted acoustic sampling rate and an adjusted acoustic pulsing rate.

    # INPUTS
    # AcSettingsfile is the path where pXXXX.mat (containing acoustic settings) can be loaded
    # Time is the time vector from biax data
    # Sync is the time vector from biax data
    # idxft is the index corresponding to the first trigger of the run
    # idxlt is the index corresponding to the last trigger of the run
    # idxref1 and idxref2 are two large triggers chosen towards the beginning
    # and the end of the run
    # runpath is the path to directory containing relevant WF files

    # OUTPUTS
    # acTime is the time vector for acoustic data obtained after synchronization
    # acPeriod_adjusted is the time between two consecutive acoustic pulses in microsec (found after adjusting using biax recorder as a reference)
    # ts_adjusted is the sampling period in microsec (i.e. 1/fs) found after adjusting using biax recorder as a reference
    # TotalNumberOfFiles is the total number of acoustic files to be analyzed 

    #------------------------------------------------------------------------------
    
    numAcFiles = len(fnmatch.filter(os.listdir(runpath), '*.ac')) # count number of files in WF directory

    acSettings = sio.loadmat(AcSettingsfile); # load acoustic settings
    acPeriod = float(acSettings['SeqControl'][0,0][1][0,0]); # time btw pulses in microsec (SeqControl(1).argument = timeBetweenpulses)
    acRate = 1e6/(acPeriod); # in Hz (number of WF per second)
    numSFpfile = float(acSettings['numFrames'][0,0]/2); # number of superframes per file
    numWFpSFpCH = float(acSettings['numAcqs'][0,0]); # number of WF per superframe and per channel
    numWFpfilepCH = numSFpfile*numWFpSFpCH;  # number of WF per file and per channel
    fs = float(acSettings['samplingFreq'][0,0]); # acoustic sampling rate in MHz

    del acSettings

    ts = 1/fs; # sampling time in microsec
    
    
    MAX = max(Sync[idxft:idxlt])
    MIN = min(Sync[idxft:idxlt])
    AMP = MAX-MIN
    Sync[idxft:idxlt] = (Sync[idxft:idxlt] - np.mean(Sync[idxft:idxlt]))/AMP
    Sync[idxft:idxlt] = Sync[idxft:idxlt] - max(Sync[idxft:idxlt])

    # theoretical time between two consecutive triggers
    t_btw_trig = numWFpfilepCH * acPeriod/1e6 # time in sec

    # build a trigger time vector using the large trigger chosen above as a reference

    # number of triggers between the two reference triggers:
    # this number should be an integer if no mismatch occurs between mechanical and acoustical recorders.
    # In the following, we take the biax recorder as the reference and adjust
    # the acoustical data to match the mechanical data.
    Nspans = (Time[idxref2]- Time[idxref1])/t_btw_trig

    # actual time between two consecutive triggers (using biax recorder as a reference)
    actual_t_btw_trig = (Time[idxref2] - Time[idxref1])/round(Nspans) # sec   (round instead of floor here...) format long
    mismatchPRF = (actual_t_btw_trig - t_btw_trig)*1e6/numWFpfilepCH
    acPeriod_adjusted = acPeriod + mismatchPRF # microsec
    acRate_adjusted = 1e6/(acPeriod_adjusted) # Hz

    print('In theory, acoustic pulses are sent every '+str(round(acPeriod/1e3,15))+' ms.\n' 
           'Using biax recorder as a reference, we find '+str(round(acPeriod_adjusted/1e3,15))+' ms.\n' 
           'Mismatch for the acoustic pulsing rate is '+str(mismatchPRF*1e3)+' ns.\n\n')

    mismatchts = (acRate/1e6)/fs * mismatchPRF; # microsec * MHz/MHz = microsec
    ts_adjusted = ts + mismatchts; # adjusted sampling period (microsec)
    fs_adjusted = 1/ts_adjusted; # adjusted sampling frequency (MHz)

    print('Said differently, sampling frequency is '+str(fs)+' MHz in theory.\n' 
          'Using biax recorder as a reference, it is adjusted to '+str(fs_adjusted)+' MHz.\n' 
          'Mismatch for sampling time is '+str(mismatchts*1e6)+' ps.\n\n')

    # idxref2 is used as a reference to build sync vectors
    # raw sync, using theoretical t_btw_trig
    trigger_time_begraw = np.flip(np.arange(Time[idxref2],Time[idxft],-1*t_btw_trig),axis=0)
    trigger_time_endraw = np.arange(Time[idxref2],Time[idxlt],t_btw_trig)


    # adjusted sync, using actual_t_btw_trig
    trigger_time_beg = np.flip(np.arange(Time[idxref2],Time[idxft],-1*actual_t_btw_trig),axis=0)
    trigger_time_end = np.arange(Time[idxref2],Time[idxlt],actual_t_btw_trig)

    # to make sure the real first trigger is not missed.
    ftmismatch = np.abs(trigger_time_beg[0] - Time[idxft])
    if ftmismatch > 0.5*actual_t_btw_trig: 
        trigger_time_beg = np.hstack([(trigger_time_beg[0] - actual_t_btw_trig), trigger_time_beg])
        print('First trigger added\n\n')

    # the sample corresponding to idxref2 is both in "beg" and "end" so we remove it from "beg".
    trigger_timeraw = np.hstack([trigger_time_begraw[0:-1], trigger_time_endraw]); # raw sync
    trigger_time = np.hstack([trigger_time_beg[0:-1], trigger_time_end]) # adjusted sync 

    check = np.diff(trigger_timeraw) # check time between triggers (should be constant...)
    check2 = np.diff(trigger_time) # check time between triggers (should be constant...)

    # this number should equal the number of acoustic files
    Ntriggerraw = len(trigger_timeraw)
    Ntrigger = len(trigger_time)

    TotalNumberOfFiles = Ntrigger
    print('The total number of acoustic files should be '+str(TotalNumberOfFiles)+'.\n\n')

    trigsraw = -0.5*np.ones(Ntriggerraw) # -0.5 to be adjusted depending on the experiment (adjust it to see both signals clearly on figure 2)
    trigs = -0.51*np.ones(Ntrigger)

    fig1 = []
    f1 = figure(title='Normalized Sync', tools='pan,box_zoom,undo,save,hover', plot_width=1200, plot_height=800) 
    f1.line(Time[idxft:idxlt], Sync[idxft:idxlt], line_width = 2, line_color='dimgray')
    f1.circle(trigger_timeraw, trigsraw, size = 10, fill_color='red', line_width = 0, legend_label='Raw Sync')
    f1.circle(trigger_time, trigs, size = 10, fill_color='gold', line_width = 0, legend_label='Adjusted Sync (using biax redorder as reference)')
    f1.square(Time[idxref1], -0.52, size = 10, fill_color='blue', line_width = 0, legend_label='1st Reference Trigger')
    f1.square(Time[idxref2], -0.52, size = 10, fill_color='navy', line_width = 0, legend_label='2nd Reference Trigger')
    f1.yaxis.axis_label = 'Normalized Sync'
    f1.xaxis.axis_label = 'Time (s)'
    fig1.append(f1)
    show(f1)

    
    # % Find peaks in the sync data (will work when recording rate is 1000Hz or higher)
    # % uncomment to see the experimental peaks...
    # % [pks,locs] = findpeaks(-Sync(idxft:idxlt),'MinPeakDistance',t_btw_trig*980); 
    # % % 980 is slightly lower than 1000Hz, the sampling rate used for this run
    # % locs = locs + idxft - 1;
    # % plot(Time(locs),Sync(locs),'sm');hold on

    fig2 = []
    f2 = figure(title='Time b/w consecutive triggers should be constant', tools='pan,box_zoom,undo,save,hover', plot_width=1200, plot_height=400) 
    f2.circle(np.arange(len(check)), check*1000, size = 10, fill_color='midnightblue', line_width = 0, legend_label='Expected')
    f2.circle(np.arange(len(check2)), check2*1000, size = 10, fill_color='dodgerblue', line_width = 0, legend_label='Adjusted (using biax recorder as reference)')
    f2.yaxis.axis_label = 'Time between consecutive triggers (ms)'
    f2.xaxis.axis_label = 'Index Number'
    fig2.append(f2)
    show(f2)

    print('Check that the time between consecutive triggers is constant, then zoom on the other figure to check that the sync is correct.\n\n')   

    
    if numAcFiles != TotalNumberOfFiles:
        print('Calculated number of files does not match number of .ac files in '+ runpath+'\n'
             'Make sure beggining and end sync triggers are correct')
    else: 
        print('Calculated number of files matches number of .ac files in '+ runpath+'\n'
             'You can now move on to futher analysis')
#     print('Check that the time between consecutive triggers is constant,\n' 
#           'then zoom on the other figure to check that the sync is correct.\n\n')                          

    acN = TotalNumberOfFiles * numWFpfilepCH; # total number of WF per channel

    # build acoustic time vector
    acTime = np.arange(0,acPeriod_adjusted*(acN),acPeriod_adjusted)/1e6; # seconds
    acTime = acTime + np.hstack(trigger_time)[0];                            # seconds

    # figure to check that 'acTime' vector is correct (it should overlap 'trigger_time' vector)
    # figure; 
    # plot(trigger_time,trigger_time,'*');hold on
    # plot(acTime(1:numWFpfilepCH:end),acTime(1:numWFpfilepCH:end),'s');
    # legend('trigger_time','acTime');
    # dcmObj = datacursormode;
    # set(dcmObj,'UpdateFcn',@GoodCursor); 
    
    return [acTime, acPeriod_adjusted, ts_adjusted, TotalNumberOfFiles]