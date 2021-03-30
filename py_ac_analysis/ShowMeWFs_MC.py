import numpy as np
import scipy.io as sio
from findidxs import *
from LoadAcFile import *
import matplotlib.pyplot as plt

from bokeh.plotting import figure, show
from bokeh.io import output_notebook
from bokeh.layouts import gridplot, row, column
output_notebook()

def ShowMeWFs_MC(WF_path,AcSettingsfile,SyncFile,AtWhichTimes,NtoStack,Offset,WhichTransTomoMode,IdxWindow):

    # This function is used to display waveforms at particular times during the run, 
    # with the options of stacking them. It is typically used to pickarrival times by hand, 
    # or decide how many waveforms need to be stacked before further processing

    # acoustic parameters
    acSettings = sio.loadmat(AcSettingsfile); # load acoustic settings
    numSFpfile = acSettings['numFrames'][0,0]/2; # number of superframes per file
    numWFpSFpCH = acSettings['numAcqs'][0,0]; # number of WF per superframe and per channel
    numWFpfilepCH = numSFpfile*numWFpSFpCH; # number of WF per file and per channel
    numCHR = np.size(acSettings['channels2save'][0]); # number of channels
    numCHT = np.size(acSettings['channels2transmit'][0]); # number of channels
    WFlength = acSettings['Nsamples'][0,0]; # waveform length
    del acSettings


    # MAYBE THIS DOES NOT NEED TO BE IMPLEMENTED AT THE MOMENT
    # if ~isequal(size(IdxWindow),[2 numCHR numCHT]) && ~isequal(size(IdxWindow),[0 0]) && ~isequal(size(IdxWindow),[1 2]) && ~isequal(size(IdxWindow),[2 1])
    #     display(['Size of IdxWindow is not correct. Size should either be [2 by ' num2str(numCHR) ' by ' num2str(numCHT) '], or empty, or a two element vector.']);

    # # if only a two element vector is provided, use it for all combinations of TR
    # if isequal(size(IdxWindow),[1 2]) || isequal(size(IdxWindow),[2 1])
    #     IdxWindowInit = IdxWindow;
    #     IdxWindow = repmat(IdxWindowInit,1,numCHR,numCHT);
    # end


    # Load sync data
    [acTime, acPeriod_adjusted, ts, TotalNumberOfFiles] = list(map(lambda x: np.load(SyncFile)[x], ['arr_0','arr_1','arr_2','arr_3']))
    fs = 1/ts; # acoustic sampling rate
    
    if WhichTransTomoMode > numCHT:    
        print('You chose to display transmitter #'+str(WhichTransTomoMode)+'. Please choose a transmitter between 1 and '+str(numCHT)+'.');
    # IT WOULD BE NICE TO RETURN AN ERROR MESSAGE HERE    

    # time vector for each waveform
    timeWF = np.arange(0,WFlength)*ts;

    # number of WFs to show
    N = np.size(AtWhichTimes)
    rangeTimes = np.zeros((N,2))

    p = figure(title='WFs to Stack', tools='pan,box_zoom,undo,save,hover', plot_width=1200, plot_height=800) # initialize plot for WFs

    for ii in range(0,N):
        idxAcTime = np.where(acTime > AtWhichTimes[ii])[0][0];      
        idxAcTimeVec = np.arange(idxAcTime-np.ceil(NtoStack/2)+1, idxAcTime+np.floor(NtoStack/2), dtype=int); # vector of indexes centered around 'idxAcTime'    
        if idxAcTimeVec[0] < 1:
            print('Either the first value of ''AtWhichTimes'' is too small or ''NtoStack'' is too large. Can not stack '+np.str(NtoStack)+' waveforms centered around '+np.str(AtWhichTimes[0])+ ' s.');
        elif idxAcTimeVec[-1] > np.size(acTime):
            print('Either the last value of ''AtWhichTimes'' is too large or ''NtoStack'' is too large. Can not stack '+np.str(NtoStack)+' waveforms centered around '+np.str(AtWhichTimes[-1])+' s.');

        rangeTimes[ii,:] = np.array([acTime[idxAcTimeVec[0]], acTime[idxAcTimeVec[-1]]]);     

        [filenumber,idxWFwithinfile,idxExactTime,ExactTime] = findidxs(acTime,rangeTimes[ii,1],WhichTransTomoMode,numCHT,numWFpfilepCH);

        fullWFref = np.zeros([WFlength,numCHR]);
        for kk in range(0,NtoStack): # number of WFs to stack                                

            if np.logical_or(idxWFwithinfile <= numCHT, kk == 0): # open new file if idxWFwithinfile is 1 or if it's the first WF to stack
                ACdata = LoadAcFile(WF_path,int(filenumber-1),numCHR,numSFpfile);

            fullWFref = fullWFref + ACdata[int(WFlength*(idxWFwithinfile-1)):int(WFlength*idxWFwithinfile),:]; # stack WFs

            if idxWFwithinfile <= numWFpfilepCH - numCHT: # stay within the same file for the next loop
                idxWFwithinfile = idxWFwithinfile + numCHT; # update to the next WF corresponding to the same transmitter 
            else:                                # use next file for the next loop
                filenumber = filenumber + 1; # go to next file
                idxWFwithinfile = WhichTransTomoMode; # start in the next file at WFs corresponding to transmitter "WhichTransTomoMode"

        fullWFref = fullWFref/NtoStack;    

        colors = np.array(['red','orange','blue'])
        if N <= 3:
            for chnum in range(numCHR):
                p.line(timeWF, fullWFref[:,chnum]-Offset*(chnum), line_color=colors[ii], line_width = 3, legend_label='WF_'+str(ii+1))
        else:
            for chnum in range(numCHR):
                p.line(timeWF, fullWFref[:,chnum]-Offset*(chnum), line_width = 3, legend_label='WF_'+str(ii+1))
    p.yaxis.axis_label = 'Amplitude'
    p.xaxis.axis_label = 'Time (us)'        
    p.legend.click_policy='hide'
    show(p)