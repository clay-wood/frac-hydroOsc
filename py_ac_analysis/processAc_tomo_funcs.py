import numpy as np
import scipy.io as sio
import psutil
import os
import fnmatch
from LoadAcFile import LoadAcFile
from findidxs import findidxs

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


def getfullWFref(NtoStack, numWFpSFpCH, WF_path, numSFpfile, WFlength, numCHR, numCHT):

    StackToNumWFs = NtoStack/(numWFpSFpCH*2)
    if StackToNumWFs <= 1:
        WFdat = LoadAcFile(WF_path,0,numCHR,numSFpfile)
        WFdat = np.reshape(WFdat,(numCHT,-1,WFlength,numCHR))
        WFdat = np.moveaxis(WFdat,-1,0)
        WFdat = np.moveaxis(WFdat,-1,0)
        fullWFref = np.sum(WFdat[:,:,:,0:NtoStack],axis = -1)
    else: 
        numRefFiles = np.ceil(StackToNumWFs).astype(int)
        WFdat = np.empty((numRefFiles,(numWFpSFpCH*2*8*WFlength),numCHR))
        for aa in range(numRefFiles):
            WFdat[aa,:,:] = LoadAcFile(WF_path,aa,numCHR,numSFpfile)
        WFdat = np.reshape(WFdat,(numCHT,-1,WFlength,numCHR))
        WFdat = np.reshape(WFdat,(numRefFiles,numCHT,-1,WFlength,numCHR))
        WFdat = np.moveaxis(WFdat,-1,2)
        WFdat = np.moveaxis(WFdat,3,4)
        WFdat = np.reshape(WFdat,(numCHT,numCHR,WFlength,-1))
        fullWFref = np.sum(WFdat[:,:,:,0:NtoStack],axis = -1)

    return fullWFref


def rms(X,ax):
    return np.sqrt(np.nanmean(X,axis=ax))


def amp(X):
    return (np.amax(X)-np.amin(X))


def folder_size(path):
    total = 0
    for entry in os.scandir(path):
        if entry.is_file():
            total += entry.stat().st_size
        elif entry.is_dir():
            total += folder_size(entry.path)
    return total


def calcWFchunks(WF_path):
    numAcFiles = len(fnmatch.filter(os.listdir(WF_path[0:-3]), '*.ac'))
    dir_size = folder_size(WF_path[0:-3])
    # int(dir_size / numAcFiles)
    WF1_size = os.path.getsize(WF_path[0:-3]+'WF_1.ac')
    print('There are '+np.str(numAcFiles)+' files in ' + str(WF_path[0:-3]))
    print('Size of Directory: '+str(np.round(dir_size/1e9,2))+' GB')

    avail_mem = psutil.virtual_memory().available
    WF_chunk = int(np.floor((0.75*avail_mem)/WF1_size))
    print('Available Memory: '+str(np.round(avail_mem / 1e9,2))+' GB')

    if dir_size > avail_mem:
        print('\nThere will be '+str(WF_chunk)+' *.ac files read into memory at a time.')
    else: 
        print('\nThe entire directory of *.ac files will be read into memory.')
    return numAcFiles, WF1_size, WF_chunk


def delay(X, ts, plot01):

    NN = (len(X)+1)/2
    tps = np.arange(-NN*ts,(NN+1)*ts,ts)
    [maxi,ind] = [np.amax(X),np.argmax(X)]
    [x1, x2, x3] = [ind-1, ind, ind+1]
    
    if x1 == 0:
        x1 = 1
        x2 = 2
        x3 = 3
    elif x3 == len(X):
        x1 = len(X) - 3 # arbitrarily choose the maximum at x2 = length(X) - 1;
        x2 = len(X) - 2
        x3 = len(X) - 1
        ind = len(X) - 1; 
        
    [y1 ,y2, y3] = X[[x1, x2, x3]]
    b = ((y1 - y3)*(x1**2 - x2**2) - (y1 - y2)*(x1**2 - x3**2))/((x1 - x3)*(x1**2 - x2**2) - (x1 - x2)*(x1**2 - x3**2))
    a = (y1 - y3 - b*(x1 - x3))/(x1**2 - x3**2)
    c = y1 - a*x1**2 - b*x1

    ind_max_corr = -b/(2*a)
    max_interpolation = a*ind_max_corr**2 +b*ind_max_corr + c

    timedelay = ts*(ind_max_corr - NN)
    interpX = np.arange(x1, x3+1, 0.01)
    interpY = a*interpX**2 + b*interpX + c
    interpT = ts*(interpX - NN+1)

    if plot01 == 1:
        fig1 = figure(title='WFs to Stack', tools='pan,box_zoom,undo,save,hover', plot_width=1200, plot_height=800) 
        fig1.circle(tps[ind], maxi, size = 10, fill_color='red', line_width = 0, legend_label='max of intercorr') 
        fig1.circle(timedelay, max_interpolation, size = 10, fill_color='limegreen', line_width = 0, legend_label='refined max of intercorr')
        fig1.circle(tps[NN], X[NN], size = 10, fill_color='blue', line_width = 0, legend_label='interpolation')
        fig1.circle(interpT, interpY, size = 10, fill_color='black', line_width = 0)
        fig1.yaxis.axis_label = 'Time between consecutive triggers (ms)'
        fig1.xaxis.axis_label = 'Index Number'
        show(fig1)
    
    return [max_interpolation, timedelay]