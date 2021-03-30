import numpy as np
import scipy.io as sio
from LoadAcFile import *


def rms(x):
    rms = np.sqrt(np.mean(y**2))
    return rms


def WF_ref(fullWFref, ):

    fullWFref = fullWFref/NtoStack/reference_NrefWF; 

    # -----------------------------------------------------------------------------
    # IMPLEMENT FILTERING IN PYTHON
    # if Filter == 'yes':
    #     fullWFrefF = filtfilt(filterparam,1,fullWFref); # filtering
    #     if Filter_view == 0:
    #         figure(765);plot(fullWFref(:,6,5));hold on; # plot pair R1-T1 unfiltered
    #         plot(fullWFrefF(:,6,5),'k');hold off; # plot pair R1-T1 filtered

    #     fullWFref = fullWFrefF;
    # -----------------------------------------------------------------------------

    # if only a two element vector is provided, use it for all combinations of TR
    if size(IdxWindow),[1 2] or size(IdxWindow),[2 1]:
        IdxWindowInit = IdxWindow;
        IdxWindow = repmat(IdxWindowInit,1,numCHR,numCHT);

    # isolate windows to be analyzed
    windows_length = IdxWindow[1,:,:] - IdxWindow[0,:,:] + 1;
    maxwindowlength = np.max[windows_length(:)];

    WFref = np.zeros((maxwindowlength,numCHR,numCHT));
    for chnumt in range(0,numCHT):
        for chnumr in range(0,numCHR):        
            WFref[0:windows_length[0,chnumr,chnumt],chnumr,chnumt] = fullWFref[IdxWindow[0,chnumr,chnumt]:IdxWindow[1,chnumr,chnumt],chnumr,chnumt]; # part of the WF to be analyzed

    for chnumr in range(0,numCHR):
        RmsAmpRef[chnumr,:] = rms(WFref[:,chnumr,:]); # RmsAmp of the reference waveform
        AmpRef[chnumr,:] = np.max(WFref[:,chnumr,:])-np.min(WFref[:,chnumr,:]); # Peak-to-Peak Amp of the reference waveform

    #matlab code creates plot

    del ACdata
return [WFref, windows_length, maxwindowlength, RmsAmpRef, AmpRef]