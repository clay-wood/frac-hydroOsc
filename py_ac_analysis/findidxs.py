import numpy as np

def findidxs(acTime,Time,TransNum,numCHT,numWFpfilepCH):

#     Use this function to find indexes corresponding to a certain Time within the run and corresponding to transmitter "TransNum".

#     Inputs:
#     acTime: Sync Time vector, output of SyncAcData.
#     Time: time at which we want the indexes (the closest indexes to that time will be found).
#     TransNum: The indexes provided corresponds to a waveform recorded when
#     Transmitter TransNum was active.
#     numCHT: number of transmitters used during the run
#     numWFpfilepCH: number of WF per receiver and per file

#     Outputs:
#     filenumber corresponding to the provided Time
#     idxWFwithinfile: idx of the WF of interest within filenumber
#     idxAcTime: idx of the WF within the acTime vector
#     ExactTime: exact time corresponding to the index found

    [mini,idxAcTime] = [np.amin(np.abs(acTime-Time)), np.argmin(np.abs(acTime-Time))] # find nearest idx within acTime

    ShiftTrans = np.mod(idxAcTime,numCHT);
    if ShiftTrans == 0:
        ShiftTrans = numCHT;  # i.e. no shift
    
    idxAcTime = idxAcTime + (TransNum - ShiftTrans); # closest WF within that file corresponding to transmitter TransNum

    ExactTime = acTime[idxAcTime];

    filenumber = np.ceil(idxAcTime/numWFpfilepCH); # corresponding file number
    idxWFwithinfile = np.mod(idxAcTime,numWFpfilepCH); # closest WF within the first file

    return [filenumber,idxWFwithinfile,idxAcTime,ExactTime]