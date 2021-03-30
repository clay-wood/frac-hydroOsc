import numpy as np

def LoadAcFile(WF_path,filenumber,numCHR,numSFpfile):

#     Load acoustic file and reshape according to the
#     number of receivers
#     WF_path: location of the files
#     filenumber: which file is being loaded
#     numCHR: number of receivers 
#     numSFpfile: number of 'superframes' per file (Verasonics jargon)

    with open(WF_path+np.str(int(filenumber)+1)+'.ac', 'rb') as fid:
            
        ACdata = np.array(np.fromfile(fid, np.int16)) 
        
        # reshape to get one column per channel
        ACdata = np.reshape(ACdata,(int(numSFpfile),numCHR,-1)) # 3D matrix with WF vs numCHR vs number of SF
        ACdata = np.rollaxis(ACdata,1,0) # put numCHR as the last dimension before reshaping
        ACdata = ACdata.reshape((numCHR,-1)) # WF vs numCHRs
        
    fid.close();
    
    return ACdata.T