import numpy as np

def projectionProcedure(t, V_HF, strain, Nharm, fn):
    
    offset = np.mean(V_HF)
    offset_strain = np.mean(strain)
    dCoverC = V_HF-offset
    dCoverC_strain = strain-offset_strain
    
    T = np.zeros((len(t),2*Nharm)) # matrix containing the basis of orthogonal functions
        
    # make functions
    wLF = 2*pi*fn; # Freq will be the frequency of the pump
    
    for iii in range(Nharm):
        T[:,2*iii-1] = np.sin(iii*wLF*t); # rows of sines and cosines
        T[:,2*iii] = np.cos(iii*wLF*t);
    

    # this loop ensures orthogonality
    for iii in range(2*Nharm):
        if iii > 1:
            for jjj in range(iii-1):
                T[:,iii] = T[:,iii] - np.dot(T[:,jjj],T[:,iii])*T[:,jjj]; # remove cross terms

            T[:,iii] = T[:,iii]/np.sqrt(np.dot(T[:,iii],T[:,iii])); # for all but the first
        else:
            T[:,iii] = T[:,iii]/np.sqrt(np.dot(T[:,iii],T[:,iii])); # for the first one

            
    # make scalar products to control orthogonality
    orthotest = np.zeros((2*Nharm,2*Nharm));
    for iii in range(2*Nharm):
        for jjj in range(iii,2*Nharm):
            orthotest[jjj,iii] = np.dot(T[:,iii],T[:,jjj]);

    print(orthotest) # we should get 1 on the diagonal, 0 otherwise
    del orthotest
    
    
    projection = 0
    projection_strain = 0

    #EVERYTHING BELOW NEEDS TRANSLATED
    
    A = np.zeros(size(T,2),1); # coefficients
    A_strain = np.zeros(size(T,2),1); # coefficients
    # norm f to 1
    # dCoverC=dCoverC/sqrt(dot(dCoverC,dCoverC));
    for jjj = 1:size(T,2)    
        A(jjj) = sum(T(:,jjj).*dCoverC); # dCoverC is the signal you want to analyze (dC/C)
        projection = projection + A(jjj)*T(:,jjj);   
        A_strain(jjj) = sum(T(:,jjj).*dCoverC_strain); # dCoverC is the signal you want to analyze (dC/C)
        projection_strain = projection_strain + A_strain(jjj)*T(:,jjj);   
    end
    
    
    # figure(11);
    # plot(t,dCoverC+mean(V_HF),'-ok');hold on;
    # plot(t,projection+mean(V_HF),'-r');

    for iii = 1:Nharm
        R(iii) = sqrt((A(2*iii-1)*max(T(:,2*iii-1)))^2 + (A(2*iii)*max(T(:,2*iii)))^2);
        R_strain(iii) = sqrt((A_strain(2*iii-1)*max(T(:,2*iii-1)))^2 + (A_strain(2*iii)*max(T(:,2*iii)))^2);
    end

    R0 = offset;
    R0_strain = offset_strain;

    # Nonlinear parameter
    beta = 2*R(1)/((max(projection_strain)-min(projection_strain))*1e-6);
    alpha = 2*R0/((max(projection_strain)-min(projection_strain))*1e-6);
    delta = 2*R(2)/((max(projection_strain)-min(projection_strain))*1e-6)^2;

    # projection signal
    signal_projection = 0;
    strain_projection = 0;

    for iii = 1:Nharm
        signal_projection = signal_projection + A(2*iii-1)*sin(iii*wLF*t)*max(T(:,2*iii-1))+A(2*iii)*cos(iii*wLF*t)*max(T(:,2*iii));
        strain_projection = strain_projection + A_strain(2*iii-1)*sin(iii*wLF*t)*max(T(:,2*iii-1))+A_strain(2*iii)*cos(iii*wLF*t)*max(T(:,2*iii));
    end

    R1 = R(1);
    R2 = R(2);

    # 
    # figure(44)
    # plot(strain_projection,attenuation_projection,'.k')
    
    return [strain_projection, signal_projection, alpha, beta, delta, R0, R1, R2]