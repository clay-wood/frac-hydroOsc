import numpy as np


def delay(X,ts):

    NN = (len(X)+1)/2;

    tps = np.arange(-NN*ts,(NN+1)*ts,ts);
    
    [maxi,ind] = [np.amax(X),np.argmax(X)];
    
    

    [x1, x2, x3] = [ind-1, ind, ind+1]
    
    if x1 == 0:
        x1 = 1; 
        x2 = 2;
        x3 = 3;
    elif x3 == len(X):
        x1 = len(X) - 3; # arbitrarily choose the maximum at x2 = length(X) - 1;
        x2 = len(X) - 2
        x3 = len(X) - 1
        ind = len(X) - 1; 
        

    [y1 ,y2, y3] = X[[x1, x2, x3]]

    b = ((y1 - y3)*(x1**2 - x2**2) - (y1 - y2)*(x1**2 - x3**2))/((x1 - x3)*(x1**2 - x2**2) - (x1 - x2)*(x1**2 - x3**2));

    a = (y1 - y3 - b*(x1 - x3))/(x1**2 - x3**2);

    c = y1 - a*x1**2 - b*x1;

    ind_max_corr = -b/(2*a);
    max_interpolation = a*ind_max_corr**2 +b*ind_max_corr + c;

    timedelay = ts*(ind_max_corr - NN);

    interpolationX = np.arange(x1,x3+1,0.01);

    interpolationY = a*interpolationX**2 + b*interpolationX + c;

    interpolationt = ts*(interpolationX - NN+1);
    
    return [max_interpolation, timedelay]



