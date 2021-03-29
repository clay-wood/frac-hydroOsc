import numpy as np
import pandas as pd
import scipy.io as sio


def DistanceSensorsVessel(sample_thickness,PZT_L,PZT_R,dy,dz):
    
    if len(dy) != len(dz):
        print('Vectors dy and dz must have same length.')
    #         break

    Nvec = len(dy); # length of the time series

    #     if size(dy,1) > size(dy,2) #
    #         dy = dy.T #size should be '1*length(vector)'
    #     if size(dz,1) > size(dz,2)  
    #         dz = dz.T # size should be '1*length(vector)'

    # steel_thickness = 4; #mm distance between piezos and the sample (steel thickness)
    steel_thickness_R = 4.1 # (mm) distance between piezos and the sample -- short block
    steel_thickness_L = 21.5367 # (mm) distance b/w piezos and the sample -- long block

    height_L = 68.453 #mm
    height_R = 50 #mm
    depth = 44.831 #mm
    NL = 12 # 12 sensors on the left
    NR = 9 # 9 sensors on the right

    d_sens_z = 15.494 # distance between sensors in the z direction (from bottom to top)
    d_sens_x = 12.7   # distance between sensors in the x direction (from front to back)

    d_edge_sens_x = 9.8;   # distance between front edge and 1st sensor in the x direction
    d_edge_sens_zL = 12.5768 # distance between bottom edge and 1st sensor in the z direction (left block)
    d_edge_sens_zR = 9.412 # distance between bottom edge and 1st sensor in the z direction (right block)

    # sample thickness time series (vector)
    sample_thickness_vec = sample_thickness - dy

    # # Size
    length_y = sample_thickness_vec + (steel_thickness_R + steel_thickness_L) # new thickness sample + 2*4 mm of steel at rest


    # ## Transducers positions (before shearing...)
    pos_L = np.zeros((3,NL,Nvec)) # xyz vs sensor number
    pos_R = np.zeros((3,NR,Nvec)) # xyz vs sensor number

    pos_L[0,0:3,:] = ((d_edge_sens_x + d_sens_x * np.arange(3)) * np.ones((3,Nvec)).T).T # large left block 
    pos_L[0,3:6,:] = ((d_edge_sens_x + d_sens_x * np.arange(3)) * np.ones((3,Nvec)).T).T
    pos_L[0,6:9,:] = ((d_edge_sens_x + d_sens_x * np.arange(3)) * np.ones((3,Nvec)).T).T
    pos_L[0,9:12,:] = ((d_edge_sens_x + d_sens_x * np.arange(3)) * np.ones((3,Nvec)).T).T
    pos_L[1,:,:] = (np.ones((NL,Nvec)) * length_y - steel_thickness_R)
    pos_L[2,[0, 3, 6, 9],:] = (d_edge_sens_zL + d_sens_z*np.arange(4) * np.ones((4,Nvec)).T - (np.ones((4,Nvec)) * dz).T).T
    pos_L[2,[1, 4, 7, 10],:] = (d_edge_sens_zL + d_sens_z*np.arange(4) * np.ones((4,Nvec)).T - (np.ones((4,Nvec)) * dz).T).T
    pos_L[2,[2, 5, 8, 11],:] = (d_edge_sens_zL + d_sens_z*np.arange(4) * np.ones((4,Nvec)).T - (np.ones((4,Nvec)) * dz).T).T

    pos_R[0,0:3,:] = ((d_edge_sens_x + d_sens_x * np.arange(3)) * np.ones((3,Nvec)).T).T # small right block
    pos_R[0,3:6,:] = ((d_edge_sens_x + d_sens_x * np.arange(3)) * np.ones((3,Nvec)).T).T
    pos_R[0,6:9,:] = ((d_edge_sens_x + d_sens_x * np.arange(3)) * np.ones((3,Nvec)).T).T
    pos_R[1,:,:] = -steel_thickness_R*np.ones((NR,Nvec)) # sample side is at 0
    pos_R[2,[0, 3, 6],:] = ((d_edge_sens_zR + d_sens_z * np.arange(3)) * np.ones((3,Nvec)).T).T
    pos_R[2,[1, 4, 7],:] = ((d_edge_sens_zR + d_sens_z * np.arange(3)) * np.ones((3,Nvec)).T).T
    pos_R[2,[2, 5, 8],:] = ((d_edge_sens_zR + d_sens_z * np.arange(3)) * np.ones((3,Nvec)).T).T

    distxz = np.full((NL, NR, Nvec), np.nan)
    disty = np.full((NL, NR, Nvec), np.nan)
    for R in range(NR):
        for L in range(NL):
            distxz[L,R,:] = np.sqrt((pos_R[0,R,:] - pos_L[0,L,:])**2 + (pos_R[2,R,:]-pos_L[2,L,:])**2)
            disty[L,R,:] = np.sqrt((pos_R[1,R,:] - pos_L[1,L,:])**2)  

    delta = np.arctan2(distxz,disty)
    dist_b = (steel_thickness_R + steel_thickness_L)/(np.cos(delta)); # dist_b is the distance of the ray in the block (distance between PZTs and the sample)
    # # dist_b = 2*steel_thickness./cosd(delta); # dist_b is the distance of the ray in the block (distance between PZTs and the sample)

    dist = np.full((NL, NR, Nvec), np.nan) # distance between pair of transducers
    for R in range(NR):
        for L in range(NL):
            dist[L,R,:] = np.sqrt((pos_R[0,R,:]-pos_L[0,L,:])**2 + (pos_R[1,R,:]-pos_L[1,L,:])**2 + (pos_R[2,R,:]-pos_L[2,L,:])**2);


    dist_s = dist - dist_b # distance within the sample

    pos_L = pos_L[:,PZT_L,:] # only keep PZTs that were used on the left
    pos_R = pos_R[:,PZT_R,:] # only keep PZTs that were used on the right
    dist_s = dist_s[PZT_L,PZT_R,:] # keep corresponding distances only
    dist_b = dist_b[PZT_L,PZT_R,:] # keep corresponding distances only

    return pos_L,pos_R,dist_s,dist_b