import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from plotter_functions import *
from legend_specs import *
#----------------------------------------------------------------
#USEFUL THINGS & FUNCTIONS:

filepath = '/Users/clay/Documents/research/jgr20/manuscript/new_plots/'

[legend_elements1N, legend_elements1P, legend_elements2N, legend_elements2P, legend_elements2Q, legend_elements2R] = legend_stuff()

#----------------------------------------------------------------
#----------------------------------------------------------------

def avgDelc_All_plot(x, y, yerr, osc_type, xlim = [-15, 20], ylim = [-0.2, 0.01], save = 'n'):
# def avgDelc_All_plot(x, y, yerr, osc_type, xlim, ylim, save):

    [x1, x2, x3, x4, x5] = list(map(lambda aa: x[aa], np.arange(5)))
    if osc_type == 'PP':
        [y1, y2, y3, y4, y5] = list(map(lambda aa: tenHz_filter(x[aa],y[aa]), np.arange(5)))
    else:
        [y1, y2, y3, y4, y5] = list(map(lambda aa: y[aa], np.arange(5)))
        
    [yerr1, yerr2, yerr3, yerr4, yerr5] = list(map(lambda aa: yerr[aa], np.arange(5)))

    plt.style.use('/Users/clay/Documents/research/jgr20/summary_data/jgr2020_style2.mplstyle')
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw={'wspace': 0.1})
    
    clr1, kdot, amp, perm = which_osc(osc_type,x1)
    for aa in range(len(y1)):
        ax1.errorbar(x1[aa,perm], y1[aa], yerr1[aa], fmt = marker_shape(x1[aa,1]), markersize=marker_size(x1[aa,1]), color=clr1, mfc=clr1)
    ax1.set_title('$Post-fracture$')
        
    clr2, kdot, amp, perm = which_osc(osc_type,x2)
    for bb in range(len(y2)): 
        ax1.errorbar(x2[bb,perm], y2[bb], yerr2[bb], fmt = marker_shape(x2[bb,1]), markersize=marker_size(x2[bb,1]), color=clr2, mfc=clr2)
    
    clr3, kdot, amp, perm = which_osc(osc_type,x3)
    for aa in range(len(y3)):
        ax2.errorbar(x3[aa,perm], y3[aa], yerr3[aa], fmt = marker_shape(x3[aa,1]), markersize=marker_size(x3[aa,1]), color=clr3, mfc=clr3) 
    clr4, kdot, amp, perm = which_osc(osc_type,x4)
    for bb in range(len(y4)): 
        ax2.errorbar(x4[bb,perm], y4[bb], yerr4[bb], fmt = marker_shape(x4[bb,1]), markersize=marker_size(x4[bb,1]), color=clr4, mfc=clr4)
    ax2.set_title('$Post-shear\ 1$')

    clr5, kdot, amp, perm = which_osc(osc_type,x5)
    for aa in range(len(y5)):
        ax3.errorbar(x5[aa,perm], y5[aa], yerr5[aa], fmt = marker_shape(x5[aa,1]), markersize=marker_size(x5[aa,1]), color=clr5, mfc=clr5)
    ax3.set_title('$Post-shear\ 2$')
        
    #ADD LEGEND TO RIGHT OF LAST PLOT

    if osc_type=='NS':
        legend1 = plt.legend(handles = legend_elements1N, title='$\sigma_{NS}\ Exp.\ \#$', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        legend2 = plt.legend(handles = legend_elements2N, title='$Osc.\ Freq.$', bbox_to_anchor=(1.05, 0.65), loc='upper left', borderaxespad=0.)
    else: 
        legend1 = plt.legend(handles = legend_elements1P, title='$P_P\ Exp.\ \#$', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        legend2 = plt.legend(handles = legend_elements2P, title='$Osc.\ Freq.$', bbox_to_anchor=(1.05, 0.65), loc='upper left', borderaxespad=0.)
       
    ax3.add_artist(legend1); ax3.add_artist(legend2)
        
    #PLOT 0 LINES
    for ax in (ax1,ax2,ax3):
        ax.plot([0, 0], [-0.2, 0.1], color = 'gray',ls = '--')
        ax.plot([-15, 20], [0, 0], color = 'gray',ls = '--')
     
    # AXES: MINOR TICKS, COLORS 
    for ax in (ax1,ax2,ax3):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('$\Delta k/k_0$ (%)',color=[0,0,0])
    ax2.tick_params(labelleft=False)
    ax3.tick_params(labelleft=False)
 
    # AXES LABELS & COLORS   
    ax1.ticklabel_format(axis='y',style='sci',scilimits=(-4,2))
    ax1.set_ylabel('$\overline{\Delta c/c_0}$ (%)',color=[0,0,0])
    plt.show()

    if save == 'y': 
#         fig.savefig(filepath+'avgDelc_All'+osc_type+'.png')
        fig.savefig(filepath+'avgDelc_All'+osc_type+'.svg')        

#----------------------------------------------------------------
#----------------------------------------------------------------

def avgDelc_All_plot2(x, y, yerr, osc_type, xlim = [-15, 20], ylim = [-0.2, 0.01], ytick = np.arange(0,-0.14,0.04), save = 'n'):
    mpl.rcParams.update({'font.size': 12})
# def avgDelc_All_plot(x, y, yerr, osc_type, xlim, ylim, save):

    [x1, x2, x3, x4, x5] = list(map(lambda aa: x[aa], np.arange(5)))
    if osc_type == 'PP':
        [y1, y2, y3, y4, y5] = list(map(lambda aa: tenHz_filter(x[aa],y[aa]), np.arange(5)))
    else:
        [y1, y2, y3, y4, y5] = list(map(lambda aa: y[aa], np.arange(5)))
        
    [yerr1, yerr2, yerr3, yerr4, yerr5] = list(map(lambda aa: yerr[aa], np.arange(5)))

    plt.style.use('/Users/clay/Documents/research/jgr20/summary_data/jgr2020_style2.mplstyle')
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw={'wspace': 0.1})
    
    clr1, kdot, amp, perm = which_osc(osc_type,x1)
    for aa in range(len(y1)):
        ax1.errorbar(x1[aa,perm], y1[aa], yerr1[aa], fmt = marker_shape(x1[aa,1]), markersize=mrk_amp_sz(x1[aa,2]), color=clr1, mfc=clr1)
    ax1.set_title('$Post-fracture$')
        
    clr2, kdot, amp, perm = which_osc(osc_type,x2)
    for bb in range(len(y2)): 
        ax1.errorbar(x2[bb,perm], y2[bb], yerr2[bb], fmt = marker_shape(x2[bb,1]), markersize=mrk_amp_sz(x2[bb,2]), color=clr2, mfc=clr2)
    
    clr3, kdot, amp, perm = which_osc(osc_type,x3)
    for aa in range(len(y3)):
        ax2.errorbar(x3[aa,perm], y3[aa], yerr3[aa], fmt = marker_shape(x3[aa,1]), markersize=mrk_amp_sz(x3[aa,2]), color=clr3, mfc=clr3) 
    clr4, kdot, amp, perm = which_osc(osc_type,x4)
    for bb in range(len(y4)): 
        ax2.errorbar(x4[bb,perm], y4[bb], yerr4[bb], fmt = marker_shape(x4[bb,1]), markersize=mrk_amp_sz(x4[bb,2]), color=clr4, mfc=clr4)
    ax2.set_title('$Post-shear\ 1$')

    clr5, kdot, amp, perm = which_osc(osc_type,x5)
    for aa in range(len(y5)):
        ax3.errorbar(x5[aa,perm], y5[aa], yerr5[aa], fmt = marker_shape(x5[aa,1]), markersize=mrk_amp_sz(x5[aa,2]), color=clr5, mfc=clr5)
    ax3.set_title('$Post-shear\ 2$')
        
    #ADD LEGEND TO RIGHT OF LAST PLOT

    if osc_type=='NS':
        legend1 = plt.legend(handles = legend_elements1N, title='$\sigma_{NS}\ Exp.\ \#$', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        legend2 = plt.legend(handles = legend_elements2N, title='$Osc.\ Freq.$', bbox_to_anchor=(1.05, 0.65), loc='upper left', borderaxespad=0.)
    else: 
        legend1 = plt.legend(handles = legend_elements1P, title='$P_P\ Exp.\ \#$', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        legend2 = plt.legend(handles = legend_elements2P, title='$Osc.\ Freq.$', bbox_to_anchor=(1.05, 0.65), loc='upper left', borderaxespad=0.)
       
    ax3.add_artist(legend1); ax3.add_artist(legend2)
        
    #PLOT 0 LINES
    for ax in (ax1,ax2,ax3):
        ax.plot([0, 0], [-0.2, 0.1], color = 'gray',ls = '--')
        ax.plot([-15, 20], [0, 0], color = 'gray',ls = '--')
        ax.set_yticks(ytick)
     
    # AXES: MINOR TICKS, COLORS 
    for ax in (ax1,ax2,ax3):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('$\Delta k/k_0$ (%)',color=[0,0,0])
    ax2.tick_params(labelleft=False)
    ax3.tick_params(labelleft=False)
 
    # AXES LABELS & COLORS   
    ax1.ticklabel_format(axis='y',style='sci',scilimits=(-4,2))
    ax1.set_ylabel('$\overline{\Delta c/c_0}$ (%)',color=[0,0,0])
    plt.show()

    if save == 'y': 
#         fig.savefig(filepath+'avgDelc_All'+osc_type+'.png')
        fig.savefig(filepath+'avgDelc_All_amps'+osc_type+'.svg')
        fig.savefig(filepath+'avgDelc_All_amps'+osc_type+'.pdf')

#----------------------------------------------------------------
#----------------------------------------------------------------

# def Delc_plot(x1, osc_type, xlim, ylim, save, name):

#     plt.style.use('/home/clay/Documents/research/jgr20/summary_data/jgr2020_style.mplstyle')
#     fig, ax = plt.subplots()

#     clr, kdot, amp, perm = which_osc(osc_type,x1)
#     for aa in range(len(x1)):
#         ax.errorbar(x1[aa,perm], x1[aa,3], yerr=None, fmt = 'o', markersize='6', color=clr, mfc=clr) 

#     # AXES: MINOR TICKS, COLORS 
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
    
#     # AXES LABELS & COLORS   
#     ax.set_xlabel('$\Delta k/k_0\ (\%)$',color=[0,0,0])
#     ax.ticklabel_format(axis='y',style='sci',scilimits=(-4,2))
#     ax.set_ylabel('$\Delta c/c_0 (\%)$',color=[0,0,0])
#     plt.show()

#     if save == 'y': 
#         fig.savefig(filepath+'Delc_'+osc_type+'_'+name+'.svg')
#         fig.savefig(filepath+'Delc_'+osc_type+'_'+name+'.png')


# #----------------------------------------------------------------


# def Delc_Both_plot(x1, x2, osc_type, xlim, ylim, save, name):
#     #order is: p4975 f(x1, y1, y1err), p4966 (x2, y2, y2err)

#     plt.style.use('/home/clay/Documents/research/jgr20/summary_data/jgr2020_style.mplstyle')
#     fig, ax = plt.subplots()
    
#     clr1, kdot, amp, perm = which_osc(osc_type,x1)
#     for aa in range(len(x1)):
#         ax.errorbar(x1[aa,perm], x1[aa,3], yerror=None, fmt = 'o', markersize='6', color=clr1, mfc=clr1) 

#     clr2, kdot, amp, perm = which_osc(osc_type,x2)
#     for bb in range(len(x2)): 
#         ax.errorbar(x2[bb,perm], x2[bb,3], yerror=None, fmt = 'o', markersize='6', color=clr2, mfc=clr2) 

#     # AXES: MINOR TICKS, COLORS 
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
 
#     # AXES LABELS & COLORS   
#     ax.set_xlabel('$\Delta k/k_0\ (\%)$',color=[0,0,0])
#     ax.ticklabel_format(axis='y',style='sci',scilimits=(-4,2))
#     ax.set_ylabel('$\Delta c/c_0 (\%)$',color=[0,0,0])
#     plt.show()

#     if save == 'y': 
#         fig.savefig(filepath+'Delc_'+osc_type+'_'+name+'.svg')
#         fig.savefig(filepath+'Delc_'+osc_type+'_'+name+'.png')


#----------------------------------------------------------------
   
def delc_amp_plots(x, y, osc_type, ylim, yerr = 'na', xlim = None, save = 'n'):
    
    mpl.rcParams.update({'font.size': 12})
    
    [x1, x2, x3, x4, x5] = list(map(lambda aa: oneHz_filter(x[aa]), np.arange(5)))
    [y1, y2, y3, y4, y5] = list(map(lambda aa: y[aa], np.arange(5)))
    runs_75 = np.array([1,1,1,1,1,1,1,2,2,2,2,2,2,2,None,None,None,None,3,3,3,3,None,None,None,None,None,None,None,None,None,4,None]);
    runs_66 = np.array([1,1,1,1,2,2,2,2,None,None,None,None,3,3,3,3,None,None,None,None]);

        
    if yerr == 'na':
        [yerr1, yerr2, yerr3, yerr4, yerr5] = [None,None,None,None,None]
    else: 
        [yerr1, yerr2, yerr3, yerr4, yerr5] = list(map(lambda aa: yerr[aa], np.arange(5)))

    plt.style.use('/Users/clay/Documents/research/jgr20/summary_data/jgr2020_style2.mplstyle')
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw={'wspace': 0.1})
    
    clr1, kdot, amp, perm = which_osc(osc_type,x1)
    for aa in range(len(y1)):
        ax1.errorbar(x1[aa,amp], y1[aa], yerr=yerr1, fmt = mrkr_run(runs_75[aa])[0], markersize=mrkr_run(runs_75[aa])[1], color=clr1, mfc=clr1)
    ax1.set_title('$Post-fracture$')
        
    clr2, kdot, amp, perm = which_osc(osc_type,x2)
    for bb in range(len(y2)): 
        ax1.errorbar(x2[bb,amp], y2[bb], yerr=yerr2, fmt = mrkr_run(runs_66[bb])[0], markersize=mrkr_run(runs_66[bb])[1], color=clr2, mfc=clr2)
    
    clr3, kdot, amp, perm = which_osc(osc_type,x3)
    for aa in range(len(y3)):
        ax2.errorbar(x3[aa,amp], y3[aa], yerr=yerr3, fmt = mrkr_run(runs_75[aa])[0], markersize=mrkr_run(runs_75[aa])[1], color=clr3, mfc=clr3) 
    clr4, kdot, amp, perm = which_osc(osc_type,x4)
    for bb in range(len(y4)): 
        ax2.errorbar(x4[bb,amp], y4[bb], yerr=yerr4, fmt = mrkr_run(runs_66[bb])[0], markersize=mrkr_run(runs_66[bb])[1], color=clr4, mfc=clr4)
    ax2.set_title('$Post-shear\ 1$')

    clr5, kdot, amp, perm = which_osc(osc_type,x5)
    for aa in range(len(y5)):
        ax3.errorbar(x5[aa,amp], y5[aa], yerr=yerr5, fmt = mrkr_run(runs_75[aa])[0], markersize=mrkr_run(runs_75[aa])[1], color=clr5, mfc=clr5)
    ax3.set_title('$Post-shear\ 2$')
        
    #ADD LEGEND TO RIGHT OF LAST PLOT
    if osc_type=='NS':
        legend1 = plt.legend(handles = legend_elements1N, title='$\sigma_{NS}\ Exp.\ \#$', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        legend2 = plt.legend(handles = legend_elements2Q, title='$Osc.\ Order$', bbox_to_anchor=(1.05, 0.65), loc='upper left', borderaxespad=0.)
    else: 
        legend1 = plt.legend(handles = legend_elements1P, title='$P_P\ Exp.\ \#$', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        legend2 = plt.legend(handles = legend_elements2R, title='$Osc.\ Order$', bbox_to_anchor=(1.05, 0.65), loc='upper left', borderaxespad=0.)
       
    ax3.add_artist(legend1); ax3.add_artist(legend2)
    
    for ax in (ax1,ax2,ax3):
        ax.plot([0.1, 1.5], [0, 0], color = 'gray',ls = '--') #Horizontal Line
    
    # AXES: MINOR TICKS, COLORS 
    for ax in (ax1,ax2,ax3):
        ax.set_xlim([0.15, 1.3])
        ax.set_ylim(ylim)
        ax.set_xlabel('$Oscillation\ Amp.\ (MPa)$',color=[0,0,0])
        ax.set_xticks(np.arange(0.2, 1.4, 0.2))
    ax2.tick_params(labelleft=False)
    ax3.tick_params(labelleft=False)
 
    # AXES LABELS & COLORS   
    ax1.ticklabel_format(axis='y',style='sci',scilimits=(-4,2))
    ax1.set_ylabel('$\Delta c/c_0\ (\%)$',color=[0,0,0])
    plt.show()

    if save == 'y': 
        fig.savefig(filepath+'delc_amp_'+osc_type+'.svg')
        fig.savefig(filepath+'delc_amp_'+osc_type+'.pdf')


#----------------------------------------------------------------


def freq_plot(x, osc_type, xlim, ylim, save, name):

    [x1, x2, x3, x4, x5] = list(map(lambda aa: oneHz_filter(x[aa]), np.arange(5)))
   
    plt.style.use('/Users/clay/Documents/research/jgr20/summary_data/jgr2020_style.mplstyle')
    fig, ax = plt.subplots()

    freq = 1; 

    clr, kdot, amp, perm = which_osc(osc_type,x1)
    for aa in range(len(x1)):
        ax.errorbar(x1[aa,freq], x1[aa,3], yerror=None, fmt = 'd', markersize=6, color = clr, mfc = clr) 

    clr, kdot, amp, permq = which_osc(osc_type,x2)
    for aa in range(len(x2)):
        ax.errorbar(x2[aa,freq], x2[aa,3], yerror=None, fmt = '*', markersize=9, color = clr, mfc = clr) 

    clr, kdot, amp, perm = which_osc(osc_type,x3)
    for aa in range(len(x3)):
        ax.errorbar(x3[aa,freq], x3[aa,3], yerror=None, fmt = 's', markersize=6, color = clr, mfc = clr) 

    clr, kdot, amp, perm = which_osc(osc_type,x4)
    for aa in range(len(x4)):
        ax.errorbar(x4[aa,freq], x4[aa,3], yerror=None, fmt = 'd', markersize=6, color = clr, mfc = clr) 

    clr, kdot, amp, perm = which_osc(osc_type,x5)
    for aa in range(len(x5)):
        ax.errorbar(x5[aa,freq], x5[aa,3], yerror=None, fmt = '*', markersize=9, color = clr, mfc = clr)

    # AXES: MINOR TICKS, COLORS 
    # ax.set_xlim(xlim)
    ax.set_xscale('log')
    ax.set_ylim(ylim)
    ax.minorticks_on()

    # AXES LABELS & COLORS   
    ax.set_xlabel('Oscillation Freq. (Hz)',color=[0,0,0])
    ax.ticklabel_format(axis='y',style='sci',scilimits=(-4,2))
    if name == 'Delc': 
        ax.set_ylabel('$\Delta c/c_0\ (\%)$',color=[0,0,0])
    else: 
        ax.set_ylabel('$c\ recovery$',color=[0,0,0])
    plt.show()

    if save == 'y': 
        fig.savefig(filepath+name+'_Freq_'+osc_type+'_All.svg')
        fig.savefig(filepath+name+'_Freq_'+osc_type+'_All.png')
        
        
        
        
        
def Delc_All_plots(x, y, yerr, osc_type, xlim = [-15, 20], ylim = [-0.2, 0.01], save = 'n'):
# def avgDelc_All_plot(x, y, yerr, osc_type, xlim, ylim, save):

    [x1, x2, x3, x4, x5] = list(map(lambda aa: x[aa], np.arange(5)))
    if osc_type == 'PP':
        [y1, y2, y3, y4, y5] = list(map(lambda aa: tenHz_filter(x[aa],y[aa]), np.arange(5)))
    else:
        [y1, y2, y3, y4, y5] = list(map(lambda aa: y[aa], np.arange(5)))
        
    [yerr1, yerr2, yerr3, yerr4, yerr5] = list(map(lambda aa: yerr[aa], np.arange(5)))

    plt.style.use('/Users/clay/Documents/research/jgr20/summary_data/jgr2020_style2.mplstyle')
    
    for cc in range(4):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw={'wspace': 0.1})

        clr1, kdot, amp, perm = which_osc(osc_type,x1)
        for aa in range(len(y1)):
            ax1.errorbar(x1[aa,amp], y1[aa,cc], yerr1[aa,cc], fmt = marker_shape(x1[aa,1]), markersize=marker_size(x1[aa,1]), color=clr1, mfc=clr1)
        ax1.set_title('$Post-fracture$')

    #     clr2, kdot, amp, perm = which_osc(osc_type,x2)
    #     for bb in range(len(y2)): 
    #         ax1.errorbar(x2[bb,amp], y2[bb], yerr2[bb], fmt = marker_shape(x2[bb,1]), markersize=marker_size(x2[bb,1]), color=clr2, mfc=clr2)

        clr3, kdot, amp, perm = which_osc(osc_type,x3)
        for aa in range(len(y3)):
            ax2.errorbar(x3[aa,amp], y3[aa,cc], yerr3[aa,cc], fmt = marker_shape(x3[aa,1]), markersize=marker_size(x3[aa,1]), color=clr3, mfc=clr3) 
    #     clr4, kdot, amp, perm = which_osc(osc_type,x4)
    #     for bb in range(len(y4)): 
    #         ax2.errorbar(x4[bb,amp], y4[bb], yerr4[bb], fmt = marker_shape(x4[bb,1]), markersize=marker_size(x4[bb,1]), color=clr4, mfc=clr4)
        ax2.set_title('$Post-shear\ 1$')

        clr5, kdot, amp, perm = which_osc(osc_type,x5)
        for aa in range(len(y5)):
            ax3.errorbar(x5[aa,amp], y5[aa,cc], yerr5[aa,cc], fmt = marker_shape(x5[aa,1]), markersize=marker_size(x5[aa,1]), color=clr5, mfc=clr5)
        ax3.set_title('$Post-shear\ 2$')

        #ADD LEGEND TO RIGHT OF LAST PLOT

    #     if osc_type=='NS':
    #         legend1 = plt.legend(handles = legend_elements1N, title='$\sigma_{NS}\ Exp.\ \#$', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    #         legend2 = plt.legend(handles = legend_elements2N, title='$Osc.\ Freq.$', bbox_to_anchor=(1.05, 0.7), loc='upper left', borderaxespad=0.)
    #     else: 
    #         legend1 = plt.legend(handles = legend_elements1P, title='$P_P\ Exp.\ \#$', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    #         legend2 = plt.legend(handles = legend_elements2P, title='$Osc.\ Freq.$', bbox_to_anchor=(1.05, 0.7), loc='upper left', borderaxespad=0.)

    #     ax3.add_artist(legend1); ax3.add_artist(legend2)

        #PLOT 0 LINES
    #     for ax in (ax1,ax2,ax3):
    #         ax.plot([0, 0], [-0.2, 0.1], color = 'gray',ls = '--')
    #         ax.plot([-15, 20], [0, 0], color = 'gray',ls = '--')

        # AXES: MINOR TICKS, COLORS 
        for ax in (ax1,ax2,ax3):
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel('$\Delta k/k_0$ (%)',color=[0,0,0])
        ax2.tick_params(labelleft=False)
        ax3.tick_params(labelleft=False)

        # AXES LABELS & COLORS   
        ax1.ticklabel_format(axis='y',style='sci',scilimits=(-4,2))
        ax1.set_ylabel('$\overline{\Delta c/c_0}$ (%)',color=[0,0,0])
        plt.show()

    if save == 'y': 
#         fig.savefig(filepath+'avgDelc_All'+osc_type+'.png')
        fig.savefig(filepath+'avgDelc_All'+osc_type+'.svg') 