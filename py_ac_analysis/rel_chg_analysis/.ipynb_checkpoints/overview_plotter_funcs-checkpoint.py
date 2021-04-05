import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
# from plotter_functions import *
from legend_specs import *
#----------------------------------------------------------------


def relChg_amp_plots(data, x, y1, y2, color1, color2, xlim, ylim, xlabel, ylabel):
    
    [legend_elements1N, legend_elements1P, legend_elements2N, legend_elements2P, legend_elements2Q, legend_elements2R] = legend_stuff(color1,color2)
    mpl.rcParams.update({'font.size': 12})
        
    ns_runs = ['$\sigma_{eff}$ = 10 MPa', '$\sigma_{eff}$ = 12.5 MPa', '$\sigma_{eff}$ = 15 MPa', '$\sigma_{eff}$ = 17.5 MPa', '$\sigma_{eff}$ = 20 MPa']

    plt.style.use('jgr2020_style3.mplstyle')
    fig, axs = plt.subplots(1, 5, gridspec_kw={'wspace': 0.05})
    
    for aa in range(5):
        axs[aa].errorbar(data[aa*2,x], data[aa*2,y1], yerr=None, fmt = 'o', markersize=6, mfc=color1)
        axs[aa].errorbar(data[(aa*2)+1,x], data[(aa*2)+1,y1], yerr=None, fmt = 's', markersize=6, mfc=color1)
#         fin_idx = np.isfinite(data[aa*2,y1])
#         axs[aa].plot(np.linspace(-0.1, 1.1, 100), np.polyval(np.polyfit(data[aa*2,x][fin_idx], data[aa*2,y1][fin_idx], 1), np.linspace(-0.1, 1.1, 100)), c=color1,ls='--')
#         print(np.polyfit(data[aa,x][fin_idx1], data[aa,y1][fin_idx1], 1)[0])
#         axs[aa].plot(np.linspace(-0.1, 1.1, 100), np.polyval(np.polyfit(data[aa+5,x][fin_idx2], data[aa+5,y1][fin_idx2], 1), np.linspace(-0.1, 1.1, 100)), c='orange',ls='--')
#         print(np.polyfit(data[aa+5,x][fin_idx2], data[aa+5,y1][fin_idx2], 1)[0])
        if y2 != None:
#             fin_idx = np.isfinite(data[aa*2,y2])
#             axs[aa].plot(np.linspace(-0.1, 1.1, 100), np.polyval(np.polyfit(data[aa*2,x][fin_idx], data[aa*2,y2][fin_idx], 1), np.linspace(-0.1, 1.1, 100)), c=color2,ls='--')
            axs[aa].errorbar(data[aa*2,x], data[aa*2,y2], yerr=None, fmt = 'o', markersize=6, mfc=color2)
            axs[aa].errorbar(data[(aa*2)+1,x], data[(aa*2)+1,y2], yerr=None, fmt = 's', markersize=6, mfc=color2)
        axs[aa].set_title(ns_runs[aa])
 
        
    for ax in axs:
        ax.plot([-100, 200], [0, 0], color = 'gray',ls = '--') #Horizontal Line
        ax.plot([0, 0], [-100, 200], color = 'gray',ls = '--') #Horizontal Line
        
    legend1 = plt.legend(handles = legend_elements1P, title='Pair', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    legend2 = plt.legend(handles = legend_elements2P, title='$Osc.\ Set$', bbox_to_anchor=(1.05, 0.65), loc='upper left', borderaxespad=0.)
       
    axs[-1].add_artist(legend1); axs[-1].add_artist(legend2)
    
    # AXES: MINOR TICKS, COLORS 
    for ax in axs:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel,color='black')
#         ax.set_xticks(np.arange(0.2, 1.2, 0.2))
    for ax in axs[1::]:
        ax.tick_params(labelleft=False)
 
    # AXES LABELS & COLORS   
    axs[0].ticklabel_format(axis='y',style='sci',scilimits=(-4,3))
    axs[0].set_ylabel(ylabel, color='black')
    plt.show()

#     if save == 'y': 
#         fig.savefig(filepath+'delc_amp_'+osc_type+'.svg')
#         fig.savefig(filepath+'delc_amp_'+osc_type+'.pdf')




def slope_plots(data, x, y1, y2, color1, color2, xlim, ylim, xlabel, ylabel):
    [legend_elements1N, legend_elements1P, legend_elements2N, legend_elements2P, legend_elements2Q, legend_elements2R] = legend_stuff(color1,color2)
    mpl.rcParams.update({'font.size': 12})
    plt.style.use('jgr2020_style3.mplstyle')
    fig, axs = plt.subplots(1, 1, figsize=(3,3))
    
    ns_levels = np.array([10, 12.5, 15, 17.5, 20])
    q = np.zeros((5,2)); q2 = np.zeros((5,2))
    err = np.zeros(5); err2 = np.zeros(5)

    for aa in range(5):
        fin_idx = np.isfinite(np.hstack(data[(aa*2):(aa*2)+2,y1])) & np.isfinite(np.hstack(data[(aa*2):(aa*2)+2,x]))
        q[aa,:] = np.polyfit(np.hstack(data[(aa*2):(aa*2)+2,x])[fin_idx], np.hstack(data[(aa*2):(aa*2)+2,y1])[fin_idx], 1)
        err[aa] = np.std(np.hstack(data[(aa*2):(aa*2)+2,y1])[fin_idx])
        if y2 != None:
            fin_idx2 = np.isfinite(np.hstack(data[(aa*2):(aa*2)+2,y2])) & np.isfinite(np.hstack(data[(aa*2):(aa*2)+2,x]))
            q2[aa,:] = np.polyfit(np.hstack(data[(aa*2):(aa*2)+2,x])[fin_idx2], np.hstack(data[(aa*2):(aa*2)+2,y2])[fin_idx2], 1)
            err2[aa] = np.std(np.hstack(data[(aa*2):(aa*2)+2,y2])[fin_idx])

    axs.errorbar(ns_levels, q[:,0], yerr=err, fmt = 'd', markersize=6, mfc=color1)
    axs.plot(ns_levels, q[:,0], c=color1, ls='--')
    axs.errorbar(ns_levels, q2[:,0], yerr=err2, fmt = '*', markersize=8, mfc=color2)
    axs.plot(ns_levels, q2[:,0], c=color2, ls='--')

    # AXES: MINOR TICKS, COLORS 
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xlabel(xlabel,color='black')
    if x == 0:
        axs.set_xticks(np.arange(10, 22.5, 2.5))

    legend1 = plt.legend(handles = legend_elements1P, title='Pair', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    # # AXES LABELS & COLORS   
    axs.ticklabel_format(axis='y',style='sci',scilimits=(-4,3))
    axs.set_ylabel(ylabel, color='black')
    plt.show()
    
    return q, err, q2, err2