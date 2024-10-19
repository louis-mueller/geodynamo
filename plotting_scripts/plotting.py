#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:09:57 2024

@author: theresabuettner
"""

from input_data import smoothing, show_rawdata, spacing, params_units, save_fig, output_name, file_type, plot_dir
from coloring import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import signal


#%% PLOTTING INPUT DATA

def plot_input(all_data, homedir, comparison, params, params_units):
    
    foldernames = list(all_data.keys())
    
    if comparison == True:
        foldernames.sort(key=str.lower)
    
    plots = []
    
    for key, value in params['i'].items():
        if value == True:
            plots.append(key)

    plt.rcParams['figure.dpi']=500
    plt.figure() 
    
    nrPlots = len(plots)
    if nrPlots == 1:
        Plx = 1; Ply = 1
        figsize=(5,4)
    elif nrPlots == 2:
        Plx = 1; Ply = 2
        figsize=(10,4)
    elif nrPlots == 3:
        Plx = 1; Ply = 3
        figsize=(15,4)
    elif nrPlots == 4:
        Plx = 2; Ply = 2
        figsize=(10,8)
    elif nrPlots == 5:
        Plx = 2; Ply = 3
        figsize=(15,8)
    elif nrPlots == 6:
        Plx = 2; Ply = 3
        figsize=(15,8)
    elif nrPlots == 7:
        Plx = 3; Ply = 3
        figsize=(15,12)
    elif nrPlots == 8:
        Plx = 3; Ply = 3
        figsize=(15,12)
    elif nrPlots == 9:
        Plx = 3; Ply = 3
        figsize=(15,12)
    elif nrPlots == 10:
        Plx = 4; Ply = 3
        figsize=(15,16)
    elif nrPlots == 11:
        Plx = 4; Ply = 3
        figsize=(15,16)
    elif nrPlots == 12:
        Plx = 4; Ply = 3
        figsize=(15,16)
    elif nrPlots == 13:
        Plx = 5; Ply = 3
        figsize=(15,20)
    else:
        print('Number of plots exceeds or undercuts possible number.')
    
    fig,ax = plt.subplots(Plx,Ply,figsize=figsize)
        
    ind = 0
    colors = colorcoding(foldernames)
    
    for i in range(Plx):
        for k in range(Ply):
                        
            if Plx==1:
                it = k
                if Ply==1:
                    ax = np.array([ax])
            else:
                it = (i,k) # index tuple
            
            count = 0
            
            for entry in foldernames:
                
                radius = [ x / 1000 for x in all_data[entry][list(params_units['i'].keys())[3]] ]

                c = count // 2

                if count % 2 == 0:
                    a = 0.5
                    im=ax[it].plot(radius, all_data[entry][plots[ind]], label=entry, color=colors[c], alpha=a)
                else:
                    a = 1
                    im=ax[it].plot(radius, all_data[entry][plots[ind]], label=entry, color=colors[c], alpha=a)
                
                count = count + 1
            
            ax[it].legend()
                
            ax[it].set_title(plots[ind], pad=10)
            ax[it].set_ylabel(params_units['i'][plots[ind]], labelpad=10)
            
            ax[it].set_xlabel(params_units['i']['Radius'], labelpad=10)
            ax[it].set_xlim((0,None))
            ax[it].legend(fontsize=5)
            
            if ind < nrPlots-1:
                ind = ind + 1
            
    
    if nrPlots == 5:
        adjust = [1, 2]
        
    if nrPlots == 7:
        adjust = [2, 1]
        
    if nrPlots == 11:
        adjust = [3, 2]
    
    if nrPlots == 13:
        adjust = [4, 1]
        
    if 'adjust' in locals():
        
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        
        ax[adjust[0]][adjust[1]].set_visible(False)
        
        chartBox1 = ax[adjust[0]][0].get_position() 
        chartBox2 = ax[adjust[0]][1].get_position()
        chartBox3 = ax[adjust[0]][2].get_position() 
        
        if nrPlots == 7 or nrPlots == 13:
            
            ax[adjust[0]][adjust[1]+1].set_visible(False)
            
            ax[adjust[0]][0].set_position([chartBox1.x0 + (chartBox2.x0-chartBox1.x0), 
                                   chartBox1.y0, 
                                   chartBox1.width, 
                                   chartBox1.height])
            
            ax[adjust[0]][1].set_position([chartBox2.x0 + (chartBox2.x0-chartBox1.x0), 
                                   chartBox2.y0, 
                                   chartBox2.width, 
                                   chartBox2.height])
            
            ax[adjust[0]][2].set_position([chartBox3.x0 + (chartBox2.x0-chartBox1.x0), 
                                   chartBox2.y0, 
                                   chartBox2.width, 
                                   chartBox2.height])
        
        else:
            
            ax[adjust[0]][0].set_position([chartBox1.x0 + 0.5 * (chartBox2.x0-chartBox1.x0), 
                                   chartBox1.y0, 
                                   chartBox1.width, 
                                   chartBox1.height])
            
            ax[adjust[0]][1].set_position([chartBox2.x0 + 0.5 * (chartBox2.x0-chartBox1.x0), 
                                   chartBox2.y0, 
                                   chartBox2.width, 
                                   chartBox2.height])
            
            ax[adjust[0]][2].set_position([chartBox3.x0 + 0.5 * (chartBox2.x0-chartBox1.x0), 
                                   chartBox2.y0, 
                                   chartBox2.width, 
                                   chartBox2.height])
    
    else:
        plt.tight_layout(pad=2.0)
        
    if save_fig == True:
        plt.savefig(plot_dir+output_name+'i'+file_type)
    else:
        plt.show()

#%% PLOTTING HEAT DATA

def plot_heat(all_data, homedir, comparison, params, params_units):
    
    foldernames = list(all_data.keys())
    
    plots = []
    
    for key, value in params['h'].items():
        if value == True:
            plots.append(key)

    plt.rcParams['figure.dpi']=500
    plt.figure() 
    
    nrPlots = len(plots)
    if nrPlots == 1:
        Plx = 1; Ply = 1
        figsize=(5,4)
    elif nrPlots == 2:
        Plx = 1; Ply = 2
        figsize=(10,4)
    elif nrPlots == 3:
        Plx = 1; Ply = 3
        figsize=(15,4)
    elif nrPlots == 4:
        Plx = 2; Ply = 2
        figsize=(10,8)
    elif nrPlots == 5:
        Plx = 2; Ply = 3
        figsize=(15,8)
    else: # 6 subplots
        Plx = 2; Ply = 3
        figsize=(15,8)
    
    fig,ax = plt.subplots(Plx,Ply,figsize=figsize)
        
    ind = 0
    colors = colorcoding(foldernames)
         
    for i in range(Plx):
        for k in range(Ply):
                        
            if Plx==1:
                it = k
                if Ply==1:
                    ax = np.array([ax])
            else:
                it = (i,k) # index tuple
            
            count = 0
            
            if plots[ind] == list(params_units['h'].keys())[4] and smoothing == True:
                
                legendlines = []
                count = 0
                
                for entry in foldernames:
                    
                    windowsize = len(all_data[entry][plots[ind]][::spacing]) // spacing
                    smooth = signal.savgol_filter(all_data[entry][plots[ind]][::spacing], windowsize, 3)
                    
                    if show_rawdata == True:
                    
                        ax[it].plot(all_data[entry][list(params_units['h'].keys())[0]], all_data[entry][plots[ind]], color=colors[count], alpha=0.5)
                    
                    p, = ax[it].plot(all_data[entry][list(params_units['h'].keys())[0]][::spacing], smooth, color=colors[count], label = entry)
                    
                    legendlines.append(p)
                    
                    count = count + 1
                    
                ax[it].legend(legendlines, foldernames)
            
            
            else:
                
                for entry in foldernames:
                    im=ax[it].plot(all_data[entry][list(params_units['h'].keys())[0]], all_data[entry][plots[ind]], label=entry)
                
                ax[it].legend()
                
            ax[it].set_title(plots[ind], pad=10)
            ax[it].set_ylabel(params_units['h'][plots[ind]], labelpad=10)
            
            ax[it].set_xlabel(params_units['h']['Time'], labelpad=10)
            ax[it].set_xlim((0,None))
            ax[it].legend(fontsize=8)
            
            if ind < nrPlots-1:
                ind = ind + 1
            
    
    if nrPlots == 5:
        
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        
        ax[1][2].set_visible(False)
        
        chartBox1 = ax[1][0].get_position() 
        chartBox2 = ax[1][1].get_position() 
        
        ax[1][0].set_position([chartBox1.x0 + 0.5 * (chartBox2.x0-chartBox1.x0), 
                               chartBox1.y0, 
                               chartBox1.width, 
                               chartBox1.height])
        
        ax[1][1].set_position([chartBox2.x0 + 0.5 * (chartBox2.x0-chartBox1.x0), 
                               chartBox2.y0, 
                               chartBox2.width, 
                               chartBox2.height])
        
    else:
        plt.tight_layout(pad=2.0)
        
    if save_fig == True:
        plt.savefig(plot_dir+output_name+'h'+file_type)
    else:
        plt.show()


#%% PLOTTING 2D DATA

def plot_2D(all_data, homedir, comparison, params, params_units):
    
    foldernames = list(all_data.keys())
    
    plots = []
    
    for key, value in params['b'].items():
        if value == True:
            plots.append(key)
    
    #def plot_results(variables, temp, visc, vel, shf):
        
    plt.rcParams['figure.dpi']=500
    plt.figure() 
    
    colormap = colorcoding_b()

    if comparison == False:
        
        nrPlots = len(plots)
        if nrPlots == 1:
            Plx = 1; Ply = 1
            figsize=(5,4)
        elif nrPlots == 2:
            Plx = 1; Ply = 2
            figsize=(10,4)
        elif nrPlots == 3:
            Plx = 1; Ply = 3
            figsize=(15,4)
        elif nrPlots == 4:
            Plx = 2; Ply = 2
            figsize=(10,8)
            
        fig,ax = plt.subplots(Plx,Ply,figsize=figsize)
        
        ind = 0

        for i in range(Plx):
            for k in range(Ply):
                
                if Plx==1:
                    it = k
                    if Ply==1:
                        ax = np.array([ax])
                else:
                    it = (i,k) # index tuple
                
                if plots[ind] == list(params['b'].keys())[3]:
                
                    for entry in foldernames:
                        
                        time = [item[0]/1e6 for item in all_data[entry][plots[ind]]]
                        data_pure = [item[1:] for item in all_data[entry][plots[ind]]]
                        lateral_cells = list(range(1,len(all_data[entry][plots[ind]][0])))
                        
                        im = ax[it].pcolormesh(lateral_cells, time, data_pure, cmap=colormap[plots[ind]])  
                        
                    ax[it].set_xlabel('Lateral cells', labelpad=10)
                    ax[it].set_ylabel('Time [Myr]', labelpad=10)
                
                else:
                    
                    for entry in foldernames:
                        
                        time = [item[0]/1e6 for item in all_data[entry][plots[ind]]]
                        data_pure = [item[1:] for item in all_data[entry][plots[ind]]]
                        data_trans = np.array(data_pure).T.tolist()
                        radial_cells = list(range(1,len(all_data[entry][plots[ind]][0])))
                        
                        if plots[ind] == list(params['b'].keys())[1]:
                            
                            im = ax[it].pcolormesh(time, radial_cells, data_trans, cmap=colormap[plots[ind]], norm=mpl.colors.LogNorm())
                        
                        else:
        
                            #if plots[ind] == list(params['b'].keys())[0]:
                            #    im = ax[it].pcolormesh(time, radial_cells, data_trans, cmap=colormap['Radial Velocity'])
                            
                            #if plots[ind] == list(params['b'].keys())[2]:
                            #    im = ax[it].pcolormesh(time, radial_cells, data_trans, cmap=colormap[plots[ind]])
                            
                            im = ax[it].pcolormesh(time, radial_cells, data_trans, cmap=colormap[plots[ind]])
                        
                    ax[it].set_xlabel('Time [Myr]', labelpad=10)
                    ax[it].set_ylabel('Radial cells', labelpad=10)
                    
                ax[it].set_title(plots[ind], pad=10)
                fig.colorbar(im, ax=ax[it], label=params_units['b'][plots[ind]])
                    
                ind = ind + 1
                    
        plt.tight_layout(pad=2.0)
        
        if save_fig == True:
            plt.savefig(plot_dir+output_name+'b.png')
        else:
            plt.show()
    

    else: 
        
        nrPlots = len(foldernames)
        if nrPlots == 1:
            Plx = 1; Ply = 1
            figsize=(6,4)
        elif nrPlots == 2:
            Plx = 1; Ply = 2
            figsize=(12,4)
        elif nrPlots == 3:
            Plx = 1; Ply = 3
            figsize=(18,4)
        elif nrPlots == 4:
            Plx = 2; Ply = 2
            figsize=(12,8)
            
            
        for plot in plots:
            
            fig,ax = plt.subplots(Plx,Ply,figsize=figsize)
            
            ind = 0
            
            for i in range(Plx):
                for k in range(Ply):
                
                    if Plx==1:
                        it = k
                        if Ply==1:
                            ax = np.array([ax])
                    else:
                        it = (i,k) # index tuple
                    
                    if plot == list(params['b'].keys())[3]:
                        
                        time = [item[0]/1e6 for item in all_data[foldernames[ind]][plot]]
                        data_pure = [item[1:] for item in all_data[foldernames[ind]][plot]]
                        lateral_cells = list(range(1,len(all_data[foldernames[ind]][plot][0])))
                        
                        im = ax[it].pcolormesh(lateral_cells, time, data_pure, cmap=colormap[plot])  
                        
                        ax[it].set_xlabel('Lateral cells', labelpad=10)
                        ax[it].set_ylabel('Time [Myr]', labelpad=10)
                        
                    else:
                        
                        time = [item[0]/1e6 for item in all_data[foldernames[ind]][plot]]
                        data_pure = [item[1:] for item in all_data[foldernames[ind]][plot]]
                        data_trans = np.array(data_pure).T.tolist()
                        radial_cells = list(range(1,len(all_data[foldernames[ind]][plot][0])))
                        
                        if plot == list(params['b'].keys())[1]:
                            
                            colormap = 'YlGnBu'
                            im = ax[it].pcolormesh(time, radial_cells, data_trans, cmap=colormap[plot], norm=mpl.colors.LogNorm())
                        
                        else:
        
                            #if plot == list(params['b'].keys())[0]:
                             #   colormap = 'summer'
                            
                            #if plot == list(params['b'].keys())[2]:
                             #   colormap = 'magma'
                            
                            im = ax[it].pcolormesh(time, radial_cells, data_trans, cmap=colormap[plot])
                        
                        ax[it].set_xlabel('Time [Myr]', labelpad=10)
                        ax[it].set_ylabel('Radial cells', labelpad=10)
                        
                    ax[it].set_title(foldernames[ind], pad=10)
                    
                    ind = ind + 1

            
            #plt.tight_layout(pad=2.0)
            fig.subplots_adjust(right=0.8)
            # colorbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(im, ax=ax[:], cmap = colormap[plot], label=params_units['b'][plot])

            if save_fig == True:
                plt.savefig(plot_dir+output_name+'b.png')
            else:
                plt.show()

#%% PLOTTING CORE DATA

def plot_core(all_data, homedir, comparison, params, params_units):
    
    foldernames = list(all_data.keys())

    if comparison == True:
        foldernames.sort(key=str.lower)
    
    plots = []
    
    for key, value in params['c'].items():
        if value == True:
            plots.append(key)

    plt.rcParams['figure.dpi']=500
    plt.figure() 
    
    nrPlots = len(plots)
    if nrPlots == 1:
        Plx = 1; Ply = 1
        figsize=(5,4)
    elif nrPlots == 2:
        Plx = 1; Ply = 2
        figsize=(10,4)
    elif nrPlots == 3:
        Plx = 1; Ply = 3
        figsize=(15,4)
    elif nrPlots == 4:
        Plx = 2; Ply = 2
        figsize=(10,8)
    elif nrPlots == 5:
        Plx = 2; Ply = 3
        figsize=(15,8)
    elif nrPlots == 6:
        Plx = 2; Ply = 3
        figsize=(15,8)
    elif nrPlots == 7:
        Plx = 3; Ply = 3
        figsize=(15,12)
    elif nrPlots == 8:
        Plx = 3; Ply = 3
        figsize=(15,12)
    elif nrPlots == 9:
        Plx = 3; Ply = 3
        figsize=(15,12)
    elif nrPlots == 10:
        Plx = 4; Ply = 3
        figsize=(15,16)
    elif nrPlots == 11:
        Plx = 4; Ply = 3
        figsize=(15,16)
    elif nrPlots == 12:
        Plx = 4; Ply = 3
        figsize=(15,16)
    elif nrPlots == 13:
        Plx = 5; Ply = 3
        figsize=(15,20)
    else:
        print('Number of plots exceeds or undercuts possible number.')
    
    fig,ax = plt.subplots(Plx,Ply,figsize=figsize)
        
    ind = 0
    colors = colorcoding(foldernames)
    
    for i in range(Plx):
        for k in range(Ply):
                        
            if Plx==1:
                it = k
                if Ply==1:
                    ax = np.array([ax])
            else:
                it = (i,k) # index tuple
            
            count = 0
            
            for entry in foldernames:
                
                time = np.array(all_data[entry][list(params_units['c'].keys())[0]])
                c = count // 2

                yyy = np.array(all_data[entry][plots[ind]])

                if plots[ind] == list(params_units['c'].keys())[10]:
                    yyy = np.array(all_data[entry][plots[ind]]) * 1e6
                elif plots[ind] == list(params_units['c'].keys())[11]:
                    yyy = np.array(all_data[entry][plots[ind]]) * 1e6
                    if entry == foldernames[0]:
                        ax[it].fill_between(time, 22, 62, color='lightblue', label='Today\'s Suface B-Field', alpha=0.5)

                if count % 2 == 0:
                    a = 0.5
                    im=ax[it].plot(time, yyy, label=entry, color=colors[c], alpha=a)
                else:
                    a = 1
                    im=ax[it].plot(time, yyy, label=entry, color=colors[c], alpha=a)
                
                count = count + 1
            
            ax[it].legend()
                
            ax[it].set_title(plots[ind])
            ax[it].set_xlabel(params_units['c']['Time'])
            ax[it].set_ylabel(params_units['c'][plots[ind]])
            ax[it].set_xlim((0,None))
            ax[it].legend(fontsize=5)
            
            if ind < nrPlots-1:
                ind = ind + 1
            
    
    if nrPlots == 5:
        adjust = [1, 2]
        
    if nrPlots == 7:
        adjust = [2, 1]

    if nrPlots == 9:
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        
    if nrPlots == 11:
        adjust = [3, 2]

    if nrPlots == 12:
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
    
    if nrPlots == 13:
        adjust = [4, 1]
        
    if 'adjust' in locals():
        
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        
        ax[adjust[0]][adjust[1]].set_visible(False)
        
        chartBox1 = ax[adjust[0]][0].get_position() 
        chartBox2 = ax[adjust[0]][1].get_position()
        chartBox3 = ax[adjust[0]][2].get_position() 
        
        if nrPlots == 7 or nrPlots == 13:
            
            ax[adjust[0]][adjust[1]+1].set_visible(False)
            
            ax[adjust[0]][0].set_position([chartBox1.x0 + (chartBox2.x0-chartBox1.x0), 
                                   chartBox1.y0, 
                                   chartBox1.width, 
                                   chartBox1.height])
            
            ax[adjust[0]][1].set_position([chartBox2.x0 + (chartBox2.x0-chartBox1.x0), 
                                   chartBox2.y0, 
                                   chartBox2.width, 
                                   chartBox2.height])
            
            ax[adjust[0]][2].set_position([chartBox3.x0 + (chartBox2.x0-chartBox1.x0), 
                                   chartBox2.y0, 
                                   chartBox2.width, 
                                   chartBox2.height])
        else:
        
            ax[adjust[0]][0].set_position([chartBox1.x0 + 0.5 * (chartBox2.x0-chartBox1.x0), 
                                    chartBox1.y0, 
                                    chartBox1.width, 
                                    chartBox1.height])
            
            ax[adjust[0]][1].set_position([chartBox2.x0 + 0.5 * (chartBox2.x0-chartBox1.x0), 
                                    chartBox2.y0, 
                                    chartBox2.width, 
                                    chartBox2.height])
            
            ax[adjust[0]][2].set_position([chartBox3.x0 + 0.5 * (chartBox2.x0-chartBox1.x0), 
                                    chartBox2.y0, 
                                    chartBox2.width, 
                                    chartBox2.height])
                    
    if save_fig == True:
        plt.savefig(plot_dir+output_name+'c'+file_type)
    else:
        plt.show()