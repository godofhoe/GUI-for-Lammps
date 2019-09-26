# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 23:59:22 2016

@author: Lina492375qW1188
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
from scipy.stats import binned_statistic_2d

class Intensity2D:
    
    def __init__(self, cmin, cmax, mag):
        
        self.cmin = cmin # minimum of colorbar
        self.cmax = cmax # maximum of colorbar
        self.mag = mag   # magnification of z-value
        
    def MakePlot(self, bs, ticks_list, Lognormal = None):
        """if need lognormal color distribution, please provide 4th argument."""
        try:    
            bz = bs.statistic
            
            # norm: normalize colorbar.
            if Lognormal is None:
                norm = (cm.colors.Normalize(vmin=np.nanmin(bz)*self.cmin,
                                            vmax=np.nanmax(bz)*self.cmax))
            else:
                norm = (LogNorm(vmin=np.nanmax(bz)*self.cmin,
                                vmax=np.nanmax(bz)*self.cmax))
                
            fig = plt.figure()
            ax = fig.add_subplot(111)
            im = (ax.imshow(bz, extent=(np.amin(bs.y_edge), np.amax(bs.y_edge),
                                    np.amin(bs.x_edge), np.amax(bs.x_edge)),
                                    interpolation='sinc', origin='lower',
                                    aspect='auto', cmap=cm.afmhot, norm=norm))
            
            fig.colorbar(im, ticks = ticks_list)

        except:
            print('bs should be <scipy.stats._binned_statistic.BinnedStatistic2dResult>!')
          
        plt.xlabel('y')
        plt.ylabel('x')

        plt.show()         
        
    def EnergyMap(self, coord, Nx, Ny):
               
        x, y, z = zip(*coord)
        z = [element**self.mag for element in z]
        
        bs = binned_statistic_2d(x, y, z, statistic='mean', bins=[Nx, Ny])
        
        return bs
        
    def ScipyBinned(self, coord, Nx, Ny): # coord = [(x,y,z)]
        
        bs = self.EnergyMap(coord, Nx, Ny)
        
        bz = bs.statistic
        d = (np.nanmax(bz)-np.nanmin(bz))/4            
        ticks_list = ([np.nanmin(bz), np.nanmin(bz)+d, np.nanmin(bz)+2*d, 
                       np.nanmin(bz)+3*d, np.nanmax(bz)])
    
        self.MakePlot(bs, ticks_list)
