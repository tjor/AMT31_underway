#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 09:41:42 2022

@author: tjor

This script contains extra plot functions for step3 data from DY151


"""

import os
import sys
import csv
import netCDF4 as nc
import pandas as pd

import numpy as np
import glob   
import pandas as pd
import ephem

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import math
import matplotlib.dates as mdates
import datetime as dt
from scipy.interpolate import interp1d
import scipy.io as io
import datetime


import pandas as pd
import glob

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def GL_map(fname):
    
        plt.figure(figsize=(14,8))#
        plt.rc('font', size=24)
        
        min_lat = np.floor(np.nanmin(lat)) - 1
        max_lat = np.ceil(np.nanmax(lat))  + 1
        min_lon  = np.min(np.nanmin(lon)) - 1
        max_lon  = np.max(np.nanmax(lon)) + 1
        extent = [min_lon, max_lon, min_lat, max_lat] 
        #extent = [-4.40, -4.10, 50.00, 50.40] 
        
        request = cimgt.Stamen(style='terrain')
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent, ccrs.PlateCarree())
        ax.add_image(request,8)
    
        gl = ax.gridlines(draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter =  LONGITUDE_FORMATTER
        gl.yformatter =  LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 18,  'rotation': 0}
        gl.ylabel_style = {'size': 18,  'rotation': 0}
        
        
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.tick_params(labelsize = 10)
    
        time_plot = [mdates.date2num(time[i]) for i in range(len(time))]
        loc = mdates.AutoDateLocator()
        sc = plt.scatter(lon,lat,s=10, c = time_plot, cmap ='plasma', vmin = time_plot[0], vmax =time_plot[-1])
       
               
        filename  = dir_figs +  fname + '.png'
        plt.savefig(filename, dpi=400)
         
       
        return
    
   
def GL_time_series(fname):
        
        plt.figure(figsize=(24,10))
        plt.rc('font', size=24)
        plt.ylabel('Chl [mg m$^{-3}$]')
        colors = cm.plasma(np.linspace(0, 1, len(time))) # heat map for color  bar
        for i in range(len(time)):
                plt.plot_date(time[i], chl[i], color=colors[i], ms=2) # assumes wl vector
        plt.gca().set_yscale('log')
        plt.ylim(0.03,16)
        plt.xlabel('Date')
        plt.title('Chlorophyll-a concentrations')
        plt.xticks(rotation=45)
        
        filename  =  dir_figs +  fname + '.png'
        plt.savefig(filename, dpi=400)
        
        return

if __name__ == '__main__':
        
    # directories 
    dir_main = '/data/datasets/cruise_data/active/DY151/'
    dir_data = dir_main + 'Processed/Underway/Step3/'
    path_data= dir_data + 'dy151_optics.mat'
    dir_figs = dir_main + 'Figures/Underway/Summary_figures_python'
   
    # unpack step3 data structure
    data = io.loadmat(path_data, squeeze_me = True, struct_as_record = False) 
    
    # unpack timestamp in datetime format, chl, lat and lon
    time_num = data['dy151'].uway.time
    first_time = dt.datetime.fromordinal(-365 + int(np.floor(time_num[60]))) # first date of data
    time = np.nan*np.ones(len(time_num), dtype=object)
    for i in range(len(time_num)):
      secs = (time_num[i] - np.floor(time_num[0]))*86400
      time[i] = first_time + datetime.timedelta(seconds=secs)
    chl = data['dy151'].acs.chl
    lat = data['dy151'].uway.lat
    lon = data['dy151'].uway.long
    
    # remove nans
    time = time[np.isnan(chl)==0]
    lat = lat[np.isnan(chl)==0]
    lon = lon[np.isnan(chl)==0]
    chl =  chl[np.isnan(chl)==0]
    
    GL_map('/map_DY151_cruise')
    GL_time_series('/TS_DY151_cruise')

    # zoom plot - currently hard coded
    time_zoom = dt.datetime(2022,6,18,12,0,0)
    tol = int(12*3600)       # time window tolerance in seconds
    delta = np.abs(np.array([(time_zoom - time[i]).total_seconds() for i in range(len(time))]))   

    lat = lat[delta < tol] #  filter in-situ rrs subject to time window mask
    lon = lon[delta < tol]
    chl = chl[delta < tol]
    time = time[delta < tol]

    GL_map('/map_18_june')
    GL_time_series('/TS_18_june')