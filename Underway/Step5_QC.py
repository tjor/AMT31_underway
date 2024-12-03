#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 15:12:03 2024

@author: tjor
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a python script version of the AMT QC script, used on James cook
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy import signal as sg
from datetime import datetime as dt
def prcrng(x):
    return (np.nanpercentile(x,84) - np.nanpercentile(x,16))/2.

import matplotlib.dates as mdates
import matplotlib.cm as cm

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def plot_timeseries(acs):
        
        plt.figure(figsize=(18,20))
        plt.rcParams.update({'font.size': 16})
        plt.subplot(3,1,1)
        
        
        colors = cm.inferno(np.linspace(0, 1, len((acs.time))))# heat map for color  bar
        #plt.plot_date(acs.time, acs.acs_chl, color='gray', ms=2) # as
        for i in range(len(acs.time)):
                plt.plot_date(acs.time[i], acs.acs_chl[i], color=colors[i], ms=2) # assumes wl vector
        
        plt.gca().set_yscale('log')
        plt.ylim(0.05,1)
        plt.ylabel('Chl  (AC-S) [mg m$^{-3}$]')
        #plt.xlabel('Date')
        ax = plt.gca()
        ax.set_xticks([])  
        #plt.title('Chlorophyll-a concentration (AC-S)')
                
           
        plt.subplot(3,1,2)
        plt.plot_date(acs.time, acs.uway_sst, ms=2, color='red' ,alpha=0.5,label='Intake') # assumes, 
        plt.plot_date(acs.time, acs.ctd_sst, ms=2, color='blue',alpha=0.5,label='Flow-through')
        plt.ylim(12,25)
        plt.ylabel('SST [$^{\circ}$C]')
        plt.legend()
        #lt.xlabel('Date')
        ax = plt.gca()
        ax.set_xticks([])  
        #plt.title('Sea surface temperature')
        
        plt.subplot(3,1,3)
        plt.plot_date(acs.time, acs.uway_sal, ms=2, color='red', alpha=0.5, label='Intake') # assumes, 
        plt.plot_date(acs.time, acs.ctd_sal, ms=2, color='blue', alpha=0.5,label='Flow-through') # assumes, 
        plt.ylim(34,37)
        plt.ylabel('Salinity [PSU]')
        plt.legend()
        plt.xlabel('Date')
        #plt.title('Sea surface temperature')
        plt.tight_layout()
 
        
        filename  =  fig_dir + '/'  + '_AMT31_timeseries.png'
        plt.savefig(filename,dpi=300)    
        
        return
    

def plot_latitudinal(acs):    
        
        plt.figure(figsize=(18,20))
        plt.rcParams.update({'font.size': 16})
    
        plt.subplot(3,1,1)
        colors = cm.inferno(np.linspace(0, 1, len((acs.time))))# heat map for color  bar
        #plt.plot_date(acs.time, acs.acs_chl, color='gray', ms=2) # as
        for i in range(len(acs.time)):
                plt.scatter(acs.uway_lat[i], acs.acs_chl[i], color=colors[i], s=5) # assumes wl vector
        plt.gca().set_yscale('log')
        plt.ylim(0.1,1)
        plt.ylabel('Chl  (AC-S) [mg m$^{-3}$]')
        #plt.xlabel('Date')
        ax = plt.gca()
        ax.set_xticks([])  
        plt.gca().invert_xaxis()
        #plt.title('Chlorophyll-a concentration (AC-S)')
        
           
        plt.subplot(3,1,2)
        plt.scatter(acs.uway_lat, acs.uway_sst, s=5, color='red' ,alpha=0.5,label='Intake') # assumes, 
        plt.scatter(acs.uway_lat, acs.ctd_sst, s=5, color='blue',alpha=0.5,label='Flow-through')
        # plt.scatter(acs.uway_lat, acs.uway_sst, s=5, color='red' ,alpha=0.5,label='Intake') # assumes, 
        # plt.scatter(acs.uway_lat, acs.ctd_sst, s=5, color='blue',alpha=0.5,label='Flow-through')
        plt.ylim(12,22)
        plt.ylabel('SST [$^{\circ}$C]')
        plt.legend()
        # plt.xlabel('Date')
        ax = plt.gca()
        ax.set_xticks([])  
        plt.gca().invert_xaxis()
        #plt.title('Sea surface temperature')
        
        plt.subplot(3,1,3)
        
        plt.scatter(acs.uway_lat, acs.uway_sal, s=5, color='red', alpha=0.5, label='Intake') # assumes, 
        plt.scatter(acs.uway_lat, acs.ctd_sal, s=5, color='blue', alpha=0.5,label='Flow-through') # assumes, 
        #  plt.scatter(acs.lat, acs.uway_sal, s=5, color='red', alpha=0.5, label='Intake') # assumes, 
        #  plt.scatter(acs.lat, acs.ctd_sal, s=5, color='blue', alpha=0.5,label='Flow-through') # assu
        plt.ylim(34,38)
        plt.ylabel('Salinity [PSU]')
        plt.legend()
        plt.xlabel('Latitude')
        # plt.title('Sea surface temperature')
        plt.tight_layout()
        plt.gca().invert_xaxis()
              
        filename  =  fig_dir + '/'  + '_AMT31_latitudinal.png'
        plt.savefig(filename,dpi=300)    

        
        return

def plot_coverage(acs):
    
         plt.figure(figsize=(14,8))#
         plt.rc('font', size=24)
         
         lat = acs.uway_lat
         lon = acs.uway_long
         time = acs.time
         
         min_lat = np.floor(np.nanmin(lat)) - 1
         max_lat = np.ceil(np.nanmax(lat))  + 1
         min_lon  = np.min(np.nanmin(lon)) - 1
         max_lon  = np.max(np.nanmax(lon)) + 1
         extent = [min_lon, max_lon, min_lat, max_lat] 
         #extent = [-4.40, -4.10, 50.00, 50.40] 
         
         request = cimgt.GoogleTiles(style='satellite')
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
         sc = plt.scatter(lon,lat,s=10, c = time_plot, cmap ='inferno', vmin = time_plot[0], vmax =time_plot[-1])

         filename  =  fig_dir + '/'  + '_Coverage.png'
         plt.savefig(filename,dpi=300)    

         return
     
def _167_122(acs):

            plt.figure()
            plt.plot(acs.acs_chl,label='122')
            plt.plot(acs.acs2_chl,label='167')
            plt.legend()
            plt.xlim(950,1100)
            plt.ylim(0,0.8)
            plt.ylabel('Chl  (AC-S) [mg m$^{-3}$]')
            plt.xlabel('Minute in day')
            
        
            index = 975
            plt.figure()
            plt.plot(acs.acs_wv, acs.acs_ap.values[index],label='167')
            plt.plot(acs.acs_wv,acs.acs2_ap.values[index],label='122')
            plt.legend()
        
        
            index = 1000
            plt.figure()
            plt.plot(acs.acs_wv,acs.acs_ap.values[index],label='167')
            plt.plot(acs.acs_wv,acs.acs2_ap.values[index],label='122')
            plt.legend()
        
        
            index = 1050
            plt.figure()
            plt.plot(acs.acs_wv,acs.acs_ap.values[index],label='167')
            plt.plot(acs.acs_wv,acs.acs2_ap.values[index],label='122')
            plt.legend()
        
        
            plt.figure()
            plt.plot(acs.acs_wv,np.nanmedian(np.array(acs.acs_ap.values[900:1100]), axis=0),label='167')
            plt.plot(acs.acs_wv,np.nanmedian(np.array(acs.acs2_ap.values[900:1100]), axis=0),label='122')
            plt.xlabel('Wavelength [nm]')
            plt.ylabel('m$^{-1}$')
            plt.legend()
        
            return
        
def flow_sal_filter(acs, MIN_FLOW_RATE = 26, MIN_SAL =35):

   # MIN_FLOW_RATE = 25
   # MIN_SAL = 35

    i2f1 = np.where((acs.ctd_sal > MIN_SAL) & (acs.flow > MIN_FLOW_RATE))[0]  
    i2f1 = np.where((acs.flow > MIN_FLOW_RATE))[0] 
    i2f2 = np.where((np.isnan(acs.uway_sal)==True) | (np.isnan(acs.flow)==True))[0]
    
    i2f = np.union1d(i2f1, i2f2)
   
    fig, ax = plt.subplots(2,1, figsize=(13, 12), sharex=True)
    ax[0].plot(acs.time, acs.flow, '.-', lw=0.5, ms=1, alpha=0.5)
    ax[0].plot(acs.time[i2f], acs.flow[i2f], 'ro', lw=0.5, ms=5, mfc='none', alpha=0.15)
    ax[0].set_ylabel('flow')
    ax[0].grid('on')
    ax[0].set_ylim([-1, 60])


    # fig0, ax0 = plt.subplots(1, figsize=(15, 5))
    ax[1].plot(acs.time, acs.uway_sal, '.-', lw=0.5, ms=1, alpha=0.5)
    ax[1].plot(acs.time[i2f], acs.uway_sal[i2f], 'r.', lw=0.1, ms=3, mfc='none', alpha=0.15)
    ax[1].set_ylabel('salinity')
    ax[1].grid('on')
    
    return i2f


def QC_filtering_2systems(acs):
        
        # acs
        # step x - # filters w.r.t. i2fn (mQ interval + previous manual spike removal) and innan
        ix = xr.DataArray(acs.time[i2f][innan], dims=['time']) 
        acs_filtered_x = acs.sel(time = ix)
        #acs_filtered_x['acs_chl'].values = sg.medfilt(acs.acs_chl[i2f][innan], kernel_size=31)# Filter-size set to 1 as default #
      
        # step y # filters w.r.t. ap (10th element) being > 0
        i2kp = np.where((acs_filtered_x.acs_ap[:,10] > 0)) [0]
        iy = xr.DataArray(acs_filtered_x.time[i2kp], dims=['time']) 
        acs_filtered_xy = acs_filtered_x.sel(time = iy)
        
        # step z - # filters w.r.t. acs chl being > 0
        i2kp = np.where((acs_filtered_xy.acs_chl[:] > 0)) [0]
        iz = xr.DataArray(acs_filtered_xy.time[i2kp], dims=['time']) 
        acs_filtered_xyz = acs_filtered_xy.sel(time = iz)

        
        # acs 2
        # step x # filters w.r.t. i2fn (mQ interval + previous manual spike removal) and innan
        ix = xr.DataArray(acs.time[i2f][innan2], dims=['time']) 
        acs2_filtered_x = acs.sel(time = ix)
        #acs2_filtered_x['acs2_chl'].values = sg.medfilt(acs.acs2_chl[i2f][innan2], kernel_size=31) # apply med filter to acx chl and acx chl debiasesd 
   
    
        # step y - # filters w.r.t. ap (10th element) being > 0
        i2kp = np.where((acs2_filtered_x.acs2_ap[:,10] > 0))[0]
        iy = xr.DataArray(acs2_filtered_x.time[i2kp], dims=['time']) 
        acs2_filtered_xy = acs2_filtered_x.sel(time = iy)
        
        # step z - # filters w.r.t. acs chl being > 0
        i2kp = np.where((acs2_filtered_xy.acs2_chl[:] > 0)) [0]
        iz = xr.DataArray(acs2_filtered_xy.time[i2kp], dims=['time']) 
        acs2_filtered_xyz = acs2_filtered_xy.sel(time = iz)
        
        
        # combinging 2 filtered files
        key_string = acs_filtered_xyz.keys()
        hplc_keys = [s for s in key_string if "hplc" in s]
        other_keys = [s for s in key_string if "hplc" not in s]
    
        
        # acs_filtered  = xr.concat([acs_filtered_xyz, acs2_filtered_xyz],dim='time') # combines t 
        # breakpoint()
        acs_filtered = xr.merge([acs_filtered_xyz, acs2_filtered_xyz],compat='no_conflicts') 
       
        return acs_filtered

if __name__ == '__main__':
    
    
    # load ACS file
    DIN_acs = "/mnt/d/AMT31/Optics_all/Processed/Step3/"
    fn_acs = "amt31_iop_flowthrough.nc"
    acs = xr.open_dataset(DIN_acs + fn_acs)
    
    fig_dir = "/mnt/d/AMT31/Optics_all/Figures/Summary/"
    
    #if "uway_long" in acs.keys():
      #  acs['uway_lon'] = acs['uway_long']
       # acs = acs.drop(labels="uway_long")

    # initial QC masks for filtered (f) and nan-padding intervals
    i2f = flow_sal_filter(acs)
    innan = np.where(~np.isnan(acs.acs_chl[i2f]))[0] # need to remove nans to prevent medfilt to be spiky near edges
    innan2 = np.where(~np.isnan(acs.acs2_chl[i2f]))[0] # need to 

    # QC filtering
    acs_filtered = QC_filtering_2systems(acs)

    #temp overwrite

    acs_filtered['acs_chl'].values[2329:] = acs_filtered['acs2_chl'].values[2329:]
     
    plot_timeseries(acs_filtered)
    plot_latitudinal(acs_filtered)
    plot_coverage(acs_filtered)

    # write to file      
    acs_filtered.to_netcdf(DIN_acs + fn_acs[:-3] + '_QC.nc')
    acs_filtered.close()
   
    plt.figure()
    plt.plot(np.log10(acs_filtered['acs_chl']))
  #  plt.figure(figsize=(12,5))
  #  import matplotlib.dates as mdates
  #  plt.rcParams.update({'font.size': 18})
  #  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
  #  plt.ylabel('Chl-a [mg m$^{-3}]$')
    #plt.plot_date(acs_filtered_x['acs_chl'].time, acs_filtered_x['acs_chl'].values, ms=5, label = 'Failed QC')
  #  plt.plot_date(acs_filtered_xyz['acs_chl'].time, acs_filtered_xyz['acs_chl'].values, ms=5, label = 'Passed QC')
  #  plt.ylim(0,0.7)
   # plt.xlim()
   # plt.xlabel('Time [UTC]')
 
   # plt.figure(figsize=(12,5))
  #  import matplotlib.dates as mdates
   # plt.rcParams.update({'font.size': 18})
   # plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
  ##  plt.ylabel('Lat [degs]')
   # plt.plot_date(acs.time, acs.uway_lat, ms=5, label = 'Passed QC')
  #  plt.xlabel('Time [UTC]')


 