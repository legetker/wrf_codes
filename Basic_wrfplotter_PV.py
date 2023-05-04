#!/usr/bin/env python
# coding: utf-8

# Basic Python program using wrf-python capabilities to generate basic meteorology plots from wrfout files
# Program is designed for further development, customization, and improvememt.  Adding computed parameters is the next step.
# 
# See https://wrf-python.readthedocs.io/en/latest/index.html, from where examples were drawn.
# 
# Gary Lackmann, for MEA716, December 2020

from __future__ import print_function
import os
#from os import listdir
#from os.path import isfile, join
import glob

import subprocess
import wrf
import numpy as np
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.cbook as cbook
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

from wrf import (to_np, interplevel, geo_bounds, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)
from matplotlib.colors import Normalize
from metpy.plots import ctables

# Specify which plots you want to generate:  slp+radar, 850, 500, and 250 all default to yes
doSLP =  False
do850 = True
do500 = False
do250 = False
doPV = False

# Set directory where wrfout files reside, and list the files for processing.  Set up for a directory with only wrfout files.
os.chdir("/gpfs_common/mea716_share/gary/Dec2020_winterstorm/WRF-4.2.1/test/em_real/")
os.getcwd()

datafiles = (glob.glob("./wrfout*"))
numfiles=len(datafiles)
print(numfiles)
#print(datafiles[0])

# Set up for saving graphics files
plotsdir = './plots/'

# create directory for plots, if needed
if os.path.isdir(plotsdir) != 1:
    subprocess.call(["mkdir","-p",plotsdir])

# Set up loop to run same plots for multiple different times (assume here that we have 1 time per wrfout file)
for j in range(0,numfiles-1):
    #for j in range(0,1):
    ncfile = Dataset(datafiles[j])
    Time=wrf.extract_times(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=False, do_xtime=False)
    timestr=(str(Time))
    # Set up one time string for plot titles, another for file names
    titletime=(timestr[0:10]+' '+timestr[11:16])
    filetime=(timestr[0:10]+'_'+timestr[11:13])
    print('WRF valid time: ',filetime)

    # Get all the variables we need
    z = getvar(ncfile, "z")
    dbz3 =getvar(ncfile, "dbz")
    p = getvar(ncfile, "pressure")
    slp = getvar(ncfile, "slp")
    mdbz = getvar(ncfile, "mdbz")
    tk = getvar(ncfile, "tk")
    ua = getvar(ncfile, "ua", units="kt")
    va = getvar(ncfile, "va", units="kt")
    pvo = getvar(ncfile, "pvo")
    wspd = getvar(ncfile, "wspd_wdir", units="kts")[0,:]

    # Do vertical interpolation to needed pressure surfaces - could do a loop over vertical levels at some point
    ht_850 = interplevel(z, p, 850.)
    ht_500 = interplevel(z, p, 500.)
    ht_250 = interplevel(z, p, 250.)
    #dbz_850 = interplevel(dbz3, p, 850.)
    u_850 = interplevel(ua, p, 850)
    v_850 = interplevel(va, p, 850)
    u_500 = interplevel(ua, p, 500)
    v_500 = interplevel(va, p, 500)
    u_250 = interplevel(ua, p, 250)
    v_250 = interplevel(va, p, 250)
    wspd_500 = interplevel(wspd, p, 500)
    wspd_250 = interplevel(wspd, p, 250)
    tk_850 = interplevel(tk, p, 850)
    PV_850 = interplevel(pvo, p, 850)

    # For first time through, get the cartopy mapping object from the netcdf file
    if j == 0:
        cart_proj = get_cartopy(wrfin=ncfile)
        print (cart_proj)
        # Get the geobounds from the netcdf file (by default, uses XLAT, XLONG)
        # You can supply a variable name to get the staggered boundaries
        bounds = geo_bounds(wrfin=ncfile)
        # Download and add the states and coastlines
        states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none", name="admin_1_states_provinces_shp")
        # Get the latitude and longitude points
        lats, lons = latlon_coords(slp)
        
    if doSLP == True:
        # Now do some cosmetic set-up to prepare for plotting
        # Smooth the sea level pressure since it tends to be noisy near complex terrain
        # The smoother can be adjusted, experiment with different values
        smooth_slp = smooth2d(slp, 3, cenweight=4)

        # Make reflectivity below 0 or -5 NaN to avoid plotting, get reflectivity color table
        mdbz_plot = np.where(mdbz > 0., mdbz, "NaN")
        ctables.registry.get_colortable('NWSReflectivity')

        # Now do the plotting, first SLP and reflectivity
        # Create figure, adjust size for higher quality or resolution
        fig = plt.figure(figsize=(16,12))

        # Set the GeoAxes to the projection used by WRF, add state boundaries and coastlines
        ax = plt.axes(projection=cart_proj)
        ax.add_feature(states, linewidth=.5, edgecolor="brown")
        ax.coastlines('50m', linewidth=0.8)

        # Make the contour outlines and filled contours for reflectivity.
        plt.contourf(to_np(lons), to_np(lats), to_np(mdbz_plot), 10, transform=crs.PlateCarree(),
             cmap=ctables.registry.get_colortable('NWSReflectivity'), norm=Normalize(0,60), vmin=0, vmax=60, alpha=.5)

        # Plot SLP with contour labels
        cs = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp), levels=range(916,1068,4), colors='black', 
                 linestyles='solid', transform=crs.PlateCarree())
        plt.clabel(cs, fmt= '%.0f', inline = True)

        # Add color bar
        m = plt.cm.ScalarMappable(cmap=ctables.registry.get_colortable('NWSReflectivity'))
        m.set_array(mdbz)
        m.set_clim(0., 60.)
        plt.colorbar(m, shrink=.75, boundaries=np.linspace(0, 60, 13), alpha=.5)

        # Set the map bounds
        ax.set_xlim(cartopy_xlim(smooth_slp))
        ax.set_ylim(cartopy_ylim(smooth_slp))

        # Add the lat/long gridlines (for Mercator projection)
        ax.gridlines(color="white", linestyle="dotted")

        # Set plot title
        plt.title("Sea Level Pressure (hPa) and composite reflectivity,"+' '+titletime+' UTC')

        # Create separate plot file and save, then show and close
        outTPlotName= 'SLP_maxdBZ'+filetime+'.png'
        plt.savefig(plotsdir+outTPlotName)
        # uncomment following line if you want plots to show inline, otherwise just check for files in ./plots
        #plt.show()
        plt.close()
        
    if do850 == True:
        # Now do 850-mb map with height, temperature, and wind
        fig = plt.figure(figsize=(16,12))
        ax = plt.axes(projection=cart_proj)
        ax.add_feature(states, linewidth=0.4, edgecolor="brown")
        ax.coastlines('50m', linewidth=0.)

        # Add the 850 hPa geopotential height contours, why is 850 in m, other levels dam?
        ht_850 = ht_850/10.  #put in dam
        levels = np.arange(100., 175., 3.)
        contours = plt.contour(to_np(lons), to_np(lats), to_np(ht_850),
                       levels=levels, colors="black",
                       transform=crs.PlateCarree())
        plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

        # temperature fill contours
        tc_850 = tk_850 - 273.15
        levels = [-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35 ]
        temp_contours = plt.contourf(to_np(lons), to_np(lats), to_np(tc_850),
                             levels=levels, cmap=get_cmap("rainbow"),
                             transform=crs.PlateCarree(), alpha=.9)
        plt.colorbar(temp_contours, ax=ax, orientation="horizontal", pad=.03, shrink=.8, alpha=.9)

        # Add wind barbs, only plotting every iskip_th data point.
        iskip = 15
        plt.barbs(to_np(lons[::iskip,::iskip]), to_np(lats[::iskip,::iskip]),
          to_np(u_850[::iskip, ::iskip]), to_np(v_850[::iskip, ::iskip]),
          transform=crs.PlateCarree(), length=6, color='gray')

        # Set the map bounds, grid lines
        ax.set_xlim(cartopy_xlim(slp))
        ax.set_ylim(cartopy_ylim(slp))
        ax.gridlines()

        plt.title("850-hPa Height (dam), Temperature (C), Barbs (kt),"+' '+titletime+' UTC')
        # Create separate plot file and save as .png, then show and close
        outTPlotName= '850Z_T'+filetime+'.png'
        plt.savefig(plotsdir+outTPlotName)
        #plt.show()
        plt.close()

    if do500 == True:
        # Now do 500-mb map with height and wind
        fig = plt.figure(figsize=(16,12))
        ax = plt.axes(projection=cart_proj)
        ax.add_feature(states, linewidth=0.4, edgecolor="brown")
        ax.coastlines('50m', linewidth=0.)

        # Add the 500 hPa geopotential height contours, change interval to 6 for standard convention, 3 for summer
        ht_500 = ht_500/10.  #put in dam
        levels = np.arange(498., 606., 3.)
        contours = plt.contour(to_np(lons), to_np(lats), to_np(ht_500),
                       levels=levels, colors="black",
                       transform=crs.PlateCarree())
        plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

        # Wind speed fill contours
        levels = [25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120]
        wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(wspd_500),
                             levels=levels,
                             cmap=get_cmap("rainbow"),
                             transform=crs.PlateCarree())
        plt.colorbar(wspd_contours, ax=ax, orientation="horizontal", pad=.03, shrink=.8)

        # Add 500 hPa wind barbs, only plotting every iskip_th data point.
        iskip = 15
        plt.barbs(to_np(lons[::iskip,::iskip]), to_np(lats[::iskip,::iskip]),
          to_np(u_500[::iskip, ::iskip]), to_np(v_500[::iskip, ::iskip]),
          transform=crs.PlateCarree(), length=6, color='gray')

        # Set the map bounds, grid lines
        ax.set_xlim(cartopy_xlim(slp))
        ax.set_ylim(cartopy_ylim(slp))
        ax.gridlines()

        plt.title("500-hPa Height (dam), Wind Speed (kt), Barbs (kt),"+' '+titletime+' UTC')
        # Create separate plot file and save as .png, then show and close
        outTPlotName= '500Z_wind'+filetime+'.png'
        plt.savefig(plotsdir+outTPlotName)
        #plt.show()
        plt.close()
        
    if doPV == True:
        # Now plot PV map
        fig = plt.figure(figsize=(16,12))
        ax = plt.axes(projection=cart_proj)
        ax.add_feature(states, linewidth=0.4, edgecolor="brown")
        ax.coastlines('50m', linewidth=0.)

        # Add the 850 hPa PV contours
        levels = np.arange(0.5, 5, 0.5)
        contours = plt.contour(to_np(lons), to_np(lats), to_np(PV_850),
                       levels=levels, colors="black", linewidth=.02,
                       transform=crs.PlateCarree())
        plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

        # Add the wind speed contours
        #levels = [40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200]
        wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(PV_850),
                             levels=levels,
                             cmap=get_cmap("rainbow"),
                             transform=crs.PlateCarree())
        plt.colorbar(wspd_contours, ax=ax, orientation="horizontal", pad=.03, shrink=.8, aspect=50)

        # Add the 850 hPa wind barbs, only plotting every iskip_th data point.
        iskip = 15
        plt.barbs(to_np(lons[::iskip,::iskip]), to_np(lats[::iskip,::iskip]),
          to_np(u_850[::iskip, ::iskip]), to_np(v_850[::iskip, ::iskip]),
          transform=crs.PlateCarree(), length=6, color='gray')

        # Set the map bounds, grid lines
        ax.set_xlim(cartopy_xlim(slp))
        ax.set_ylim(cartopy_ylim(slp))
        ax.gridlines()

        plt.title("850-hPa Potential vorticity (PVU), Wind Barbs (kt),"+' '+titletime+' UTC')
        # Create separate plot file and save as .png, then show and close
        outTPlotName= '850PV_'+filetime+'.png'
        plt.savefig(plotsdir+outTPlotName)
        #plt.show()
        plt.close()
        
    if do250 == True:
        # Now plot 250 mb map
        # Create the figure
        fig = plt.figure(figsize=(16,12))
        ax = plt.axes(projection=cart_proj)
        ax.add_feature(states, linewidth=0.4, edgecolor="brown")
        ax.coastlines('50m', linewidth=0.)

        # Add the 250 hPa geopotential height contours, change interval to 6 for standard convention, 3 for summer
        ht_250 = ht_250/10.  #put in dam
        levels = np.arange(948., 1104., 6.)
        contours = plt.contour(to_np(lons), to_np(lats), to_np(ht_250),
                       levels=levels, colors="black",
                       transform=crs.PlateCarree())
        plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

        # Add the wind speed contours
        levels = [40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200]
        wspd_contours = plt.contourf(to_np(lons), to_np(lats), to_np(wspd_250),
                             levels=levels,
                             cmap=get_cmap("rainbow"),
                             transform=crs.PlateCarree())
        plt.colorbar(wspd_contours, ax=ax, orientation="horizontal", pad=.03, shrink=.8, aspect=50)

        # Add the 250 hPa wind barbs, only plotting every iskip_th data point.
        iskip = 18
        plt.barbs(to_np(lons[::iskip,::iskip]), to_np(lats[::iskip,::iskip]),
          to_np(u_250[::iskip, ::iskip]), to_np(v_250[::iskip, ::iskip]),
          transform=crs.PlateCarree(), length=7, color='gray')

        # Set the map bounds, grid lines
        ax.set_xlim(cartopy_xlim(slp))
        ax.set_ylim(cartopy_ylim(slp))
        ax.gridlines()

        plt.title("250-hPa Height (dam), Wind Speed (kt), Barbs (kt),"+' '+titletime+' UTC')
        # Create separate plot file and save as .png, then show and close
        outTPlotName= '250Z_wind'+filetime+'.png'
        plt.savefig(plotsdir+outTPlotName)
        #plt.show()
        plt.close()


# In[ ]:




