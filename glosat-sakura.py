#!/usr/bin/env python

#-----------------------------------------------------------------------
# PROGRAM: glosat-sakura.py
#-----------------------------------------------------------------------
# Version 0.1
# 28 March, 2021
# Dr Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
#-----------------------------------------------------------------------

# Dataframe libraries:
import numpy as np
import pandas as pd
import xarray as xr

# Stats libraries:
from loess import loess

# Datetime libraries:
from datetime import datetime
import nc_time_axis
import cftime

# Plotting libraries:
import matplotlib
#matplotlib.use('agg')
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import image
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

# Maths libraries:
from scipy.interpolate import griddata
from scipy import spatial
from math import radians, cos, sin, asin, sqrt

# OS libraries:
import os, sys
from  optparse import OptionParser
#import argparse

# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

#----------------------------------------------------------------------------
# DARK BACKGROUND THEME
#----------------------------------------------------------------------------

matplotlib.rcParams['text.usetex'] = False
rcParams['font.family'] = ['DejaVu Sans']
rcParams['font.sans-serif'] = ['Avant Garde']
plt.rc('text',color='white')
plt.rc('lines',color='white')
plt.rc('patch',edgecolor='white')
plt.rc('grid',color='lightgray')
plt.rc('xtick',color='white')
plt.rc('ytick',color='white')
plt.rc('axes',edgecolor='lightgray')
plt.rc('axes',facecolor='black')
plt.rc('axes',labelcolor='white')
plt.rc('figure',facecolor='black')
plt.rc('figure',edgecolor='black')
plt.rc('savefig',edgecolor='black')
plt.rc('savefig',facecolor='black')

# Calculate current time

now = datetime.now()
currentmn = str(now.month)
if now.day == 1:
    currentdy = str(cal.monthrange(now.year,now.month-1)[1])
    currentmn = str(now.month-1)
else:
    currentdy = str(now.day-1)
if int(currentdy) < 10:
    currentdy = '0' + currentdy    
currentyr = str(now.year)
if int(currentmn) < 10:
    currentmn = '0' + currentmn
titletime = str(currentdy) + '/' + currentmn + '/' + currentyr

#---------------
#-----------------------------------------------------------------------------
# SETTINGS
#-----------------------------------------------------------------------------

fontsize = 16
use_anomalies = False
station_code = '477590' # Kyoto City
station_latlon = 'lon_135.7_lat_35.0'
img = image.imread("IMAGES/sakura4.jpeg")

#-----------------------------------------------------------------------------
# LOAD: GloSAT timseries
#-----------------------------------------------------------------------------

if use_anomalies == True:
    df_in = pd.read_pickle('DATA/df_anom.pkl', compression='bz2')        
    dirstr = 'DATA/ANOMALIES/' + station_latlon
else:
    df_in = pd.read_pickle('DATA/df_temp.pkl', compression='bz2')            
    dirstr = 'DATA/ABSOLUTES/' + station_latlon                
df = df_in[df_in['stationcode']==station_code].reset_index(drop=True)

t_glosat = pd.date_range(start=str(df['year'][0]), periods=len(df), freq='A')   
ts_glosat = np.mean(np.array(df.groupby('year').mean().iloc[:,0:11]),axis=1) 

#-----------------------------------------------------------------------------
# LOAD: 20CRv3 gridcell timseries
#-----------------------------------------------------------------------------

ensemble_mean = dirstr+'/'+'ensemble_mean'+'_'+station_latlon+'.csv'
ensemble_sd = dirstr+'/ensemble_sd'+'_'+station_latlon+'.csv'
ensemble_min = dirstr+'/ensemble_min'+'_'+station_latlon+'.csv'
ensemble_max = dirstr+'/ensemble_max'+'_'+station_latlon+'.csv'
ensemble_pctl_05 = dirstr+'/ensemble_pctl_05'+'_'+station_latlon+'.csv'
ensemble_pctl_10 = dirstr+'/ensemble_pctl_10'+'_'+station_latlon+'.csv'
ensemble_pctl_25 = dirstr+'/ensemble_pctl_25'+'_'+station_latlon+'.csv'
ensemble_pctl_50 = dirstr+'/ensemble_pctl_50'+'_'+station_latlon+'.csv'
ensemble_pctl_75 = dirstr+'/ensemble_pctl_75'+'_'+station_latlon+'.csv'
ensemble_pctl_90 = dirstr+'/ensemble_pctl_90'+'_'+station_latlon+'.csv'
ensemble_pctl_95 = dirstr+'/ensemble_pctl_95'+'_'+station_latlon+'.csv'

#-----------------------------------------------------------------------------
# LOAD: Sakura DOY timeseries
#-----------------------------------------------------------------------------

# STNNo. A.D.     FiFD FuFD WORK TYPE Name of reference
# 47759 2008   PJ 0328 0404    6    0 NEWS-PAPER(ARASHIYAMA)

sakura_file = 'DATA/aono-cherry-blossom-kyoto.txt'
nheader = 26
f = open(sakura_file)
lines = f.readlines()
years = []
FiFD_months = []
FuFD_months = []
FiFD_days = []
FuFD_days = []
for i in range(nheader,len(lines)):
        chars = lines[i]   
        year = chars[6:10].strip()
        FiFD = chars[16:20].strip()
        if not FiFD: 
            FiFD_month = '999'
            FiFD_day = '999'
        else: 
            FiFD_month = FiFD[0:2]
            FiFD_day = FiFD[2:4]
        FuFD = chars[21:26].strip()
        if not FuFD: 
            FuFD_month = '999'
            FuFD_day = '999'
        else: 
            FuFD_month = FuFD[0:2]
            FuFD_day = FuFD[2:4]
        years.append(year)        
        FiFD_months.append(FiFD_month)        
        FuFD_months.append(FuFD_month)        
        FiFD_days.append(FiFD_day)        
        FuFD_days.append(FuFD_day)        
f.close()    

FiFD_DOYs = []
FuFD_DOYs = []
for i in range(len(years)):
    if FiFD_months[i] != "" and FiFD_months[i] != '999':
        date = cftime.DatetimeGregorian(int(years[i]), int(FiFD_months[i]), int(FiFD_days[i]))
        year = cftime.DatetimeGregorian(int(years[i]), 1, 1)
        doy = (date-year).days
    else:
        doy = np.nan
    FiFD_DOYs.append(doy)
    if FuFD_months[i] != "" and FuFD_months[i] != '999':
        date = cftime.DatetimeGregorian(int(years[i]), int(FuFD_months[i]), int(FuFD_days[i]))
        year = cftime.DatetimeGregorian(int(years[i]), 1, 1)
        doy = (date-year).days
    else:
        doy = np.nan
    FuFD_DOYs.append(doy)
        
t_sakura = xr.cftime_range(start=years[0].zfill(4), periods=len(years), freq='A', calendar='gregorian')     
ts_sakura_Fi = FiFD_DOYs
ts_sakura_Fu = FuFD_DOYs

df_sakura = pd.DataFrame()
df_sakura['t_sakura'] = t_sakura
df_sakura['ts_sakura_Fi'] = ts_sakura_Fi
df_sakura['ts_sakura_Fu'] = ts_sakura_Fu

#-----------------------------------------------------------------------------
# LOAD: Pages2K (Ed Hawkins)
#-----------------------------------------------------------------------------

# 1945	0.014	0.006	-0.2553	0.2013	-0.0692	-0.0692	-0.3131	0.0901
# Year CE | raw instrumental target data | reconstruction ensemble 50th | 2.5th | 97.5th percentiles | 31-year butterworth filtered instrumental target data | 31-year butterworth filtered reconstruction 50th | 2.5th | 97.5th percentiles

pages2k_file = 'DATA/Full_ensemble_median_and 95pct_range.txt'
nheader = 5
f = open(pages2k_file)
lines = f.readlines()
years = []
obs = []
for i in range(nheader,len(lines)):
        words = lines[i].split()   
        year = words[0].zfill(4)
        val = (len(words)-1)*[None]            
        for j in range(len(val)):                                
            try: val[j] = float(words[j+1])                
            except: 
                pass                                 
        years.append(year)                                     
        obs.append(val)            
f.close()    
obs = np.array(obs)

t_pages2k = xr.cftime_range(start=years[0], periods=len(years), freq='A', calendar='gregorian')[800:]
ts_pages2k_instr = pd.to_numeric(obs[:,0][800:], errors='coerce')
ts_pages2k_recon = pd.to_numeric(obs[:,5][800:], errors='coerce')
ts_pages2k = np.append(ts_pages2k_recon[0:-36],ts_pages2k_instr[-36:],axis=None)

df_pages2k = pd.DataFrame()
df_pages2k['t_pages2k'] = t_pages2k
df_pages2k['ts_pages2k'] = ts_pages2k

#-----------------------------------------------------------------------------
# STATS: LOESS fit
#-----------------------------------------------------------------------------

if use_loess == True:

    df_sakura_nonan_Fu = df_sakura.dropna( how='any', subset=['ts_sakura_Fu'])
    df_sakura_nonan_Fi = df_sakura[df_sakura['t_sakura'].index>1000].dropna( how='any', subset=['ts_sakura_Fi'])
    regsDF_Fu, evalDF_Fu = loess(df_sakura_nonan_Fu['t_sakura'].index, df_sakura_nonan_Fu['ts_sakura_Fu'], alpha=.5, poly_degree=2)
    regsDF_Fi, evalDF_Fi = loess(df_sakura_nonan_Fi['t_sakura'].index, df_sakura_nonan_Fi['ts_sakura_Fi'], alpha=.9, poly_degree=1)
    #l_x  = evalDF['v'].values
    l_x_Fu  = df_sakura_nonan_Fu['t_sakura']
    l_x_Fi  = df_sakura_nonan_Fi['t_sakura']
    l_y_Fu  = evalDF_Fu['g'].values[1:]
    l_y_Fi  = evalDF_Fi['g'].values[1:]

#-----------------------------------------------------------------------------
# PLOT: GloSAT temperatures and Sakura (Cherry Blossom) timeseries on same xaxis
#-----------------------------------------------------------------------------


figstr = 'glosat-sakura.png'
#titlestr = 'GloSAT temperatures at Kyoto and Sakura timeseries'
titlestr = 'Phenological data for flowering dates of Prunus jamasakura in Kyoto City'

fig,ax = plt.subplots(figsize=(15,10))

plt.plot(df_sakura['t_sakura'],df_sakura['ts_sakura_Fu'], color='cyan', marker='o', linestyle='None', alpha=0.4, label='Prunus jamasakura (full flowering)')
plt.plot(df_sakura['t_sakura'],df_sakura['ts_sakura_Fi'], color='red', marker='s', linestyle='None', alpha=0.4, label='Prunus jamasakura (first flowering)')
plt.plot(l_x_Fu, l_y_Fu, color='cyan', linestyle='-', linewidth=5, alpha=1.0, label='LOESS fit: '+r'$\alpha$=0.5' + ', quadratic')
plt.plot(l_x_Fi, l_y_Fi, color='red', linestyle='-', linewidth=5, alpha=1.0, label='LOESS fit: '+r'$\alpha$=0.1' + ', linear')
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
aspect = img.shape[0] / img.shape[1] * (xmax - xmin)/(ymax - ymin)
plt.imshow(img, zorder=0, extent=[xmin, xmax, ymin, ymax], aspect=aspect, alpha=0.8)
ax.tick_params(labelsize=fontsize)
ax.xaxis.grid(True, which='major')        
ax.yaxis.grid(True, which='major')      
#plt.ylim(40,140)  
leg = plt.legend(loc=2, ncol=1, fontsize=fontsize)
leg.set_frame_on(False) # make it transparent
ax.set_xlabel("Year A.D.", fontsize=fontsize)
ax.set_ylabel("Day of the Year", fontsize=fontsize)    

ax2 = ax.twinx()
ax2.plot(df_pages2k['t_pages2k'], df_pages2k['ts_pages2k'], color='pink', marker='o', linestyle='-', alpha=0.4, label=r'Temperature anomaly (from 1991-1990) [$^{\circ}$C]')
ax2.spines['right'].set_color('pink')
ax2.tick_params(axis='y', color='pink', labelsize=fontsize)
ax2.set_ylabel(r'Temperature anomaly (from 1961-1990) [$^{\circ}$C]', color='pink', fontsize=fontsize)    

#if use_anomalies == True:
#    plt.ylabel("Temperature anomaly, $\mathrm{\degree}C$", fontsize=fontsize)
#else:
#    plt.ylabel("Absolute temperature, $\mathrm{\degree}C$", fontsize=fontsize)    
#plt.title(titlestr, fontsize=fontsize)

datastr1 = r'$\bf{Dataset 1}$' + ': Aono and Kazui, 2008; Aono and Saito, 2010; Aono, 2012'        
sourcestr1 = r'$\bf{Source 1}$' + ': http://atmenv.envi.osakafu-u.ac.jp/aono/kyophenotemp4/'        
datastr2 = r'$\bf{Dataset 2}$' + ': PAGES2k (& HadCRUT4.6.0.0 for 2001-2017): courtesy of Ed Hawkins, University of Reading)'        
sourcestr2 = r'$\bf{Source 2}$' + ': http://www.climate-lab-book.ac.uk/2020/2019-years/'        
authorstr = r'$\bf{Graphic}$' + ': Michael Taylor, CRU/UEA' + ' -- ' + titletime

fig.suptitle(titlestr, fontsize=24, color='white', fontweight='bold')        
plt.annotate(datastr1, xy=(100,125), xycoords='figure pixels', color='white', fontsize=fontsize) 
plt.annotate(sourcestr1, xy=(100,100), xycoords='figure pixels', color='white', fontsize=fontsize) 
plt.annotate(datastr2, xy=(100,75), xycoords='figure pixels', color='white', fontsize=fontsize) 
plt.annotate(sourcestr2, xy=(100,50), xycoords='figure pixels', color='white', fontsize=fontsize) 
plt.annotate(authorstr, xy=(100,20), xycoords='figure pixels', color='white', fontsize=fontsize, bbox=dict(boxstyle="square, pad=0.3", fc='black', edgecolor='white', linewidth=0.2))    
fig.subplots_adjust(left=None, bottom=0.3, right=None, top=None, wspace=None, hspace=None)
plt.savefig(figstr)
plt.close('all')

# -----------------------------------------------------------------------------
print('** END')

