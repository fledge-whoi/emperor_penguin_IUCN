# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy.interpolate
#import metpy.calc as mpcalc
import math
import netCDF4
import cmocean
import pop_tools 
import matplotlib.colors as mcolors
import matplotlib as mpl



# ----------------------- CHOOSE PARAMTERS--------------------------------

# CHOOSE yend
yend1 = 636 + 5*12 #2073
yend2 = 1020 - (5*12) #2095

model = 3 #1=CMR, 2=IPM, 3=SAT, 4=3 gathered

N_ind = 10 # can be 10 or 100

# ------------------------------------------------------------------------


# --------------------- CALCULATE MEAN OF ENSEMBLES ----------------------

# Change directory to Codes_EP
ds = xr.open_mfdataset('/Users/aliceeparvier/Desktop/Codes_EP/Figures/Fig_july/Maps_python/data_for_AliceP/IFRAC_CESM2LE_SSP370_member1.nc')

ice=ds.IFRAC
ystart=108 #2024

## Calculate the variation on 2024-2073 or 2024-2100
# 10-year average centered on ystart and yend to capture a decadal trend

######## 10-year mean for 2024 - 2073
SICa1= (ice[yend1-(5*12):yend1+12+(5*12)].mean(dim='time') - ice[ystart-(5*12):ystart+12+(5*12)].mean(dim='time')) / ice[yend1-(5*12):yend1+(5*12)+12].mean(dim='time')
SIC_tot1=SICa1
    
######## 10-year mean for 2024 - 2100
SICa2= (ice[yend2-(5*12):yend2+12+(5*12) ].mean(dim='time') - ice[ystart-(5*12):ystart+12+(5*12)].mean(dim='time')) / ice[yend2-(5*12):yend2+(5*12)+12].mean(dim='time')        
SIC_tot2=SICa2


for ens in range(2,51): #climate ensembles
    ds = xr.open_mfdataset(f'/Users/aliceeparvier/Desktop/Codes_EP/Figures/Fig_july/Maps_python/data_for_AliceP/IFRAC_CESM2LE_SSP370_member{ens}.nc')
    ice=ds.IFRAC #(time, long, lat)
        
    # 2024 - 2073
    SICa1= (ice[yend1-(5*12):yend1+12+(5*12)].mean(dim='time') - ice[ystart-(5*12):ystart+12+(5*12)].mean(dim='time')) / ice[yend1-(5*12):yend1+(5*12)+12].mean(dim='time')
    SICa1 = SICa1*100 # percentage
    SIC_tot1=SIC_tot1 + SICa1
    
    # 2024-2100
    SICa2= (ice[yend2-(5*12):yend2+12+(5*12)].mean(dim='time') - ice[ystart-(5*12):ystart+12+(5*12)].mean(dim='time')) / ice[yend2-(5*12):yend2+(5*12)+12].mean(dim='time')
    SICa2 = SICa2*100 # percentage
    SIC_tot2=SIC_tot2 + SICa2

SIC_tot1 = SIC_tot1/50 #mean
SIC_tot2 = SIC_tot2/50 #mean

lon=ds.TLONG #longitude
lat=ds.TLAT #latitude
    
    
# ------------------------------ FIGURE ---------------------------------


i=1
fig =plt.figure(figsize=(20,30))

for i in range(1,7):

    if i==1:
        ax = fig.add_subplot(321, projection=ccrs.SouthPolarStereo())
        ax.set_title('SPCMR 2073', fontweight='bold', size=30)
        ax.set_title('a', loc="left", fontweight='bold', size=34, fontfamily='sans-serif', pad=17)
    elif i==2:
        ax = fig.add_subplot(322, projection=ccrs.SouthPolarStereo())
        ax.set_title('SPCMR 2100', fontweight='bold', size=30)
        ax.set_title('b', loc="left", fontweight='bold', size=34, fontfamily='sans-serif', pad=17)
    elif i==3:
        ax = fig.add_subplot(323, projection=ccrs.SouthPolarStereo())
        ax.set_title('SPIPM 2073', fontweight='bold', size=30)
        ax.set_title('c', loc="left", fontweight='bold', size=34, fontfamily='sans-serif', pad=17)
    elif i==4:
        ax = fig.add_subplot(324, projection=ccrs.SouthPolarStereo())
        ax.set_title('SPIPM 2100', fontweight='bold', size=30)
        ax.set_title('d', loc="left", fontweight='bold', size=34, fontfamily='sans-serif', pad=17)
    elif i==5:
        ax = fig.add_subplot(325, projection=ccrs.SouthPolarStereo())
        ax.set_title('SCDSAT 2073', fontweight='bold', size=30)
        ax.set_title('e', loc="left", fontweight='bold', size=34, fontfamily='sans-serif', pad=17)
    else:
        ax = fig.add_subplot(326, projection=ccrs.SouthPolarStereo())
        ax.set_title('SCDSAT 2100', fontweight='bold', size=30)
        ax.set_title('f', loc="left", fontweight='bold', size=34, fontfamily='sans-serif', pad=17)
        
    ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
    ax.coastlines('10m',linewidth=1.5)
    ax.add_feature(cartopy.feature.LAND, color='#F6F1D9')
    ax.set_facecolor('#D0D0D0')
    


    ##### ADD THE ICE DATA

    # Center white on 0
    norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=-50, vmax=10)
    
    cmap = plt.colormaps['PRGn']

    if i==1 or i==3 or i==5:
        pc = ax.pcolormesh(lon, lat, SIC_tot1, norm=norm, transform=ccrs.PlateCarree(), cmap=cmap)
    else:
        pc = ax.pcolormesh(lon, lat, SIC_tot2, norm=norm, transform=ccrs.PlateCarree(), cmap=cmap)
    #cbar1 = fig.colorbar(pc, extend='both')
    #cbar1.set_label(label='Ice concentration difference (%)', size=20)
    #cbar1.ax.tick_params(labelsize=20) 


    cmap2= plt.cm.jet

    cmaplist=[cmap2(i) for i in range(60, cmap2.N)]
    cmap2=mpl.colors.LinearSegmentedColormap.from_list('Custom map', cmaplist, cmap2.N)
    bounds=np.linspace(0,100,11)
    norm=mpl.colors.BoundaryNorm(bounds, cmap2.N)
    



    ##### Import status of colonies, observed sizes in 2009-2018

    #csv_file = '/Users/aliceeparvier/Desktop/Codes_EP/Figures/Figure_Pext-T/Pext_col_multiLE.csv'
    couleur = ['#00AE1A', '#FFD800', '#FF8700', '#FF0009', '#AE092B']

    if i==1 or i==2:
        csv_file = '/Users/aliceeparvier/Desktop/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_critE/Pext_col_CMR_SSP370_critE.csv'
    elif i==3 or i==4:
        csv_file = '/Users/aliceeparvier/Desktop/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_critE/Pext_col_IPM_SSP370_critE.csv'
    elif i==5 or i==6:
        csv_file = '/Users/aliceeparvier/Desktop/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_critE/Pext_col_SAT_SSP370_critE.csv'
    else: #ecological models gathered
        csv_file = '/Users/aliceeparvier/Desktop/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_critE/Pext_col_ecomod.csv'


    EP_df = pd.read_csv(csv_file, encoding = "ISO-8859-1")
    #EP_df
    df2 = EP_df[['Long', 'Lat', 'Pext_10_2073', 'Pext_100_2073', 'Pext_10_2100', 'Pext_100_2100', 'N']]
    lats = df2.Lat
    lons = df2.Long
    Pext_100_2073=df2.Pext_100_2073
    Pext_100_2100=df2.Pext_100_2100
    Pext_10_2073=df2.Pext_10_2073
    Pext_10_2100=df2.Pext_10_2100
    
    if i==1 or i==3 or i==5:
        Pext_100=Pext_100_2073
        Pext_10=Pext_10_2073
    else:
        Pext_100=Pext_100_2100
        Pext_10=Pext_10_2100
        
            
    if N_ind==10:
        Pext=Pext_10
    else:
        Pext=Pext_100


    # Colony size
    N = df2.N

    for col in range(1, 66):
        if N[col] < 1000:
            taille=70
        elif N[col] >= 1000 and N[col]<5000:
            taille=180
        elif N[col] >= 5000 and N[col]<10000:
            taille=330
        elif N[col] >= 10000 and N[col]<15000:
            taille=480
        elif N[col] >= 15000 and N[col]<20000:
            taille=630
        elif N[col] >= 20000 and N[col]<25000:
            taille=830
        elif N[col] >= 1000 and N[col]<5000:
            taille=1130
        elif N[col] > 25000:
            taille=1530
            
    
        #sc = ax.scatter(lons[col], lats[col], marker='o', c=couleur[status[col-1] -1], vmin=0, vmax =1, s=taille, edgecolors='k', linewidth=2, transform=ccrs.PlateCarree(), zorder=3)
        sc = ax.scatter(lons[col], lats[col], marker='o', c=Pext[col]*100, s=taille, edgecolors='k', linewidth=2, transform=ccrs.PlateCarree(), zorder=3, cmap=cmap2, norm=norm)

        
    i=i+1 
    #end of loop on panels

fig.subplots_adjust(right=0.8, bottom=0.2, top=0.9, left=0.1)
cbar_ax=fig.add_axes([0.85, 0.4, 0.03, 0.3])

cbar2 = fig.colorbar(sc, extend='both', cax=cbar_ax)
cbar2.set_label(label='Extinction probability (%)', size=34, fontweight='bold')
cbar2.ax.tick_params(labelsize=30) 



cbar_ax=fig.add_axes([0.20, 0.15, 0.5, 0.019])

cbar1 = fig.colorbar(pc, extend='both', cax=cbar_ax, orientation='horizontal')
cbar1.set_label(label='Ice concentration difference (%)', size=34, fontweight='bold')
cbar1.ax.tick_params(labelsize=30) 


# Save figure
# 2073

# if N_ind==10:
    # plt.savefig('Fig_Pext_10.eps', format='eps')
# else:
    # plt.savefig('Fig_Pext_100.eps', format='eps')


