#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 11:43:00 2022

@author: yue
"""

#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from maskout import shp2clip

#open the netCDF file and give it a label
bc_08 = Dataset('/Users/yue/Desktop/project/data/2008/ECLIPSE_V6b_CLE_monthly_anthropogenic_BC_emissions_2008_timeslice.nc','r')
bc_16 = Dataset('/Users/yue/Desktop/project/data/2016/ECLIPSE_V6b_CLE_monthly_anthropogenic_BC_emissions_2016_timeslice.nc','r')
#print(bc_08.variables)

#extract data
lat = bc_08.variables['latitude'][:]
lon = bc_08.variables['longitude'][:]
time= bc_08.variables['time'][:]        #shape(12,)

var = ['agricultural_waste_burning','residential__commercial__other','energy','industrial','international_shipping','transportation','waste','venting_and_flaring']
total = 0
#use for loops to extract all bc emission data of different departments 
for i in var:
    value_08 = bc_08.variables[i]
    value_16 = bc_16.variables[i]
    #calculate the annual mean value
    avvar_08 = np.mean(value_08,axis=0)
    avvar_16 = np.mean(value_16,axis=0)
    dvar = avvar_16-avvar_08
    #calculate the total annual bc emissions from all eight departments
    total +=dvar
    
    #olot figures
    fig = plt.figure(figsize=(10,8))
    proj = ccrs.PlateCarree()
    extent = [70,140,15,55]
    #shape file of China's map
    shp_path = '/Users/yue/Desktop/project/code/cnmap/province.shp'
    read = Reader(shp_path)
    shpf = cfeature.ShapelyFeature(read.geometries(), proj, edgecolor='k',facecolor='none')
    ax = fig.subplots(1,1,subplot_kw={'projection':proj})

    ax.add_feature(shpf)

    ax.set_extent(extent,crs = proj)
    
    #set levels of the data
    levels = [-10**-9,-10**-10,-10**-11,-10**-12,-10**-13,0,10**-13,10**-12,10**-11,10**-10,10**-9]
    cmap1 = plt.get_cmap('RdBu_r')
    norm = colors.BoundaryNorm(levels,cmap1.N,clip=True)
    cf = ax.contourf(lon,lat,dvar,levels=levels, norm=norm,cmap=cmap1,transform=proj)
    #delete the '_' in variables name and capital letters
    varname = ' '.join([name.capitalize() for name in i.split('_')])
    #print(varname)
   
    # Add a title
    ax.set_title('Difference of Annual '+varname+' BC Emission in China Between 2008 and 2016',fontsize=10)
    cbar = plt.colorbar(cf,  orientation='horizontal',fraction=0.06,pad=0.06)

    cbar.set_label('emission(kgm$^{-2}$s$^{-1}$)')
    
    # create the province list of China
    list = [110000,120000,130000,140000,150000,210000,220000,230000,310000,320000,330000,340000,350000,360000,370000,410000,420000,430000,
            440000,450000,460000,500000,510000,520000,530000,540000,610000,620000,630000,640000,650000,710000,810000,820000]
    clip = shp2clip(cf,ax,'/Users/yue/Desktop/project/code/cnmap/province.shp',list)
    plt.savefig('/Users/yue/Desktop/project/figures/bcdiff1/'+ varname +'.png',dpi=200)
    plt.close()

#plot total difference of bc emissions from eight departments in 2008 and 2016
fig = plt.figure(figsize=(10,8))
proj = ccrs.PlateCarree()#等经纬度投影
extent = [70,140,15,55]
read = Reader(shp_path)
shpf = cfeature.ShapelyFeature(read.geometries(), proj, edgecolor='k',facecolor='none')
ax = fig.subplots(1,1,subplot_kw={'projection':proj})

ax.add_feature(shpf)

ax.set_extent(extent,crs = proj) #设置经纬度范围


levels = [-10**-9,-10**-10,-10**-11,-10**-12,-10**-13,0,10**-13,10**-12,10**-11,10**-10,10**-9]
cmap1 = plt.get_cmap('RdBu_r')
norm = colors.BoundaryNorm(levels,cmap1.N,clip=True)
cf = ax.contourf(lon,lat,total,levels=levels, norm=norm,cmap=cmap1,transform=proj)
ax.set_title('Difference of Total Annual BC Emission in China Between 2008 and 2016',fontsize=10)
cbar = plt.colorbar(cf,  orientation='horizontal',fraction=0.06,pad=0.06)

cbar.set_label('emission(kgm$^{-2}$s$^{-1}$)')

list = [110000,120000,130000,140000,150000,210000,220000,230000,310000,320000,330000,340000,350000,360000,370000,410000,420000,430000,
         440000,450000,460000,500000,510000,520000,530000,540000,610000,620000,630000,640000,650000,710000,810000,820000]
clip = shp2clip(cf,ax,'/Users/yue/Desktop/project/code/cnmap/province.shp',list)
plt.savefig('/Users/yue/Desktop/project/figures/bcdiff1/total.png',dpi=200)

