#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 13:34:04 2024

@author: gruma-r
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import os
import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
import geopandas as gpd
from datetime import datetime

data_directory = '/home/gruma-r/Área de Trabalho/nc_era5_data/vort_geop'
files = os.listdir(data_directory)


netcdf_files = [f for f in files if os.path.isfile(os.path.join(data_directory, f))]
netcdf_files.sort()

geop_var = 'z'
vort_var = 'vo'

#Shapefile path
shapefile_path = '/home/gruma-r/Área de Trabalho/Lim_america_do_sul_2021.shp'
# Load Brazilian states shapefile
brazil_states = gpd.read_file(shapefile_path)
brazil_feature = ShapelyFeature(brazil_states['geometry'], ccrs.PlateCarree(), edgecolor='grey', facecolor='none')

# Folder where the figures will be saved
output_folder = '/home/gruma-r/Imagens/monique-pessoal/cpam'  # Change this to your desired output folder
os.makedirs(output_folder, exist_ok=True)

start_date_range = datetime(2024, 1, 15)
end_date_range = datetime(2024, 1, 17)
print(start_date_range)
print(end_date_range)

for file_name in netcdf_files:
    file_path = os.path.join(data_directory, file_name)
   
    with nc.Dataset(file_path, 'r') as netcdf_file:
        latitude = netcdf_file.variables['latitude'][:]
        longitude = netcdf_file.variables['longitude'][:]
        geopotential_data = netcdf_file.variables[geop_var]
        vorticity_data = netcdf_file.variables[vort_var]
        time_variable = netcdf_file.variables['time']

        # Convert units to datetime
        time_values = time_variable[:]
        time_units = time_variable.units
        calendar = 'gregorian'  # Ensure 'gregorian' calendar is used
        custom_time = nc.num2date(time_values, units=time_units, calendar=calendar)
        
        start_index = np.argmin(np.abs(custom_time - start_date_range))
        end_index = np.argmin(np.abs(custom_time - end_date_range))

        # Define gravitational force constant
        g = 98.0665  # m/s^2 (corrected value)
        
        for time_step in range(start_index, end_index + 1):
            z_data_timestep = geopotential_data[time_step, :, :] # Extract all pressure levels
            vorticity_data_timestep = vorticity_data[time_step, :, :]
            lon, lat = np.meshgrid(longitude, latitude)
                        
            # Format the datetime information for the filename
            formatted_time = custom_time[time_step].strftime("%d-%m-%Y_%H:%M:%S")

            # Find the index of the pressure level closest to 500 hPa
            pressure_levels = netcdf_file.variables['level'][:]  # Assuming 'level' is the name of the pressure levels variable
            idx_500hPa = np.abs(pressure_levels - 500).argmin()

            # Extract data at 500 hPa level
            geopotential_height_500hPa = z_data_timestep[idx_500hPa, :, :] / g
            vorticity_500hPa = vorticity_data_timestep[idx_500hPa, :, :]
            
            # Plot the data using Cartopy
            plt.title(f'Vorticidade e geopotencial às {formatted_time}')
            plt.figure(figsize=(10, 8))
            ax = plt.axes(projection=ccrs.PlateCarree())
            #ax.set_extent([longitude.min(), longitude.max(), latitude.min(), latitude.max()])
            ax.set_extent([-30, -90, -20, -55])
            
            #Gridlines (lat e lon)
            gl = ax.gridlines(draw_labels=True, linewidth=0)
            gl.top_labels = False
            gl.right_labels = False
            gl.xlocator = plt.MultipleLocator(15)
            gl.ylocator = plt.MultipleLocator(15)

            cs = ax.contour(lon, lat, geopotential_height_500hPa, levels=10, zorder=3, colors='black')
            plt.clabel(cs, fmt='%d')
            
            vort_arange = np.arange(-12,14,2)
            vorticity_plot = ax.contourf(lon, lat, vorticity_500hPa* 1e5, vort_arange[vort_arange != 0], cmap='RdBu_r', extend='both', zorder=1, transform=ccrs.PlateCarree())
            colorbar = plt.colorbar(vorticity_plot, ax=ax, orientation='horizontal', ticks=vort_arange, label='Relative Vorticity (10$^{-5}$ s$^{-1}$)', shrink = 0.76, pad = 0.08)

            # Draw coastlines and country/state borders
            ax.coastlines(linewidth=0.5, color='grey', zorder=2)
            ax.add_feature(brazil_feature, zorder=2)

            # Construct the output filename
            output_filename = os.path.join(output_folder, f'geopotential_height_500hPa_{formatted_time}.png')

            # Save the figure
            plt.savefig(output_filename, format='png', bbox_inches='tight')
