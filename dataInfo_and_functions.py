mswep_info = dict(full_path = ['/g/data/fj4/SatellitePrecip/MSWEP_V280/Past/Monthly/', '/g/data/fj4/SatellitePrecip/MSWEP_V280/NRT/Monthly/'], 
                  file_name = '*.nc',
                  varname = 'precipitation',
                  lat_slice = slice(-10, -44),
                  lon_slice = slice(112, 154),
                  units = 'mm month-1',
                  rename_latlon = False,
                  land_mask = False,
                   mon_file = 'precipitation_mswep_monthly_1979_2021.nc')
chirps_info = dict(full_path='/g/data/w97/ad9701/CHIRPS-2.0/global_monthly/netcdf/',
                    file_name = 'chirps-v2.0.monthly.nc',
                    varname = 'precip',
                    lat_slice = slice(-44, -10),
                    lon_slice = slice(112, 154),
                    units = 'mm month-1',
                    rename_latlon = True,
                    land_mask = True,
                     mon_file = 'chirps_monthly_1981_2021.nc')
agcd_info = dict(full_path = '/g/data/zv2/agcd/v1/precip/total/r005/01month/',
                file_name = 'agcd_v1_precip_total_r005_monthly_*.nc',
                varname = 'precip',
                lat_slice = slice(-44, -10),
                lon_slice = slice(112, 154),
                units = 'mm month-1',
                rename_latlon = False,
                land_mask = False,
                mon_file = 'agcd_monthly_1900_2020.nc')
# Using PET to calculate SPEI, so variable is set to PET for now
gleam_info = dict(full_path = '/g/data/ua8/GLEAM_v3-5/v3-5a/monthly/', 
                  file_name = 'Ep_1980-2020_GLEAM_v3.5a_MO.nc',
                  varname = 'Ep',
                  lat_slice = slice(-10, -44),
                  lon_slice = slice(112, 154),
                  units = 'mm month-1',
                  rename_latlon = False,
                  mon_file = 'PminusPET_gleam_monthly_1980_2020.nc')
awra_info = dict(full_path = '/g/data/fj8/BoM/AWRA/DATA/SCHEDULED-V6/', 
                  file_name = 'e0_*.nc',
                  varname = 'e0',
                  lat_slice = slice(-10, -44),
                  lon_slice = slice(112, 154),
                  units = 'mm month-1',
                  rename_latlon = True,
                  mon_file = 'PminusPET_awra_monthly_1911_2020.nc')
esacci_info = dict(full_path = ['/g/data/w97/ad9701/ESACCI/3_unzip/' + str(i) + '/' for i in range(1979, 2021)],
                   file_name = 'ESACCI-SOILMOISTURE*.nc',
                   varname = 'sm',
                   lat_slice = slice(-10, -44),
                   lon_slice = slice(112, 154),
                   units = 'm3 m-3',
                   rename_latlon = False,
                   mon_file = 'esacci_monthly_1979_2020.nc')

alldata_dict = dict(mswep = mswep_info, chirps = chirps_info, agcd = agcd_info, gleam = gleam_info, awra = awra_info, esacci = esacci_info)

import xarray as xr
import numpy as np
import sys
import glob

def save_monthly_data_byyear(alldata_dict, data_name, out_dir, calc_from_daily = False, calc_fun = 'sum'):
    '''
    Used in cases where there are too may files that need to be processed directory where each directory is a year.
    '''
    if type(data_name) == list:
        data_list = data_name
    elif type(data_name) == str:
        data_list = [data_name]
    else:
        sys.exit("data_name should be a list or a string")
        
    latlon_rename = {'latitude': 'lat', 'longitude': 'lon'}
        
    for d in data_list:
        for path in alldata_dict[d]['full_path']:
            print('working on directory ' + str(path))
            da_AU = get_da_simple(path,
                                  alldata_dict[d]['file_name'],
                                  alldata_dict[d]['varname'],
                                  alldata_dict[d]['rename_latlon'],
                                  alldata_dict[d]['lat_slice'],
                                  alldata_dict[d]['lon_slice'])                     
            if calc_from_daily:
                if calc_fun == 'sum':
                    da_AU_mon = da_AU.resample(time="M").sum()
                elif calc_fun == 'mean':
                    da_AU_mon = da_AU.resample(time="M").mean()
                out_file = alldata_dict[d]['varname'] + '_' + d + '_monthly_' + str(da_AU_mon['time.year'].min().values) + '.nc'
                da_AU_mon.to_netcdf(out_dir + out_file)
            else:
                out_file = alldata_dict[d]['varname'] + '_' + d + '_monthly_' + str(da_AU['time.year'].min().values) + '.nc'
                da_AU.to_netcdf(out_dir + out_file)
    return None

def get_da_simple(full_path, file_name, varname, rename_latlon, lat_slice, lon_slice):
    latlon_rename = {'latitude': 'lat', 'longitude': 'lon'}
    star_loc = file_name.find('*')
    if star_loc == -1:
        ds = xr.open_dataset(full_path + file_name)
    else:
        ds = xr.open_mfdataset(full_path + file_name)

    if rename_latlon:
        da_AU = ds[varname].rename(latlon_rename).sel(lat = lat_slice, lon = lon_slice)
    else:
        da_AU = ds[varname].sel(lat = lat_slice, lon = lon_slice)
    return da_AU

def save_monthly_data(alldata_dict, data_name, out_dir, calc_from_daily = False, calc_fun = 'sum'):
    if type(data_name) == list:
        data_list = data_name
    elif type(data_name) == str:
        data_list = [data_name]
    else:
        sys.exit("data_name should be a list or a string")
        
    latlon_rename = {'latitude': 'lat', 'longitude': 'lon'}
        
    for d in data_list:
        da_AU = get_da(alldata_dict, d)
    if calc_from_daily:
        if calc_fun == 'sum':
            da_AU_mon = da_AU.resample(time="M").sum()
        elif calc_fun == 'mean':
            da_AU_mon = da_AU.resample(time="M").mean()
        print('Calculated monthly data')
        out_file = alldata_dict[d]['varname'] + '_' + d + '_monthly_' + str(da_AU_mon['time.year'].min().values) + '_' + str(da_AU_mon['time.year'].max().values) + '.nc'
        da_AU_mon.to_netcdf(out_dir + out_file)
    else:
        out_file = alldata_dict[d]['varname'] + '_' + d + '_monthly_' + str(da_AU['time.year'].min().values) + '_' + str(da_AU['time.year'].max().values) + '.nc'
        da_AU.to_netcdf(out_dir + out_file)
    return None

def get_da(alldata_dict, d):
    latlon_rename = {'latitude': 'lat', 'longitude': 'lon'}
    fnames = None
    if type(alldata_dict[d]['full_path']) == list:    # data has to be read from multiple locations
        fnames = []
        for p in alldata_dict[d]['full_path']:
            fnames.extend(glob.glob(p + alldata_dict[d]['file_name']))
        # print(fnames)
                
    if fnames is None:
        star_loc = alldata_dict[d]['file_name'].find('*')
        if star_loc == -1:
            ds = xr.open_dataset(alldata_dict[d]['full_path'] + alldata_dict[d]['file_name'])
        else:
            ds = xr.open_mfdataset(alldata_dict[d]['full_path'] + alldata_dict[d]['file_name'])
            # print('completed reading daily files')
    else:
        ds = xr.open_mfdataset(fnames)
    if alldata_dict[d]['rename_latlon']:
        da_AU = ds[alldata_dict[d]['varname']].rename(latlon_rename).sel(lat = alldata_dict[d]['lat_slice'], lon = alldata_dict[d]['lon_slice'])
    else:
        da_AU = ds[alldata_dict[d]['varname']].sel(lat = alldata_dict[d]['lat_slice'], lon = alldata_dict[d]['lon_slice'])
    # print('returning daily data to main func')
    return da_AU

# regridding to common resolution for averaging
import xesmf as xe

def regrid_all_from_list(da_list, lat = np.arange(-10.125, -44.125, -0.25), lon = np.arange(112.125, 154.125, 0.25)):
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], lat),
            "lon": (["lon"], lon),
        }
    )
    shape_tuple = (len(lat), len(lon))    
    
    da_list_regrid = []  
    for da in da_list:
        da = da.chunk({'lat':-1, 'lon':-1})
        if (da.values.shape[1:] == shape_tuple) | (da.values.shape == shape_tuple):     # assumed that same shape means data is already at the desired resolution
            da_list_regrid.append(da)
        else:
            regridder = xe.Regridder(da, ds_out, 'bilinear')
            da_reg = regridder(da)
            da_list_regrid.append(da_reg)        
    return da_list_regrid

import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point

def draw_spatial_plot_3panels(ds_list, cmap, levels, subplot_title, main_title, out_dir, out_figname, add_cbar = True, cbar_label = ''):
    ds = ds_list[0]
    fig, axs = plt.subplots(nrows=1,ncols=3,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(24,9)) #width, height

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat

    for i in np.arange(len(ds_list)):
        cs=axs[i].contourf(ds_list[i]['lon'],ds_list[i]['lat'],ds_list[i],levels,
                              transform = ccrs.PlateCarree(),
                              cmap=cmap) #,extend='both')   #cmap options: coolwarm,

        # Draw the coastines for each subplot
        axs[i].coastlines()
        axs[i].add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')

        axs[i].set_title(subplot_title[i])

        # Longitude labels
        axs[i].set_xticks(xticks, crs=ccrs.PlateCarree())
        lon_formatter = cticker.LongitudeFormatter()
        axs[i].xaxis.set_major_formatter(lon_formatter)
        axs[i].set_xlim(xlim)

        # Latitude labels
        axs[i].set_yticks(yticks, crs=ccrs.PlateCarree())
        lat_formatter = cticker.LatitudeFormatter()
        axs[i].yaxis.set_major_formatter(lat_formatter)
        axs[i].set_ylim(ylim)

    # Delete the unwanted axes
    # for i in [5]:
    #     fig.delaxes(axs[i])

    # # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.3, top=0.85, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.3, 0.15, 0.4, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)

    plt.suptitle(main_title)
    plt.savefig(out_dir + out_figname)

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

def draw_spatial_plot(ds, cmap, levels, main_title=None, out_dir=None, out_figname=None, add_cbar = True, cbar_label = ''):
    fig, axs = plt.subplots(nrows=1,ncols=1,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(12,9)) #width, height

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat
    
    # xticks = np.arange(113,156,1)  #lon
    # yticks = np.arange(-40,-10,1)   #lat

    cs=axs.contourf(ds['lon'],ds['lat'],ds,levels,
                          transform = ccrs.PlateCarree(),
                          cmap=cmap) #,extend='both')   #cmap options: coolwarm,

    # Draw the coastines for each subplot
    axs.coastlines()
    axs.add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
    if main_title is not None:
        axs.set_title(main_title)

    gl = axs.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right=True
    gl.xlines = True
    gl.ylines = True
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
 
########################################################
#     # Longitude labels
#     axs.set_xticks(xticks, crs=ccrs.PlateCarree())
#     lon_formatter = cticker.LongitudeFormatter()
#     axs.xaxis.set_major_formatter(lon_formatter)
#     axs.set_xlim(xlim)

#     # Latitude labels
#     axs.set_yticks(yticks, crs=ccrs.PlateCarree())
#     lat_formatter = cticker.LatitudeFormatter()
#     axs.yaxis.set_major_formatter(lat_formatter)
#     axs.set_ylim(ylim)

    # Delete the unwanted axes
    # for i in [5]:
    #     fig.delaxes(axs[i])

    # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.3, top=0.95, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.15, 0.15, 0.7, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)
        
    if out_dir is not None:
        if out_figname is not None:
            plt.savefig(out_dir + out_figname)
        else:
            plt.savefig(out_dir + 'figure.png')


def draw_spatial_plot_addcontours(ds, cmap, levels, main_title, out_dir, out_figname, 
                                 ds_contour_list, contour_level, contour_labels, contour_colors, add_cbar = True, cbar_label = ''):
    fig, axs = plt.subplots(nrows=1,ncols=1,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(12,9)) #width, height

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat

    cs=axs.contourf(ds['lon'],ds['lat'],ds,levels,
                          transform = ccrs.PlateCarree(),
                          cmap=cmap) #,extend='both')   #cmap options: coolwarm,
    
    contour_elements_list = []
    for i in range(len(ds_contour_list)):
        cs_contour=axs.contour(ds_contour_list[i]['lon'],ds_contour_list[i]['lat'],ds_contour_list[i],contour_level,
                    transform = ccrs.PlateCarree(),
                    colors=contour_colors[i]) #,extend='both')   #cmap options: coolwarm,
        h1,_ = cs_contour.legend_elements()
        contour_elements_list.append(h1[0])
    
    axs.legend(contour_elements_list, contour_labels, loc='lower left')
    # plt.clabel(cs_contour, inline=1, fontsize=10)
    # cs_contour.collections[0].set_label(contour_label)
    # plt.legend(loc='lower left')
    
    # Draw the coastines for each subplot
    axs.coastlines()
    axs.add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
    axs.set_title(main_title)

    # Longitude labels
    axs.set_xticks(xticks, crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    axs.xaxis.set_major_formatter(lon_formatter)
    axs.set_xlim(xlim)

    # Latitude labels
    axs.set_yticks(yticks, crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    axs.yaxis.set_major_formatter(lat_formatter)
    axs.set_ylim(ylim)

    # Delete the unwanted axes
    # for i in [5]:
    #     fig.delaxes(axs[i])

    # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.3, top=0.95, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.15, 0.15, 0.7, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)

    plt.savefig(out_dir + out_figname)
    
def draw_spatial_plot_12panels(ds_list, cmap, levels, subplot_title=None, main_title=None, out_dir=None, out_figname=None, 
                               add_cbar = True, cbar_label = ''): #, unwanted_axes=None):
    npanels = len(ds_list)
    if npanels < 12:
        unwanted_axes = list(range(npanels, 12))
    else:
        unwanted_axes = None
        
    ds = ds_list[0]
    fig, axs = plt.subplots(nrows=2,ncols=6,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(22,8)) #width, height
    axs = axs.flatten()

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat

    for i in np.arange(len(ds_list)):
        #print(i)
        cs=axs[i].contourf(ds_list[i]['lon'],ds_list[i]['lat'],ds_list[i],levels,
                              transform = ccrs.PlateCarree(),
                              cmap=cmap, extend='both')   #cmap options: coolwarm,

        # Draw the coastines for each subplot
        axs[i].coastlines()
        axs[i].add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
        if subplot_title is not None:
            axs[i].set_title(subplot_title[i], pad = 2)
        
        axs[i].set_ylim(ylim)
        axs[i].set_xlim(xlim)
        
        if i in [6, 7, 8, 9, 10, 11]:
            # Longitude labels
            axs[i].set_xticks(xticks, crs=ccrs.PlateCarree())
            lon_formatter = cticker.LongitudeFormatter()
            axs[i].xaxis.set_major_formatter(lon_formatter)

        if i in [0, 6]:
            # Latitude labels
            axs[i].set_yticks(yticks, crs=ccrs.PlateCarree())
            lat_formatter = cticker.LatitudeFormatter()
            axs[i].yaxis.set_major_formatter(lat_formatter)

    # Delete the unwanted axes
    if unwanted_axes is not None:
        for i in unwanted_axes:
            fig.delaxes(axs[i])

    # # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.25, top=0.9, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.3, 0.15, 0.4, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)

    if main_title is not None:
        plt.suptitle(main_title)
    if out_dir is not None:
        if out_figname is not None:
            plt.savefig(out_dir + out_figname)
        else:
            plt.savefig(out_dir + 'figure.png')
   
import geopandas as gpd

def draw_spatial_plot_addsh(ds, cmap, levels, main_title=None, out_dir=None, out_figname=None, add_cbar = True, cbar_label = '', shapefile=None):
    fig, axs = plt.subplots(nrows=1,ncols=1,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(12,9)) #width, height

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat
    
    cs=axs.contourf(ds['lon'],ds['lat'],ds,levels,
                          transform = ccrs.PlateCarree(),
                          cmap=cmap) #,extend='both')   #cmap options: coolwarm,

    # Draw the coastines for each subplot
    axs.coastlines()
    axs.add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
    if main_title is not None:
        axs.set_title(main_title)

    # xticks = np.arange(113,156,1)  #lon
    # yticks = np.arange(-40,-10,1)   #lat
    # gl = axs.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
    # gl.xlabels_top = False
    # gl.ylabels_left = False
    # gl.ylabels_right=True
    # gl.xlines = True
    # gl.ylines = True
    # gl.xlocator = mticker.FixedLocator(xticks)
    # gl.ylocator = mticker.FixedLocator(yticks)
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER
 
########################################################
    # Longitude labels
    axs.set_xticks(xticks, crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    axs.xaxis.set_major_formatter(lon_formatter)
    axs.set_xlim(xlim)

    # Latitude labels
    axs.set_yticks(yticks, crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    axs.yaxis.set_major_formatter(lat_formatter)
    axs.set_ylim(ylim)
    
    if shapefile is not None:
        shdf = gpd.read_file(shapefile)
        #shdf = shdf.to_crs("EPSG:4326")
        projection = ccrs.PlateCarree() #central_longitude=0) 
        axs.add_geometries(shdf.geometry,
                  projection,
                  facecolor='none',
                  edgecolor='red', zorder=3)
        
        

    # Delete the unwanted axes
    # for i in [5]:
    #     fig.delaxes(axs[i])

    # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.25, top=0.95, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.15, 0.15, 0.7, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)
        
    if out_dir is not None:
        if out_figname is not None:
            plt.savefig(out_dir + out_figname)
        else:
            plt.savefig(out_dir + 'figure.png')