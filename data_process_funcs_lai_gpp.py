import numpy as np
import xarray as xr
import pandas as pd
import os
import glob
import dask

def calc_monmean_modis(data_dir, fname_prefix, fname_suffix, start_yr, end_yr, out_dir, varname, start_mon = '01', end_mon = '12'):
    '''
    Function to calculate monthly mean from MODIS data.
    The data is assumed to be saved in separate files by year
    start_mon: string; the data starts from start_mon of start_yr
    '''
    for year in range(start_yr, end_yr+1):
        print('Working on year ' + str(year))
        fname =  fname_prefix + str(year) + fname_suffix
        ds_modis = xr.open_mfdataset(data_dir + fname, chunks = {'time': 1})#chunks = {'lat':800, 'lon':1000})
        da_temp = ds_modis[varname].sel(time = slice(str(year), None)).groupby('time.month').mean('time')

        time_index = pd.date_range(str(year) + "-01", freq="M", periods=12)
        
        if year == start_yr:
            time_index = pd.date_range(str(year) + "-" + start_mon, freq="M", periods=(12-int(start_mon)+1))
        
        if year == end_yr:
            time_index = pd.date_range(str(year) + "-" + end_mon, freq="M", periods=int(end_mon))
        
        da_temp = da_temp.rename({'month':'time'}).assign_coords({'time':time_index})
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        da_temp.encoding['zlib'] = True
        da_temp.encoding['complevel'] = 1
        
        #print(da_temp['Gpp_500m'].encoding)
        # for var in da_temp.data_vars:
        #     da_temp[var].encoding['zlib'] = True
        #     da_temp[var].encoding['complevel'] = 1
        #     #del(da_temp[var].encoding['contiguous'])
        #     #del(da_temp[var].encoding['chunksizes'])
        
        out_file_path = out_dir + fname
        da_temp.to_netcdf(out_file_path)
        
def calc_seasmean_modis(data_dir, fname_prefix, fname_suffix, start_yr, end_yr, out_dir, varname):
    '''
    Function to calculate seasonal mean from MODIS data.
    The data is assumed to be saved in separate files by year
    This function can process only full years of data. 
    The data_dir should also contain contain a file for the year start_yr-1, december data from that file will be used.
    '''
    for year in range(start_yr, end_yr):
        print('Working on year ' + str(year))
        fnames =  [data_dir + fname_prefix + str(year-1) + fname_suffix, data_dir + fname_prefix + str(year) + fname_suffix]
        ds_modis = xr.open_mfdataset(fnames, chunks = {'time': 1})#chunks = {'lat':800, 'lon':1000})
        da_temp = ds_modis[varname].sel(time = slice(str(year-1)+ "-12", str(year)+ "-11")).groupby('time.season').mean('time')
        # time_index = pd.date_range(str(year) + "-01", freq="3MS", periods=4). The season are in order DJF, JJA, MAM, SON - so can't assign this directly.

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        da_temp.encoding['zlib'] = True
        da_temp.encoding['complevel'] = 1
        
        out_file_path = out_dir + fname_prefix + str(year) + fname_suffix
        da_temp.to_netcdf(out_file_path)
        
def combine_terra_aqua(terra_dir, terra_prefix, terra_suffix, 
                       aqua_dir, aqua_prefix, aqua_suffix, 
                       start_yr, end_yr, 
                       out_dir, out_prefix, out_suffix, varname = None, units = None):
    '''
    Function to calculate mean of estimates from TERRA and AQUA (Used for GPP; LAI already has a combined product). 
    The data is assumed to be saved in separate files by year
    varname 
    '''
    for year in range(start_yr, end_yr):
        print('Working on year ' + str(year))
        fname_terra =  terra_dir + terra_prefix + str(year) + terra_suffix
        ds_terra = xr.open_mfdataset(fname_terra)
        fname_aqua =  aqua_dir + aqua_prefix + str(year) + aqua_suffix
        ds_aqua = xr.open_mfdataset(fname_aqua)
        
        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            ds_comb = xr.concat([ds_terra, ds_aqua], dim = 'datasource').assign_coords({'datasource': ['terra', 'aqua']})
            ds_comb_mean = ds_comb.mean('datasource')

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            
        for var in ds_comb_mean.data_vars:
            ds_comb_mean[var].encoding['zlib'] = True
            ds_comb_mean[var].encoding['complevel'] = 1
            
        if varname is not None:
            ds_comb_mean[varname].assign_attrs({'units': units})
            
        out_file_path = out_dir + out_prefix + str(year) + out_suffix
        ds_comb_mean.to_netcdf(out_file_path)
        
def remove_fillvals(in_dir, fname_prefix, varname, fillthresh, out_dir):
    '''
    Function to remove fill values above a threshold
    Missed this while processsing the mon & seas mean data
    '''
    fnames = glob.glob(in_dir + fname_prefix + '*.nc')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for f in fnames:
        ds = xr.open_dataset(f)
        da_nofill = ds[varname].where(ds[varname] < fillthresh)
        da_nofill.encoding['zlib'] = True
        da_nofill.encoding['complevel'] = 1
        da_nofill.to_netcdf(out_dir + f.split('/')[-1])

def calc_overallMean(data_dir, fname_prefix, fname_suffix, varname, start_yr, end_yr, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    ds_sample = xr.open_dataset(data_dir + fname_prefix + start_yr + fname_suffix)
    weights = np.cos(np.deg2rad(ds_sample.lat))
    weights.name = "weights"
    
    for year in range(start_yr, end_yr+1):
        fname = data_dir + fname_prefix + str(year) + fname_suffix
        ds = xr.open_dataset(fname)
        da_mean = ds[varname].weighted(weights).mean(['lat', 'lon'])
        
        out_file = out_dir + 'overallMean_' + fname_prefix + str(year) + fname_suffix
        da_mean.encoding['zlib'] = True
        da_mean.encoding['complevel'] = 1
        da_mean.to_netcdf(out_file)

def calc_regMean(data_dir, fname_prefix, fname_suffix, varname, start_yr, end_yr, out_dir,
                 reg_file = '/g/data/w97/ad9701/drought_2017to2020/Huthcinson_vegetation_cover_map_6_classes_modis500m.nc',
                 reg_varname = 'land_cover',
                 reg_name = ['tropics', 'savanna', 'warm_temperate', 'cool_temperate', 'mediterranean', 'desert']):
    '''
    Function to calculate regional mean using a regional netcdf file
    reg_file: string; ncfile containing the regions. This file marks the regions using sequential numbers. The order should correspond to the region sequence in the reg file.
    reg_varname: string; the name of the variable in reg_file 
    reg_name: list; the names of the regions to be used for naming the dimension in the output file. 
              If set to None, the reg_nos taken from the reg_file will be used for naming the dimension.
    '''         
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    ds_reg = xr.open_dataset(reg_file)
    ds_sample = xr.open_dataset(data_dir + fname_prefix + start_yr + fname_suffix)
    weights = np.cos(np.deg2rad(ds_sample.lat))
    weights.name = "weights"

    for year in range(start_yr, end_yr+1):
        fname = data_dir + fname_prefix + str(year) + fname_suffix
        ds = xr.open_dataset(fname)

        list_regMean = []
        reg_nos_temp = np.unique(ds_reg[reg_varname].values)
        reg_nos_temp.sort()
        reg_nos = reg_nos_temp[~np.isnan(reg_nos_temp)]
        
        for ireg in reg_nos:
            da_ireg = ds[varname].where(ds_reg[reg_varname] == ireg).weighted(weights).mean(['lat', 'lon'])
            list_regMean.append(da_ireg)
        if reg_name is None:
            da_regMean = xr.concat(list_regMean, dim = 'region').assign_coords({'region': reg_nos})
        else:
            da_regMean = xr.concat(list_regMean, dim = 'region').assign_coords({'region': reg_name})
        out_file = out_dir + 'regMean_' + fname_prefix + str(year) + fname_suffix
        
        da_regMean.encoding['zlib'] = True
        da_regMean.encoding['complevel'] = 1
        da_regMean.to_netcdf(out_file)
        
if __name__ == "main":
    
    # for LAI (Quality control already done)
    data_dir = '/g/data/w97/ad9701/MODIS/MCD15A2H.006/after_QC/'
    out_dir = '/g/data/w97/ad9701/MODIS/MCD15A2H.006/after_QC/mon_mean/'
    fname_prefix = 'MCD15A2H.006_500m_aid0001_'
    fname_suffix = '_QC.nc'
    start_yr = 2002
    end_yr = 2021
    # calc_monmean_modis
    # calc_seasmean_modis
    
    # for GPP
    for prod_name in ['MOD17A2HGF.061', 'MYD17A2HGF.061']:
        data_dir = '/g/data/w97/ad9701/MODIS/' + prod_name + '/'
        out_dir = data_dir + '/mon_mean/'
        fname_prefix = prod_name + '_500m_aid0001_'
        fname_suffix = '.nc'
        start_yr = 2003
        end_yr = 2020
        varname = 'Gpp_500m'
    # calc_monmean_modis
    # calc_seasmean_modis
    
    # an exmaple of removing the fill values (same for MYD17A2HGF.061)
    in_dir = '/g/data/w97/ad9701/MODIS/MOD17A2HGF.061/mon_mean_temp/'
    out_dir = '/g/data/w97/ad9701/MODIS/MOD17A2HGF.061/mon_mean/'
    fname_prefix = 'MOD17A2HGF.061'
    varname = 'Gpp_500m'
    fillthresh = 3.2   # values above this threshold are various types of fill values. So subset for values less than this threshold.
    # remove_fillVals
    # Gpp_FILLVALUE_DOC :
    # FILL VALUE LEGEND
    # 32767 : _Fillvalue: not-computed or outside projection...
    # 32766 : water (ocean or inland)
    # 32765 : barren, very sparsely vegetated
    # 32764 : perennial snow,ice on pixel
    # 32763 : permanant wetlands,marshes
    # 32762 : urban,built-up 
    # 32761 : unclassified
    
    
    
