import numpy as np
import xarray as xr
import climtas
import os
import datetime
import itertools

if __name__ == '__main__':
    
    data_dir = '/g/data/w97/ad9701/MODIS/MCD15A2H.006/'
    out_dir = '/g/data/w97/ad9701/MODIS/MCD15A2H.006/after_QC/'
    fname_prefix = 'MCD15A2H.006_500m_aid0001_'

    for year in range(2003, 2022):
        fname =  fname_prefix + str(year) + '.nc'
        ds_lai = xr.open_dataset(data_dir + fname)
        ds_lai['FparLai_QC'].load()
        ds_lai['Lai_500m'].load()
        da_lai_fnl = ds_lai['Lai_500m'].where(ds_lai['Lai_500m']<24.9).where(ds_lai['FparLai_QC'] <= 2)

        out_file = out_dir + fname_prefix + str(year) + '_QC.nc'
        da_lai_fnl.to_netcdf(out_file)