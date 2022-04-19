# read data
#===============================================================================

library(ncdf4) # package for netcdf manipulation
library(ncdf4.helpers)
library(SPEI)

outfile = "SPEI3_gleam_monthly_1981_2020.nc"
infile = "PminusPET_gleam_monthly_1980_2020.nc"

# outfile = "SPEI3_awra_monthly_1911_2020.nc"
# infile = "PminusPET_awra_monthly_1911_2020.nc"

infile <- nc_open(infile)

varname = "PminusPET"
variable = ncvar_get(infile, varname)
time = ncvar_get(infile, "time")
time_atts = ncatt_get(infile, "time")
lat = ncvar_get(infile, "lat")
lon = ncvar_get(infile, "lon")
Time_in_Date_format <- nc.get.time.series(f = infile,
                                          time.dim.name = "time")
nc_close(infile)

# calculate SPEI3
#===============================================================================

year = 1900 + as.POSIXlt(Time_in_Date_format)$year
mon = 1 + as.POSIXlt(Time_in_Date_format)$mon

spei3_gridded = array(NA, dim(variable))
for (lon_i in 1:length(lon)) {
  for (lat_i in 1:length(lat)) {
    if (all(!is.na(variable[lon_i, lat_i, ]))){
      ts_pt <- ts(variable[lon_i, lat_i, ], start=c(year[1], mon[1]), end=c(year[length(year)], mon[length(mon)]), frequency=12)
      spei3_pt = SPEI::spei(ts_pt, scale = 3, ref.start=c(1981,1), ref.end=c(2020,5))
      spei3_gridded[lon_i, lat_i, ] = spei3_pt$fitted
    }
  }
}


# write to output netcdf file
#===============================================================================

#----------------
# Make dimensions
#----------------
londim <- ncdim_def( name = 'lon', units = 'degrees_east', vals = lon, longname = "Longitude")
latdim <- ncdim_def( name = 'lat', units = 'degrees_north', vals = lat, longname = "Latitude" )
timedim <- ncdim_def( name = 'time', units = time_atts$units, calendar = time_atts$calendar, vals = time)

#---------
# Make var
#---------
spei3_ncvar <- ncvar_def( name = 'SPEI3', units = 'standardised units', list(londim, latdim, timedim), missval = NA, prec="double")

#---------------------
# Make new output file
#---------------------
ncid_out <- nc_create( outfile, vars = spei3_ncvar) #, force_v4 = TRUE)

#-------------------------------
# Put some test data in the file
#-------------------------------
ncvar_put( nc = ncid_out, varid = spei3_ncvar, vals = spei3_gridded, start=NA, count=NA, verbose=TRUE )
nc_close(ncid_out)
