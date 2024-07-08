# standard
from collections import OrderedDict
import datetime
import glob
import os
import sys
import time
# 3rd party
from configobj import ConfigObj
import netCDF4
import numpy
import cftime
import pytz
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline
import xlrd
# check the scripts directory is present
if not os.path.exists("../scripts/"):
    # CME: if not os.path.exists("../../../PyFluxPro/scripts/"):
    print("era52nc: the scripts directory is missing")
    sys.exit()
# since the scripts directory is there, try importing the modules
sys.path.append('../scripts')
# PFP
import constants as c
import meteorologicalfunctions as mf
import pysolar
import pfp_io
import pfp_log
import pfp_utils
import pfp_compliance
"""
This script reads hourly ERA5 data and picks out the tower site data as a timeline 
it then converts these to tower time step using an Akima univariate spline interpolation.
Akima univariante spline interpolation insures the spline is going through the data point.
"""
now = datetime.datetime.now()
log_file_name = 'process_era52nc_' + now.strftime("%Y%m%d%H%M") + ".log"
log_file_name = os.path.join("./logfiles", log_file_name)
logger = pfp_log.init_logger("pfp_log", log_file_name, to_file=True, to_screen=True)

def read_site_master(xl_file_path, sheet_name):
    """
    Reads the site master file with entries of site name and site location (latitude, longitude).
    site_master.xls can have more than one sheet. In era52nc.txt the sheetname is specified.
    """
    xl_book = xlrd.open_workbook(xl_file_path)
    xl_sheet = xl_book.sheet_by_name(sheet_name)
    last_row = int(xl_sheet.nrows)
    # find the header and first data rows
    for i in range(last_row):
        if xl_sheet.cell(i,0).value == "Site":
            header_row = i
            first_data_row = header_row + 1
            break
    # read the header row
    header_row_values = xl_sheet.row_values(header_row)
    # read the site data from the master Excel spreadsheet
    site_info = OrderedDict()
    for n in range(first_data_row,last_row):
        site_name = xl_sheet.cell(n,0).value
        site_name = site_name.replace(" ","")
        site_info[site_name] = OrderedDict()
        for item in header_row_values[1:]:
            i = header_row_values.index(item)
            site_info[site_name][item] = xl_sheet.cell(n,i).value
    return site_info

# read the controlfile file
""" era52nc.txt controlfile
    [Files]
        xl_file_path   = <location for the site_master.xls file for input>/Sites/site_master.xls
        xl_sheet_name  = <sheet_name; e.g. Active>
        erai_path      = <location for the ECMWF downloaded files>/Sites/ERAI/AUS/*.nc
        out_base_path  = <location for the site specific output results>/Sites/
        out_process    = <new [default](process all existing ERAI files)/ append (process only the last file in the list)>
        concat_control = <path for writing controlfile to concatenate files>
    [Options]
        site_sa_limit = 5 #(normalise the ERA-I downwelling shortwave by the solar altitude, set minimum)
"""
if (__name__ == '__main__'):
    # get the control file name
    if len(sys.argv) == 1:
        # not on the command line, so ask the user
        cfg_file_path = input("Enter the control file name: ")
        # exit if nothing selected
        if len(cfg_file_path) == 0:
            sys.exit()
    else:
        # control file name on the command line
        if not os.path.exists(sys.argv[1]):
            # control file doesn't exist
            logger.error("Control file %s does not exist", sys.argv[1])
            sys.exit()
        else:
            cfg_file_path = sys.argv[1]
# read the control file file
"""
The controlfiles organises the input and output data path, the sheet 
from the site_master (e.g. all sites or individual sites) and the sa_limit.
"""
cf = pfp_io.get_controlfilecontents(cfg_file_path, mode="verbose")
xl_file_path        = cf["Files"]["xl_file_path"]
xl_sheet_name       = cf["Files"]["xl_sheet_name"]
era5_path           = cf["Files"]["era5_path"]
out_base_path       = cf["Files"]["out_base_path"]
out_process         = cf["Files"]["out_process"]
concat_control_path = cf["Files"]["concat_control"]
site_sa_limit       = cf["Options"]["site_sa_limit"]
# get the site information from the site master spreadsheet
site_info = read_site_master(xl_file_path, xl_sheet_name)
# get a list of sites
site_list = list(site_info.keys())
# and a list of the ERA5 files to be processed
era5_list = sorted(glob.glob(era5_path))
if len(era5_list) == 0:
    msg = 'no files found in ERA5 directory'
    logger.warning(msg)
    sys.exit()
if out_process == 'append':
    era5_files = era5_list[-1:]
else:
    era5_files = era5_list
# construct a dictionary of concatenation control files, this will be
# used to concatenate the monthly ERA5 files for each site
cf_dict = OrderedDict()
# ===
# Set flag_codes for ERA5 data
# ===
flag_codes = {"era5": 800, "era5t": 810, "akima_int": 820}
# ===
for site_name in site_list:
    # construct the output file path
    out_file_path = os.path.join(out_base_path,site_name,"Data","ERA5",site_name+"_ERA5.nc")
    # initialise the concatenation control file
    cf_dict[site_name] = ConfigObj(indent_type="    ")
    cf_dict[site_name]["Options"] = {"NumberOfDimensions":1,
                                     "MaxGapInterpolate":1,
                                     "FixTimeStepMethod":"round",
                                     "Truncate":"No",
                                     "TruncateThreshold":50,
                                     "SeriesToCheck":[]}
    cf_dict[site_name]["Files"] = {"Out":{"ncFileName":out_file_path},"In":{}}
    # Add sites' ERA5 datafile to the list of files to be concatenated
    # Check out_file_path for existing file <site>_ERA5.nc
    add_exist = 0
    if out_process == 'append':
        if os.path.isfile(out_file_path):
            # rename existing file and set code to append new processed data
            this_day = datetime.datetime.today().strftime('%Y%m%d')
            this_day_file_path = os.path.join(out_base_path,site_name,"Data","ERA5",this_day+site_name+"_ERA5.nc")
            # Here I need a safeguard, check if already renamed and if target file exists
            if os.path.isfile(this_day_file_path):
                logger.warning(this_day_file_path + " exists, "+out_file_path+" NOT renamed")
            else:
                logger.info(this_day_file_path + " does not exist, "+out_file_path+" renamed to "+this_day_file_path)
                os.rename(out_file_path,this_day_file_path)
        else:
            logger.warning(out_file_path + " does not exists, use new for out_process in control file")
        cf_dict[site_name]["Files"]["In"][str(0)] = this_day_file_path
        add_exist = 1

for n, era5_name in enumerate(era5_files):
    logger.info("Processing ERA5 file "+era5_name)
    era5_timestep = 60
    era5_file = netCDF4.Dataset(era5_name)
    latitude = era5_file.variables["latitude"][:]
    longitude = era5_file.variables["longitude"][:]
    lat_resolution = abs(latitude[-1]-latitude[0])/(len(latitude)-1)
    lon_resolution = abs(longitude[-1]-longitude[0])/(len(longitude)-1)
    # get the time and convert to Python datetime object
    era5_time = era5_file.variables["time"][:]
    time_units = getattr(era5_file.variables["time"],"units")
    try:
        dt_era5 = cftime.num2pydate(era5_time,time_units)
    except:
        dt_era5 = cftime.num2date(era5_time,time_units)
    start_date_era5 = dt_era5[0]
    end_date_era5 = dt_era5[-1]
    hour_utc = numpy.array([dt.hour for dt in dt_era5])
    # get the datetime in the middle of the accumulation period
    era5_offset = datetime.timedelta(minutes=float(era5_timestep)/2)
    dt_era5_cor = [x - era5_offset for x in dt_era5]
    # get a series of time, corrected for the offset
    # NOTE: netCDF4.date2num doesn't handle timezone-aware datetimes
    #era5_time_1hr = netCDF4.date2num(dt_era5_cor,time_units)
    era5_time_1hr = cftime.date2num(dt_era5_cor,time_units)
    # make utc_dt timezone aware so we can generate local times later
    dt_era5_utc_cor = [x.replace(tzinfo=pytz.utc) for x in dt_era5_cor]
    #site_list = ["Tumbarumba"]
    # now loop over the sites
    for site_name in site_list:
        # get the output file name
        out_site_path = os.path.join(out_base_path, site_name, "Data", "ERA5")
        if not os.path.exists(out_site_path):
            os.makedirs(out_site_path)
        out_file_name = site_name+"_ERA5_"+start_date_era5.strftime("%Y%m%d")
        out_file_name = out_file_name+"_"+end_date_era5.strftime("%Y%m%d")+".nc"
        out_file_path = os.path.join(out_site_path, out_file_name)
        # get the metadata from the control file
        logger.info("Processing "+site_name)
        # get the metadata from the site master file information
        site_latitude = site_info[site_name]["Latitude"]
        site_longitude = site_info[site_name]["Longitude"]
        site_timezone = site_info[site_name]["Time zone"]
        site_timestep = int(round(float(site_info[site_name]["Time step"])))
        # index of the site in latitude dimension
        site_lat_index = int(((latitude[0]-site_latitude)/lat_resolution)+0.5)
        era5_latitude = latitude[site_lat_index]
        # index of the site in longitude dimension
        #print(site_lat_index,era5_latitude)
        #if site_longitude<0: site_longitude = float(360) + site_longitude
        site_lon_index = int(((site_longitude-longitude[0])/lon_resolution)+0.5)
        #print(site_lon_index,site_longitude)
        era5_longitude = longitude[site_lon_index]
        logger.info("Site coordinates: "+str(site_latitude)+" "+str(site_longitude))
        logger.info("ERA5 grid: "+str(latitude[site_lat_index])+" "+str(longitude[site_lon_index]))
        # get an instance of the Datastructure
        ds_era5 = pfp_io.DataStructure()
        ds_era5.series["DateTime"] = {}
        ds_era5.globalattributes["site_name"] = site_name
        ds_era5.globalattributes["time_zone"] = site_timezone
        ds_era5.globalattributes["latitude"] = site_latitude
        ds_era5.globalattributes["longitude"] = site_longitude
        ds_era5.globalattributes["time_step"] = site_timestep
        ds_era5.globalattributes["sa_limit"] = site_sa_limit
        ds_era5.globalattributes['xl_datemode'] = str(0)
        ds_era5.globalattributes["nc_level"] = "L1"
        ds_era5.globalattributes["data_set"] = "ERA5"
        # get the UTC and local datetime series
        site_tz = pytz.timezone(site_timezone)
        # now we get the datetime series at the tower time step
        tdts = datetime.timedelta(minutes=site_timestep)
        # get the start and end datetimes rounded to the nearest time steps
        # that lie between the first and last times
        start_date = pfp_utils.rounddttots(dt_era5_utc_cor[0],ts=site_timestep)
        if start_date<dt_era5_utc_cor[0]: start_date = start_date+tdts
        end_date = pfp_utils.rounddttots(dt_era5_utc_cor[-1],ts=site_timestep)
        if end_date>dt_era5_utc_cor[-1]: end_date = end_date-tdts
        msg = "Data: "+start_date.strftime("%Y-%m-%d %H:%M")+" UTC to "
        msg = msg+end_date.strftime("%Y-%m-%d %H:%M")+" UTC"
        logger.info(msg)
        #print site_name,end_date,dt_era5_utc_cor[-1]
        # UTC datetime series at the tower time step
        dt_era5_utc_tts = [x for x in pfp_utils.perdelta(start_date,end_date,tdts)]
        # UTC netCDF time series at tower time step for interpolation
        tmp = [x.replace(tzinfo=None) for x in dt_era5_utc_tts]
        #era5_time_tts = netCDF4.date2num(tmp,time_units)
        era5_time_tts = cftime.date2num(tmp,time_units)
        # local datetime series at tower time step
        dt_era5_loc_tts = [x.astimezone(site_tz) for x in dt_era5_utc_tts]
        # NOTE: will have to disable daylight saving at some stage, towers stay on Standard Time
        # PRI hopes that the following line will do this ...
        dt_era5_loc_tts = [x-x.dst() for x in dt_era5_loc_tts]
        # make the datetime series timezone naive and put it in data structure
        dt_era5_loc_tts = [x.replace(tzinfo=None) for x in dt_era5_loc_tts]
        ds_era5.series["DateTime"]["Data"] = dt_era5_loc_tts
        ds_era5.series["DateTime"]["Flag"] = numpy.zeros(len(dt_era5_loc_tts))
        ds_era5.globalattributes["nc_nrecs"] = len(dt_era5_loc_tts)
        ds_era5.globalattributes["start_datetime"] = str(dt_era5_loc_tts[0])
        ds_era5.globalattributes["end_datetime"] = str(dt_era5_loc_tts[-1])
        # get the Excel datetime
        pfp_utils.get_xldatefromdatetime(ds_era5)
        # get the year, month, day, hour, minute and second
        pfp_utils.get_ymdhmsfromdatetime(ds_era5)
        # get the solar altitude, we will use this later to interpolate the ERA 5 solar
        # data from the ERA-5 1 hour time step to the tower time step.
        # NOTE: alt_solar is in degrees
        alt_solar_1hr = numpy.array([pysolar.GetAltitude(era5_latitude,era5_longitude,dt) for dt in dt_era5_utc_cor])
        # get the solar altitude at the tower time step
        alt_solar_tts = numpy.array([pysolar.GetAltitude(era5_latitude,era5_longitude,dt) for dt in dt_era5_utc_tts])
        idx = numpy.where(alt_solar_tts<=0)[0]
        alt_solar_tts[idx] = float(0)

        # check the number of dimensions - expvar added in dimension due to temporary release of data to 5 days delay from previous 3 month dalay
        nDims = len(era5_file.variables["ssrd"].shape)

        # === DOWNWELLING SHORTWAVE Fsd === #
        # Interpolate the 1 hourly accumulated downwelling shortwave to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude] or if temproray data included [time,expvar,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Fsd_3d = era5_file.variables["ssrd"][:,:,:]
            Fsd_accum = Fsd_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Fsd_3d = era5_file.variables["ssrd"][:,:,:,:]
            data0 = Fsd_3d[:,0,site_lat_index,site_lon_index]
            data1 = Fsd_3d[:,1,site_lat_index,site_lon_index]
            Fsd_accum = numpy.ma.where(data0.mask == False,data0,data1)
        #Fsd_3d = era5_file.variables["ssrd"][:,:,:]
        #Fsd_accum = Fsd_3d[:,site_lat_index,site_lon_index]
        # Downwelling shortwave in ERA-5 is a cummulative value since the previous post processing (archiving), so for
        #                                              HRES: accumulations are in the hour ending at the forecast step.
        # Here we convert the cummulative values to 1 hourly values.
        #Fsd_era5_1hr = numpy.ediff1d(Fsd_accum,to_begin=0)
        #idx = numpy.where((hour_utc==7)|(hour_utc==19))[0]
        #Fsd_era5_1hr[idx] = Fsd_accum[idx]
        #Fsd_era5_1hr = Fsd_era5_1hr/(era5_timestep*60)
        Fsd_era5_1hr = Fsd_accum/(era5_timestep*60)
        # normalise the ERA-5 downwelling shortwave by the solar altitude
        # clamp solar altitude to a minimum value to avoid numerical problems
        # when alt_solar is close to 0
        alt_solar_limit = float(site_sa_limit)*numpy.ones(len(alt_solar_1hr))
        sa = numpy.where(alt_solar_1hr<=float(site_sa_limit),alt_solar_limit,alt_solar_1hr)
        coef_1hr = Fsd_era5_1hr/numpy.sin(numpy.deg2rad(sa))

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, coef_1hr) 
        # get the coefficient at the tower time step
        coef_tts = int_fn(era5_time_tts) 

        # ==== old = UnivariateSpline ==== 
        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, coef_1hr, k=1)
        # get the coefficient at the tower time step
        #coef_tts = s(era5_time_tts)

        # get the downwelling solar radiation at the tower time step
        Fsd_era5_tts = coef_tts*numpy.sin(numpy.deg2rad(alt_solar_tts))
        flag = numpy.zeros(len(Fsd_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Down-welling shortwave radiation", 
                                                 group_name = "radiation", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fsd",Fsd_era5_tts,flag,attr)

        # === NET-SHORTWAVE Fn_sw === #
        # Interpolate the 1 hourly accumulated net shortwave to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Fn_sw_3d = era5_file.variables["ssr"][:,:,:]
            Fn_sw_accum = Fn_sw_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Fn_sw_3d = era5_file.variables["ssr"][:,:,:,:]
            data0 = Fn_sw_3d[:,0,site_lat_index,site_lon_index]
            data1 = Fn_sw_3d[:,1,site_lat_index,site_lon_index]
            Fn_sw_accum = numpy.ma.where(data0.mask == False,data0,data1)
        #Fn_sw_3d = era5_file.variables["ssr"][:,:,:]
        #Fn_sw_accum = Fn_sw_3d[:,site_lat_index,site_lon_index]
        # Net shortwave in ERA-5 is a cummulative value that is reset to 0 at 0300 and 1500 UTC.
        # Here we convert the cummulative values to 3 hourly values.
        #Fn_sw_era5_1hr = numpy.ediff1d(Fn_sw_accum,to_begin=0)
        # deal with the reset times at 0300 and 1500
        #idx = numpy.where((hour_utc==3)|(hour_utc==15))[0]
        #Fn_sw_era5_1hr[idx] = Fn_sw_accum[idx]
        # get the average value over the 1 hourly period
        #Fn_sw_era5_1hr = Fn_sw_era5_1hr/(era5_timestep*60)
        Fn_sw_era5_1hr = Fn_sw_accum/(era5_timestep*60)
        # normalise the ERA-5 et shortwave by the solar altitude
        coef_1hr = Fn_sw_era5_1hr/numpy.sin(numpy.deg2rad(sa))

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, coef_1hr) 
        # get the coefficient at the tower time step
        coef_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, coef_1hr, k=1)
        # get the coefficient at the tower time step
        #coef_tts = s(era5_time_tts)

        # get the downwelling solar radiation at the tower time step
        Fn_sw_era5_tts = coef_tts*numpy.sin(numpy.deg2rad(alt_solar_tts))
        flag = numpy.zeros(len(Fn_sw_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Net shortwave radiation", 
                                                 group_name = "radiation", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fn_sw",Fn_sw_era5_tts,flag,attr)

        # === UPWELLING-SHORTWAVE Fsu === #
        Fsu_era5_tts = Fsd_era5_tts - Fn_sw_era5_tts
        flag = numpy.zeros(len(Fsu_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Up-welling shortwave radiation", 
                                                 group_name = "radiation", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fsu",Fsu_era5_tts,flag,attr)

        # === DOWNWELLING LONGWAVE Fld === #
        # Interpolate the 1 hourly accumulated downwelling longwave to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Fld_3d = era5_file.variables["strd"][:,:,:]
            Fld_accum = Fld_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Fld_3d = era5_file.variables["strd"][:,:,:,:]
            data0 = Fld_3d[:,0,site_lat_index,site_lon_index]
            data1 = Fld_3d[:,1,site_lat_index,site_lon_index]
            Fld_accum = numpy.ma.where(data0.mask == False,data0,data1)
        #Fld_3d = era5_file.variables["strd"][:,:,:]
        #Fld_accum = Fld_3d[:,site_lat_index,site_lon_index]
        # Downwelling longwave in ERA-5 is a cummulative value that is reset to 0 at 0300 and 1500 UTC.
        # Here we convert the cummulative values to 1 hourly values.
        #Fld_era5_1hr = numpy.ediff1d(Fld_accum,to_begin=0)
        #idx = numpy.where((hour_utc==7)|(hour_utc==19))[0]
        #Fld_era5_1hr[idx] = Fld_accum[idx]
        #Fld_era5_1hr = Fld_era5_1hr/(era5_timestep*60)
        Fld_era5_1hr = Fld_accum/(era5_timestep*60)

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Fld_era5_1hr) 
        # get the coefficient at the tower time step
        Fld_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, Fld_era5_1hr, k=1)
        # get the downwelling longwave at the tower time step
        #Fld_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(Fld_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Down-welling longwave radiation", 
                                                 group_name = "radiation", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fld",Fld_era5_tts,flag,attr)

        # === NET-LONGWAVE Fn_lw === #
        # Interpolate the 1 hourly accumulated net longwave to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Fn_lw_3d = era5_file.variables["str"][:,:,:]
            Fn_lw_accum = Fn_lw_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Fn_lw_3d = era5_file.variables["str"][:,:,:,:]
            data0 = Fn_lw_3d[:,0,site_lat_index,site_lon_index]
            data1 = Fn_lw_3d[:,1,site_lat_index,site_lon_index]
            Fn_lw_accum = numpy.ma.where(data0.mask == False,data0,data1)
        #Fn_lw_3d = era5_file.variables["str"][:,:,:]
        #Fn_lw_accum = Fn_lw_3d[:,site_lat_index,site_lon_index]
        # Net longwave in ERA-5 is a cummulative value that is reset at 0300 and 1500 UTC.
        # Here we convert the cummulative values to 3 hourly values.
        #Fn_lw_era5_1hr = numpy.ediff1d(Fn_lw_accum,to_begin=0)
        # deal with the reset times at 0300 and 1500
        #idx = numpy.where((hour_utc==3)|(hour_utc==15))[0]
        #Fn_lw_era5_1hr[idx] = Fn_lw_accum[idx]
        # get the average value over the 1 hourly period
        #Fn_lw_era5_1hr = Fn_lw_era5_1hr/(era5_timestep*60)
        Fn_lw_era5_1hr = Fn_lw_accum/(era5_timestep*60)

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Fn_lw_era5_1hr) 
        # get the coefficient at the tower time step
        Fn_lw_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, Fn_lw_era5_1hr, k=1)
        # get the net longwave at the tower time step
        #Fn_lw_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(Fn_lw_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Net longwave radiation", 
                                                 group_name = "radiation", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fn_lw",Fn_lw_era5_tts,flag,attr)

        # === UPWELLING-LONGWAVE Flu === #
        Flu_era5_tts = Fld_era5_tts - Fn_lw_era5_tts
        flag = numpy.zeros(len(Flu_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Up-welling longwave radiation", 
                                                 group_name = "radiation", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Flu",Flu_era5_tts,flag,attr)

        # === NET RADIATION Fn === #
        Fn_era5_tts = (Fsd_era5_tts - Fsu_era5_tts) + (Fld_era5_tts - Flu_era5_tts)
        flag = numpy.zeros(len(Fn_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Net radiation", 
                                                 group_name = "radiation", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fn",Fn_era5_tts,flag,attr)

        # === SENSIBLE HEAT FLUX Fh === #
        # Interpolate the 1 hourly accumulated sensible heat flux to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Fh_3d = era5_file.variables["sshf"][:,:,:]
            Fh_accum = Fh_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Fh_3d = era5_file.variables["sshf"][:,:,:,:]
            data0 = Fh_3d[:,0,site_lat_index,site_lon_index]
            data1 = Fh_3d[:,1,site_lat_index,site_lon_index]
            Fh_accum = numpy.ma.where(data0.mask == False,data0,data1)
        #Fh_3d = era5_file.variables["sshf"][:,:,:]
        #Fh_accum = float(-1)*Fh_3d[:,site_lat_index,site_lon_index]
        # Sensible heat flux in ERA-5 is a cummulative value that is reset at 0300 and 1500 UTC.
        # Here we convert the cummulative values to 3 hourly values.
        #Fh_era5_1hr = numpy.ediff1d(Fh_accum,to_begin=0)
        # deal with the reset times at 0300 and 1500
        #idx = numpy.where((hour_utc==3)|(hour_utc==15))[0]
        #Fh_era5_1hr[idx] = Fh_accum[idx]
        # get the average value over the 1 hourly period
        #Fh_era5_1hr = Fh_era5_1hr/(era5_timestep*60)
        Fh_era5_1hr = Fh_accum/(era5_timestep*60)

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Fh_era5_1hr) 
        # get the coefficient at the tower time step
        Fh_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, Fh_era5_1hr, k=1)
        # get the net longwave at the tower time step
        #Fh_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(Fh_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Sensible heat flux", 
                                                 group_name = "flux", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fh",Fh_era5_tts,flag,attr)

        # === LATENT HEAT FLUX Fh === #
        # Interpolate the 1 hourly accumulated latent heat flux to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Fe_3d = era5_file.variables["slhf"][:,:,:]
            Fe_accum = Fe_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Fe_3d = era5_file.variables["slhf"][:,:,:,:]
            data0 = Fe_3d[:,0,site_lat_index,site_lon_index]
            data1 = Fe_3d[:,1,site_lat_index,site_lon_index]
            Fe_accum = numpy.ma.where(data0.mask == False,data0,data1)
        #Fe_3d = era5_file.variables["slhf"][:,:,:]
        #Fe_accum = float(-1)*Fe_3d[:,site_lat_index,site_lon_index]
        # Latent heat flux in ERA-5 is a cummulative value that is reset at 0300 and 1500 UTC.
        # Here we convert the cummulative values to 3 hourly values.
        #Fe_era5_1hr = numpy.ediff1d(Fe_accum,to_begin=0)
        # deal with the reset times at 0300 and 1500
        #idx = numpy.where((hour_utc==3)|(hour_utc==15))[0]
        #Fe_era5_1hr[idx] = Fe_accum[idx]
        # get the average value over the 1 hourly period
        #Fe_era5_1hr = Fe_era5_1hr/(era5_timestep*60)
        Fe_era5_1hr = Fe_accum/(era5_timestep*60)

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Fe_era5_1hr) 
        # get the coefficient at the tower time step
        Fe_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, Fe_era5_1hr, k=1)
        # get the net longwave at the tower time step
        #Fe_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(Fe_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Latent heat flux", 
                                                 group_name = "flux", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fe",Fe_era5_tts,flag,attr)

        # === GROUND HEAT FLUX Fg = Fn - Fh - Fe === #
        # get Fg as residual
        Fg_era5_tts = Fn_era5_tts - Fh_era5_tts - Fe_era5_tts
        flag = numpy.zeros(len(Fg_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Ground heat flux", 
                                                 group_name = "soil", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fg",Fg_era5_tts,flag,attr)
        # === AVAILABLE ENERGY Fa = Fn - Fg === #
        # and then Fa
        Fa_era5_tts = Fn_era5_tts - Fg_era5_tts
        flag = numpy.zeros(len(Fa_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Available energy", 
                                                 group_name = "flux", units="W/m^2")
        pfp_utils.CreateSeries(ds_era5,"Fa",Fa_era5_tts,flag,attr)

        # === AIR PRESSURE ps === #
        # Interpolate the 1 hourly air pressure to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            ps_3d = era5_file.variables["sp"][:,:,:]
            ps_era5_1hr = ps_3d[:,site_lat_index,site_lon_index]/float(1000)
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            ps_3d = era5_file.variables["sp"][:,:,:,:]
            data0 = ps_3d[:,0,site_lat_index,site_lon_index]/float(1000)
            data1 = ps_3d[:,1,site_lat_index,site_lon_index]/float(1000)
            ps_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #ps_3d = era5_file.variables["sp"][:,:,:]
        #ps_era5_1hr = ps_3d[:,site_lat_index,site_lon_index]/float(1000)

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, ps_era5_1hr) 
        # get the coefficient at the tower time step
        ps_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, ps_era5_1hr, k=1)
        # get the air pressure at the tower time step
        #ps_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(ps_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Surface air pressure", 
                                                 group_name = "meteorology", units="kPa")
        pfp_utils.CreateSeries(ds_era5,"ps",ps_era5_tts,flag,attr)

        # === TEMPERATURE Ta === #
        # Interpolate the 1 hourly air temperature to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Ta_3d = era5_file.variables["t2m"][:,:,:]
            Ta_era5_1hr = Ta_3d[:,site_lat_index,site_lon_index] - 273.15
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Ta_3d = era5_file.variables["t2m"][:,:,:,:]
            data0 = Ta_3d[:,0,site_lat_index,site_lon_index] - 273.15
            data1 = Ta_3d[:,1,site_lat_index,site_lon_index] - 273.15
            Ta_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Ta_3d = era5_file.variables["t2m"][:,:,:]
        #Ta_era5_1hr = Ta_3d[:,site_lat_index,site_lon_index] - 273.15

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Ta_era5_1hr) 
        # get the coefficient at the tower time step
        Ta_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, Ta_era5_1hr, k=1)
        # get the air temperature at the tower time step
        #Ta_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(Ta_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Air temperature", 
                                                 group_name = "meteorology", units="degC")
        pfp_utils.CreateSeries(ds_era5,"Ta",Ta_era5_tts,flag,attr)

        # === DEWPOINT TEMPERATURE Td === #
        # Interpolate the 1 hourly dew point temperature to the tower time step
        # and convert to Ah, RH and q
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Td_3d = era5_file.variables["d2m"][:,:,:]
            Td_era5_1hr = Td_3d[:,site_lat_index,site_lon_index] - 273.15
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Td_3d = era5_file.variables["d2m"][:,:,:,:]
            data0 = Td_3d[:,0,site_lat_index,site_lon_index] - 273.15
            data1 = Td_3d[:,1,site_lat_index,site_lon_index] - 273.15
            Td_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Td_3d = era5_file.variables["d2m"][:,:,:]
        #Td_era5_1hr = Td_3d[:,site_lat_index,site_lon_index] - 273.15

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Td_era5_1hr) 
        # get the coefficient at the tower time step
        Td_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, Td_era5_1hr, k=1)
        # get the dew point temperature at the towespeedr time step
        #Td_era5_tts = s(era5_time_tts)

        # === RELATIVE HUMIDITY RH === #
        # get the relative humidity
        es_era5_tts = mf.VPsat(Ta_era5_tts)
        e_era5_tts = mf.VPsat(Td_era5_tts)
        VPD_era5_tts = es_era5_tts - e_era5_tts
        flag = numpy.zeros(len(VPD_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Vapour pressure deficit", 
                                                 group_name = "meteorology", units="kPa")
        pfp_utils.CreateSeries(ds_era5,"VPD",VPD_era5_tts,flag,attr)
        RH_era5_tts = float(100)*e_era5_tts/es_era5_tts
        flag = numpy.zeros(len(RH_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Relative humidity", 
                                                 group_name = "meteorology", units="percent")
        pfp_utils.CreateSeries(ds_era5,"RH",RH_era5_tts,flag,attr)

        # === ABSOLUTE HUMIDITY Ah === #
        # get the absolute humidity
        Ah_era5_tts = mf.absolutehumidityfromrelativehumidity(Ta_era5_tts,RH_era5_tts)
        flag = numpy.zeros(len(Ah_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Absolute humidity", 
                                                 group_name = "meteorology", units="g/m^3")
        pfp_utils.CreateSeries(ds_era5,"AH",Ah_era5_tts,flag,attr)

        # === SPECIFIC HUMIDITY Ah === #
        # get the specific humidity
        q_era5_tts = mf.specifichumidityfromRH(RH_era5_tts,Ta_era5_tts,ps_era5_tts)
        flag = numpy.zeros(len(q_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Specific humidity", 
                                                 group_name = "meteorology", units="kg/kg")
        pfp_utils.CreateSeries(ds_era5,"SH",q_era5_tts,flag,attr)

        # === ATMOSPHERIC BOUNDARY LAYER HEIGHT Habl === #
        # Interpolate the 1 hourly boundary layer height to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Habl_3d = era5_file.variables["blh"][:,:,:]
            Habl_era5_1hr = Habl_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Habl_3d = era5_file.variables["blh"][:,:,:,:]
            data0 = Habl_3d[:,0,site_lat_index,site_lon_index]
            data1 = Habl_3d[:,1,site_lat_index,site_lon_index]
            Habl_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Habl_3d = era5_file.variables["blh"][:,:,:]
        #Habl_era5_1hr = Habl_3d[:,site_lat_index,site_lon_index]

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Habl_era5_1hr) 
        # get the coefficient at the tower time step
        Habl_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, Habl_era5_1hr, k=1)
        # get the boundary layer height at the tower time step
        #Habl_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(Habl_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Boundary layer height", 
                                                 group_name = "meteorology", units="m")
        pfp_utils.CreateSeries(ds_era5,"Habl",Habl_era5_tts,flag,attr)

        # === PRECIPITATION Precip === #
        # Spread the 1 hourly accumulated precipitation to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Precip_3d = era5_file.variables["tp"][:,:,:]
            Precip_accum = Precip_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Precip_3d = era5_file.variables["tp"][:,:,:,:]
            data0 = Precip_3d[:,0,site_lat_index,site_lon_index]
            data1 = Precip_3d[:,1,site_lat_index,site_lon_index]
            Precip_accum = numpy.ma.where(data0.mask == False,data0,data1)
        #Precip_3d = era5_file.variables["tp"][:,:,:]
        #Precip_accum = Precip_3d[:,site_lat_index,site_lon_index]
        
        #Precip_era5_1hr = numpy.ediff1d(Precip_accum,to_begin=0)
        #idx = numpy.where((hour_utc==3)|(hour_utc==15))[0]
        #Precip_era5_1hr[idx] = Precip_accum[idx]
        #Precip_era5_1hr = Precip_era5_1hr*float(1000)
        Precip_era5_1hr = Precip_accum*float(1000)
        Precip_era5_tts = numpy.zeros(len(dt_era5_loc_tts))
        idx = pfp_utils.FindIndicesOfBInA(dt_era5_utc_cor,dt_era5_utc_tts)
        if len(idx) != 0:
           Precip_era5_tts[idx] = Precip_era5_1hr
        flag = numpy.zeros(len(Precip_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Precipitation", 
                                                 group_name = "meteorology", units="mm")
        pfp_utils.CreateSeries(ds_era5,"Precip",Precip_era5_tts,flag,attr)

        # === SOIL MOISTURE Sws === #
        # LEVEL1: 
        # Interpolate the 1 hourly soil moisture to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Sws_3d = era5_file.variables["swvl1"][:,:,:]
            Sws_era5_1hr = Sws_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Sws_3d = era5_file.variables["swvl1"][:,:,:,:]
            data0 = Sws_3d[:,0,site_lat_index,site_lon_index]
            data1 = Sws_3d[:,1,site_lat_index,site_lon_index]
            Sws_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Sws_3d = era5_file.variables["swvl1"][:,:,:]
        #Sws_era5_1hr = Sws_3d[:,site_lat_index,site_lon_index]
        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Sws_era5_1hr) 
        # get the coefficient at the tower time step
        Sws_era5_tts = int_fn(era5_time_tts) 
        flag = numpy.zeros(len(Sws_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Soil water content, layer 1: 0 to 7 cm", 
                                                 group_name = "soil", units="m^3/m^3")
        pfp_utils.CreateSeries(ds_era5,"Sws",Sws_era5_tts,flag,attr)
        # LEVEL2: 
        if nDims==3:
            # 3 dimensions
            Sws_3d = era5_file.variables["swvl2"][:,:,:]
            Sws_era5_1hr = Sws_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Sws_3d = era5_file.variables["swvl2"][:,:,:,:]
            data0 = Sws_3d[:,0,site_lat_index,site_lon_index]
            data1 = Sws_3d[:,1,site_lat_index,site_lon_index]
            Sws_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Sws_3d = era5_file.variables["swvl2"][:,:,:]
        #Sws_era5_1hr = Sws_3d[:,site_lat_index,site_lon_index]
        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Sws_era5_1hr) 
        # get the coefficient at the tower time step
        Sws_era5_tts = int_fn(era5_time_tts) 
        flag = numpy.zeros(len(Sws_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Soil water content, layer 2: 7 to 28 cm", 
                                                 group_name = "soil", units="m^3/m^3")
        pfp_utils.CreateSeries(ds_era5,"Sws2",Sws_era5_tts,flag,attr)
        # LEVEL3: 
        if nDims==3:
            # 3 dimensions
            Sws_3d = era5_file.variables["swvl3"][:,:,:]
            Sws_era5_1hr = Sws_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Sws_3d = era5_file.variables["swvl3"][:,:,:,:]
            data0 = Sws_3d[:,0,site_lat_index,site_lon_index]
            data1 = Sws_3d[:,1,site_lat_index,site_lon_index]
            Sws_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Sws_3d = era5_file.variables["swvl3"][:,:,:]
        #Sws_era5_1hr = Sws_3d[:,site_lat_index,site_lon_index]
        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Sws_era5_1hr) 
        # get the coefficient at the tower time step
        Sws_era5_tts = int_fn(era5_time_tts) 
        flag = numpy.zeros(len(Sws_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Soil water content, layer 3: 28 to 100 cm", 
                                                 group_name = "soil", units="m^3/m^3")
        pfp_utils.CreateSeries(ds_era5,"Sws3",Sws_era5_tts,flag,attr)
        # LEVEL4: 
        if nDims==3:
            # 3 dimensions
            Sws_3d = era5_file.variables["swvl4"][:,:,:]
            Sws_era5_1hr = Sws_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Sws_3d = era5_file.variables["swvl4"][:,:,:,:]
            data0 = Sws_3d[:,0,site_lat_index,site_lon_index]
            data1 = Sws_3d[:,1,site_lat_index,site_lon_index]
            Sws_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Sws_3d = era5_file.variables["swvl4"][:,:,:]
        #Sws_era5_1hr = Sws_3d[:,site_lat_index,site_lon_index]
        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Sws_era5_1hr) 
        # get the coefficient at the tower time step
        Sws_era5_tts = int_fn(era5_time_tts) 
        flag = numpy.zeros(len(Sws_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Soil water content, layer 4: 100 to 289 cm", 
                                                 group_name = "soil", units="m^3/m^3")
        pfp_utils.CreateSeries(ds_era5,"Sws4",Sws_era5_tts,flag,attr)

        # === SOIL TEMPERATURE Ts === #
        # LEVEL1: 
        # Interpolate the 1 hourly soil temperature to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Ts_3d = era5_file.variables["stl1"][:,:,:]
            Ts_era5_1hr = Ts_3d[:,site_lat_index,site_lon_index] - 273.15
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Ts_3d = era5_file.variables["stl1"][:,:,:,:]
            data0 = Ts_3d[:,0,site_lat_index,site_lon_index] - 273.15
            data1 = Ts_3d[:,1,site_lat_index,site_lon_index] - 273.15
            Ts_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Ts_3d = era5_file.variables["stl1"][:,:,:]
        #Ts_era5_1hr = Ts_3d[:,site_lat_index,site_lon_index] - 273.15
        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Ts_era5_1hr) 
        # get the coefficient at the tower time step
        Ts_era5_tts = int_fn(era5_time_tts) 
        flag = numpy.zeros(len(Ts_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Soil temperature, layer 1: 0 to 7 cm", 
                                                 group_name = "soil", units="degC")
        pfp_utils.CreateSeries(ds_era5,"Ts",Ts_era5_tts,flag,attr)
        # LEVEL2: 
        # Interpolate the 1 hourly soil temperature to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Ts_3d = era5_file.variables["stl2"][:,:,:]
            Ts_era5_1hr = Ts_3d[:,site_lat_index,site_lon_index] - 273.15
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Ts_3d = era5_file.variables["stl2"][:,:,:,:]
            data0 = Ts_3d[:,0,site_lat_index,site_lon_index] - 273.15
            data1 = Ts_3d[:,1,site_lat_index,site_lon_index] - 273.15
            Ts_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Ts_3d = era5_file.variables["stl2"][:,:,:]
        #Ts_era5_1hr = Ts_3d[:,site_lat_index,site_lon_index] - 273.15
        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Ts_era5_1hr) 
        # get the coefficient at the tower time step
        Ts_era5_tts = int_fn(era5_time_tts) 
        flag = numpy.zeros(len(Ts_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Soil temperature, layer 2: 7 to 28 cm", 
                                                 group_name = "soil", units="degC")
        pfp_utils.CreateSeries(ds_era5,"Ts2",Ts_era5_tts,flag,attr)
        # LEVEL3: 
        # Interpolate the 1 hourly soil temperature to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Ts_3d = era5_file.variables["stl3"][:,:,:]
            Ts_era5_1hr = Ts_3d[:,site_lat_index,site_lon_index] - 273.15
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Ts_3d = era5_file.variables["stl3"][:,:,:,:]
            data0 = Ts_3d[:,0,site_lat_index,site_lon_index] - 273.15
            data1 = Ts_3d[:,1,site_lat_index,site_lon_index] - 273.15
            Ts_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Ts_3d = era5_file.variables["stl3"][:,:,:]
        #Ts_era5_1hr = Ts_3d[:,site_lat_index,site_lon_index] - 273.15
        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Ts_era5_1hr) 
        # get the coefficient at the tower time step
        Ts_era5_tts = int_fn(era5_time_tts) 
        flag = numpy.zeros(len(Ts_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Soil temperature, layer 3: 28 to 100 cm", 
                                                 group_name = "soil", units="degC")
        pfp_utils.CreateSeries(ds_era5,"Ts3",Ts_era5_tts,flag,attr)
        # LEVEL4: 
        # Interpolate the 1 hourly soil temperature to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        if nDims==3:
            # 3 dimensions
            Ts_3d = era5_file.variables["stl4"][:,:,:]
            Ts_era5_1hr = Ts_3d[:,site_lat_index,site_lon_index] - 273.15
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            Ts_3d = era5_file.variables["stl4"][:,:,:,:]
            data0 = Ts_3d[:,0,site_lat_index,site_lon_index] - 273.15
            data1 = Ts_3d[:,1,site_lat_index,site_lon_index] - 273.15
            Ts_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #Ts_3d = era5_file.variables["stl4"][:,:,:]
        #Ts_era5_1hr = Ts_3d[:,site_lat_index,site_lon_index] - 273.15
        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, Ts_era5_1hr) 
        # get the coefficient at the tower time step
        Ts_era5_tts = int_fn(era5_time_tts) 
        flag = numpy.zeros(len(Ts_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Soil temperature, layer 4: 100 to 289 cm", 
                                                 group_name = "soil", units="degC")
        pfp_utils.CreateSeries(ds_era5,"Ts4",Ts_era5_tts,flag,attr)

        # === U-WIND COMPONENT U === #
        # Interpolate the 1 hourly U and V components to the tower time step
        # NOTE: ERA-5 variables are dimensioned [time,latitude,longitude]
        # U first ...
        if nDims==3:
            # 3 dimensions
            U_3d = era5_file.variables["u10"][:,:,:]
            U_era5_1hr = U_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            U_3d = era5_file.variables["u10"][:,:,:,:]
            data0 = U_3d[:,0,site_lat_index,site_lon_index]
            data1 = U_3d[:,1,site_lat_index,site_lon_index]
            U_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #U_3d = era5_file.variables["u10"][:,:,:]
        #U_era5_1hr = U_3d[:,site_lat_index,site_lon_index]

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, U_era5_1hr) 
        # get the coefficient at the tower time step
        U_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, U_era5_1hr, k=1)
        # get the soil moisture at the tower time step
        #U_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(U_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Along wind velocity component", 
                                                 group_name = "meteorology", units="m/s")
        pfp_utils.CreateSeries(ds_era5,"U",U_era5_tts,flag,attr)

        # === V-WIND COMPONENT V === #
        # ... then V
        if nDims==3:
            # 3 dimensions
            V_3d = era5_file.variables["v10"][:,:,:]
            V_era5_1hr = V_3d[:,site_lat_index,site_lon_index]
        elif nDims==4:
            # 4 dimensions, this is specific to ERA5T data where 0 is ERA5 and 1 is ERA5T data
            V_3d = era5_file.variables["v10"][:,:,:,:]
            data0 = V_3d[:,0,site_lat_index,site_lon_index]
            data1 = V_3d[:,1,site_lat_index,site_lon_index]
            V_era5_1hr = numpy.ma.where(data0.mask == False,data0,data1)
        #V_3d = era5_file.variables["v10"][:,:,:]
        #V_era5_1hr = V_3d[:,site_lat_index,site_lon_index]

        # get the Akima interpolator function
        int_fn = scipy.interpolate.Akima1DInterpolator(era5_time_1hr, V_era5_1hr) 
        # get the coefficient at the tower time step
        V_era5_tts = int_fn(era5_time_tts) 

        # get the spline interpolation function
        #s = InterpolatedUnivariateSpline(era5_time_1hr, V_era5_1hr, k=1)
        # get the soil moisture at the tower time step
        #V_era5_tts = s(era5_time_tts)
        flag = numpy.zeros(len(V_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Across wind velocity component", 
                                                 group_name = "meteorology", units="m/s")
        pfp_utils.CreateSeries(ds_era5,"V",V_era5_tts,flag,attr)

        # === WIND SPEED Ws and WIND DIRECTION Wd === #
        # now get the wind speed and direction
        Ws_era5_tts = numpy.sqrt(U_era5_tts*U_era5_tts + V_era5_tts*V_era5_tts)
        flag = numpy.zeros(len(Ws_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Wind speed", 
                                                 group_name = "meteorology", units="m/s")
        pfp_utils.CreateSeries(ds_era5,"Ws",Ws_era5_tts,flag,attr)
        Wd_era5_tts = float(270) - numpy.arctan2(V_era5_tts,U_era5_tts)*float(180)/numpy.pi
        idx = numpy.where(Wd_era5_tts>360)[0]
        if len(idx)>0: Wd_era5_tts[idx] = Wd_era5_tts[idx] - float(360)
        flag = numpy.zeros(len(Wd_era5_tts),dtype=numpy.int32)
        attr = pfp_utils.MakeAttributeDictionary(long_name="Wind direction", 
                                                 group_name = "meteorology", units="degrees")
        pfp_utils.CreateSeries(ds_era5,"Wd",Wd_era5_tts,flag,attr)

        # === WRITE OUTPUT FILE FOR THE SITE
        # write the yearly file for this site
        ncfile = pfp_io.nc_open_write(out_file_path)
        pfp_io.nc_write_series(ncfile,ds_era5,ndims=1)
        # add this yearly file to the control file dictionary for this site
        cf_dict[site_name]["Files"]["In"][str(n+add_exist)] = out_file_path
        # tell the user we have finished this site
        logger.info("Finished "+site_name)
        logger.info("")
        
# now we need to loop over the contents of the concatenate control file dictionary
for site_name in site_list:
    cf_concat = cf_dict[site_name]
    #print('CE: cf_concat = ', cf_concat)
    cf_concat.filename = os.path.join(concat_control_path,site_name+"_concatenate.txt")
    cf_concat.write()
    msg = "Concatenating monthly files for "+site_name
    logger.info(msg)
    #print("filename",cf_concat)
    #nc_concatenate(cf_concat)
    #info = cf_concat
    #info["NetCDFConcatenate"] = {"OK": True}
    info = pfp_compliance.ParseConcatenateControlFile(cf_concat)
    pfp_io.NetCDFConcatenate(info)
    msg = "Finish concatenating monthly files for "+site_name
    logger.info(msg)
    
