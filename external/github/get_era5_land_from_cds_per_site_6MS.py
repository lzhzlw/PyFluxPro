#!/usr/bin/env python
# ==================================================================================
# Program: get_era5_land_from_cds_per_site.py ../controlfiles/ERA5/get_era5_land.txt
# CME: Sep 2019
# ==================================================================================
import cdsapi
# standard
from collections import OrderedDict
import datetime
import dateutil.parser
import glob
import os
import sys
import time
from decimal import *
# 3rd party
from configobj import ConfigObj
import netCDF4
import numpy
import pandas
import pytz
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline
import xlrd
# check the scripts directory is present
if not os.path.exists("../scripts/"):
    print("era5land2nc: the scripts directory is missing")
    sys.exit()
# since the scripts directory is there, try importing the modules
sys.path.append('../scripts')
# PFP
import pfp_io
import pfp_log
import pfp_utils
"""
This script reads hourly ERA5-land data from the CDS server and picks out the tower site data 
as a timeline.
"""
now = datetime.datetime.now()
log_file_name = 'era5land2_get_nc_' + now.strftime("%Y%m%d%H%M") + ".log"
logger = pfp_log.init_logger("pfp_log", log_file_name, to_file=True, to_screen=True)

def read_site_master(xl_file_path, sheet_name):
    """
    Reads the site master file with entries of site name and site location (latitude, longitude).
    site_master.xls can have more than one sheet. In era5land2nc.txt the sheetname is specified.
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
xl_file_path         = cf["Files"]["xl_file_path"]
xl_sheet_name        = cf["Files"]["xl_sheet_name"]

target_directory     = cf["Files"]["target_directory"]
start_date           = cf["Files"]["start_date"]
end_date             = cf["Files"]["end_date"]
era5_land_resolution = float(cf["Files"]["era5_land_resolution"])

# get the site information from the site master spreadsheet
site_info = read_site_master(xl_file_path, xl_sheet_name)
# get a list of sites
site_list = list(site_info.keys())

for site_name in site_list:
    # calculate the era5 area to be downloaded, use the 4 grid points around the tower site
    # grid resolution is 0.1 deg, hence 
    era5_info = {}
    
    ## get the metadata from the site master file information
    site_latitude  = site_info[site_name]["Latitude"]
    site_longitude = site_info[site_name]["Longitude"]
    if era5_land_resolution > 0.0 and era5_land_resolution < 1.0:
        print(site_latitude,site_longitude)
        nearlat = Decimal(site_latitude).quantize(Decimal('0.1'))
        nearlon = Decimal(site_longitude).quantize(Decimal('0.1'))
        print(nearlat,nearlon)
        north = Decimal(site_latitude).quantize(Decimal('0.1'), rounding=ROUND_UP)
        west  = Decimal(site_longitude).quantize(Decimal('0.1'), rounding=ROUND_DOWN)
        south = Decimal(site_latitude).quantize(Decimal('0.1'), rounding=ROUND_DOWN)
        east  = Decimal(site_longitude).quantize(Decimal('0.1'), rounding=ROUND_UP)
        print(east, west, north, south)
        era5_info["area"] = str(nearlat)+"/"+str(nearlon)+"/"+str(nearlat)+"/"+str(nearlon) # # North, West, South, East. Default: global 
        #era5_info["area"] = str(north)+"/"+str(west)+"/"+str(south)+"/"+str(east) # # North, West, South, East. Default: global 
    # Loxton = "-34.4/140.6/-34.5/140.7"
    print(era5_info["area"])
    # ======================================================================================
    # downloading ECWMF ERA5 data monthly via Copernicus CDS (Copernicus Climate Data Store)
    # ======================================================================================
    # Contact:
    # copernicus-support@ecmwf.int
    # ERA5-land hourly data on single levels from 2001 to present
    #start_date = "2019-05-31"
    #end_date = "2019-06-30"
    '''
    Evaluate start date and end date to download data monthly. The month is specified to start 
    'last day of previous month' and end 'last day of month', this makes sure no data gaps are 
    in the concatenated files.
    '''
    sd = dateutil.parser.parse(start_date)
    ed = dateutil.parser.parse(end_date)
    #date_rng = pandas.date_range(start=sd, end=ed, freq='M')
    # === end date range specification ===
    date_rng = pandas.date_range(start=sd,end=ed,freq='6MS')
    drlen = len(date_rng)
    fstdate = date_rng[0]
    lstdate = date_rng[-1]
    print('DATE: ', sd, fstdate, lstdate, ed)
    sds = []
    eds = []
    if sd < date_rng[0]:
        sds.append(sd.strftime("%Y-%m-%d"))
        eds.append(date_rng[0].strftime("%Y-%m-%d"))
    for i in numpy.arange(0, drlen-1):
        sds.append(date_rng[i].strftime("%Y-%m-%d"))
        eds.append(date_rng[i+1].strftime("%Y-%m-%d"))
    if ed > date_rng[-1]:
        sds.append(date_rng[-1].strftime("%Y-%m-%d"))
        eds.append(ed.strftime("%Y-%m-%d"))
    # === end date range specification ===
    #sds = sd.strftime("%Y-%m-%d")
    #eds = ed.strftime("%Y-%m-%d")
    # Specify product type ------------------
    era5_info["product_type"] = "reanalysis"
    # Specify format (grib/netcdf)-----------
    era5_info["format"] = "netcdf"
    # Specify variables ---------------------
    era5_info["variable"] = [
                '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
                '2m_temperature','evapotranspiration','soil_temperature_level_1',
                'soil_temperature_level_2','soil_temperature_level_3','soil_temperature_level_4',
                'surface_latent_heat_flux','surface_net_solar_radiation','surface_net_thermal_radiation',
                'surface_pressure','surface_sensible_heat_flux','surface_solar_radiation_downwards',
                'surface_thermal_radiation_downwards','total_precipitation','volumetric_soil_water_layer_1',
                'volumetric_soil_water_layer_2','volumetric_soil_water_layer_3','volumetric_soil_water_layer_4'
            ]
    # Specify resolution------------------
    era5_info["grid"] = [era5_land_resolution,era5_land_resolution] 
    # Specify hours of the day------------
    era5_info["time"] = [
                '00:00','01:00','02:00','03:00','04:00','05:00',
                '06:00','07:00','08:00','09:00','10:00','11:00',
                '12:00','13:00','14:00','15:00','16:00','17:00',
                '18:00','19:00','20:00','21:00','22:00','23:00'
            ]
    # Call the CDS server --
    server = cdsapi.Client()
    for counter, value in enumerate(sds):

        era5_info["date"] = sds[counter] + "/" + eds[counter]
        print(era5_info["date"])
        era5_info["target"] = target_directory+site_name+"ERA5land_"+str(sds[counter])+"_to_"+str(eds[counter]).zfill(2)+".nc"
        print(era5_info["target"])
        
        server.retrieve(
            'reanalysis-era5-land',
            {
            'format'       : era5_info["format"],
            'variable'     : era5_info["variable"],
            'grid'         : era5_info["grid"],
            'time'         : era5_info["time"],
            'date'         : era5_info["date"],
            'area'         : era5_info["area"]
            },
            era5_info["target"])
        
