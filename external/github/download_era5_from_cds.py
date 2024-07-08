#!/usr/bin/env python3
"""
  Purpose:
   Download ECWMF ERA5 data fixed/monthly via Copernicus CDS (Copernicus Climate Data Store)
   Contact: copernicus-support@ecmwf.int; ERA5 hourly data on single levels from 2000 to present
  Usage:
   python3 download_era5_from_cds.py arg1 arg2 arg3 arg4
   e.g.
   case 1: start and end date provided
   python3 download_era5_from_cds.py australia ../../Sites/ERA5/AUS/ 2020-04-30 2020-05-31
   case 2: start and end date not provided, calculated as end of last month, end of current month
   where;
      arg1 = country, e.g. "australia"
      arg2 = target_directory
      arg3 = start_date, e.g. "2020-04-30"
      arg4 = end_date, e.g. "2020-05-31"
  Side effects:
    (1) needs an account at CDS and a key, where url and key are stored in .cdsapirc
        url: https://cds.climate.copernicus.eu/api/v2
        key: 0123456789012345678901234567890123456789
    (2) for the area of Australia the data limit restricts requests to about one month of data
  Author: CME
  Date: May 2020
"""
import cdsapi
import datetime
import dateutil.parser
import sys
import pandas as pd
import numpy as np
era5_info = {}
if len(sys.argv)<3:
    print("Command line syntax is:")
    print(" python3 download_era5_from_cds.py <country> <target_directory> <start_date> <end_date>")
    print("where <start/end_date> is YYYY-MM-DD or")
    print(" python3 download_era5_from_cds.py <country> <target_directory>; start/end date calculated")
    sys.exit
elif len(sys.argv)==3:
    # evaluate last month start (last day month before and end (last day last month) days
    today = datetime.date.today()
    fstofm = today.replace(day=1)
    endlm = fstofm - datetime.timedelta(days=1)
    fstoflm = endlm.replace(day=1)
    endmb = fstoflm - datetime.timedelta(days=1)
    era5_info["start_date"] = endmb.strftime("%Y-%m-%d")
    era5_info["end_date"]   = endlm.strftime("%Y-%m-%d")
elif len(sys.argv)==5:
    era5_info["start_date"] = sys.argv[3]
    era5_info["end_date"]   = sys.argv[4]
else:
    print("Command line syntax is:")
    print(" python3 download_era5_from_cds.py <country> <target_directory> <start_date> <end_date>")
    print("where <start/end_date> is YYYY-MM-DD or")
    print(" python3 download_era5_from_cds.py <country> <target_directory>; start/end date calculated")
    sys.exit

# === Which area, see index for country given
list_country = ["australia", "indonesia", "brazil", "nz","mara"]
list_area    = ["-10/110/-45/155", "-2/114/-3/115", "-5/-45/-15/-35", "-30/165/-50/180", "-39.5/176.5/-39.75/176.75"]

era5_info["country"] = sys.argv[1]          #"australia"
era5_info["target_directory"] = sys.argv[2] #"/rdsi/market/erai_processing/ERA5/AUS/"

index = list_country.index(era5_info["country"])
era5_info["area"] = list_area[index]        #"-10/110/-45/155"
era5_info["date"] = era5_info["start_date"] + "/" + era5_info["end_date"]
era5_info["target"] = era5_info["target_directory"]+"ERA5_"+era5_info["start_date"]+"_to_"+era5_info["end_date"]+".nc"
print(era5_info["start_date"],era5_info["end_date"],era5_info["target"])

# "australia": "-10/110/-45/155"   "/rdsi/market/erai_processing/ERA5/AUS/"
# "usa":       "70/229.5/30/300"   "/mnt/AmeriFlux/ERA5/"
# "brazil":    "-5/-45/-15/-35"    "/run/media/cilli/data190208/1_Work/1_OzFlux/Sites/ERA5/BRA/ # lat = -9.05; long = -40.316667; Jan 1 2013, Jan 1 2014
# "nz":        "-30/165/-50/180"   "/run/media/cilli/cillidata/cilli/1_Work/1_OzFlux/Sites/ERAI/NZ/"
# "maraekakaho":                   "/data/cilli/OzFlux/Sites/ERA5/NZ/Maraekakaho/
# "china":     "60/80/20/140"      "/home/peter/ChinaFlux/ERAI/"
# "borneo":    "8.25/108/3.75/120" "/run/media/cilli/cillidata/cilli/1_Work/1_OzFlux/Sites/ERAI/BORNEO/"
# "indonesia": "-2/114/-3/115" target_directory = "/run/media/cilli/data190208/1_Work/1_OzFlux/Sites/ERA5/INDONESIA/"

# Specify product type
era5_info["product_type"] = "reanalysis"
# Specify format (grib/netcdf)
era5_info["format"] = "netcdf"
# Specify variables
era5_info["variable"] = [
            '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
            '2m_temperature','boundary_layer_height','soil_temperature_level_1',
            'soil_temperature_level_2','soil_temperature_level_3','soil_temperature_level_4',
            'surface_latent_heat_flux','surface_net_solar_radiation','surface_net_thermal_radiation',
            'surface_pressure','surface_sensible_heat_flux','surface_solar_radiation_downwards',
            'surface_thermal_radiation_downwards','total_precipitation','volumetric_soil_water_layer_1',
            'volumetric_soil_water_layer_2','volumetric_soil_water_layer_3','volumetric_soil_water_layer_4'
                        ]
# Specify resolution
era5_info["grid"] = [0.25,0.25] 
# Specify hours of the day
era5_info["time"] = [
            '00:00','01:00','02:00','03:00','04:00','05:00',
            '06:00','07:00','08:00','09:00','10:00','11:00',
            '12:00','13:00','14:00','15:00','16:00','17:00',
            '18:00','19:00','20:00','21:00','22:00','23:00'
                    ]
#print(era5_info)
# Call the CDS server
server = cdsapi.Client()
# Retrieve data
server.retrieve(
    'reanalysis-era5-single-levels',
    {
    'product_type' : era5_info["product_type"],
    'format'       : era5_info["format"],
    'variable'     : era5_info["variable"],
    'grid'         : era5_info["grid"],
    'time'         : era5_info["time"],
    'date'         : era5_info["date"],
    'area'         : era5_info["area"]
    },
    era5_info["target"])
