#!/usr/bin/env python
import cdsapi
import datetime
import dateutil.parser
import sys
import pandas as pd
import numpy as np
# 
# downloading ECWMF ERA5 data monthly via Copernicus CDS (Copernicus Climate Data Store)
# 

# Contact:
# copernicus-support@ecmwf.int
# ===================
# dates to download: individual = yes (set start and end), no (calculate monthly)
individual = "yes"
# ERA5 hourly data on single levels from 2000 to present
era5_info = {}
start_date = "2019-11-01"
end_date   = "2019-11-02"
#start_date = "1999-12-31"
#end_date = "2018-12-31"
'''
Evaluate start date and end date to download data monthly. The month is specified to start 
'last day of previous month' and end 'last day of month', this makes sure no data gaps are 
in the concatenated files.
'''
if individual != "yes":
   sd = dateutil.parser.parse(start_date)
   ed = dateutil.parser.parse(end_date)
   date_rng = pd.date_range(start=sd, end=ed, freq='M')
   drlen = len(date_rng)
   sds = []
   eds = []
   if drlen > 0:
       fstdate = date_rng[0]
       lstdate = date_rng[-1]
       print('DATE: ', sd, fstdate, lstdate, ed)
       if sd < date_rng[0]:
           sds.append(sd.strftime("%Y-%m-%d"))
           eds.append(date_rng[0].strftime("%Y-%m-%d"))
       for i in np.arange(0, drlen-1):
           sds.append(date_rng[i].strftime("%Y-%m-%d"))
           eds.append(date_rng[i+1].strftime("%Y-%m-%d"))
       if ed > date_rng[-1]:
           sds.append(date_rng[-1].strftime("%Y-%m-%d"))
           eds.append(ed.strftime("%Y-%m-%d"))
   else:
       sds.append(sd.strftime("%Y-%m-%d"))
       eds.append(ed.strftime("%Y-%m-%d"))

if len(sys.argv)==1:
    print("Command line syntax is:")
    print(" python get_era5.py <country>")
    print("where <country> can be;")
    print(" Australia")
    print(" USA")
    sys.exit

if sys.argv[1].lower()=="australia":
    era5_info["area"] = "-10/110/-45/155", # North, West, South, East. Default: global 
    target_directory = "/data/cilli/OzFlux/Sites/ERA5/AUS/"
elif sys.argv[1].lower()=="usa":
    era5_info["area"] = "70/229.5/30/300"
    target_directory = "/mnt/AmeriFlux/ERA5/"
elif sys.argv[1].lower()=="brazil":  # latitude = -9.05; longitude = -40.316667; Jan 1 2013, Jan 1 2014
    era5_info["area"] = "-5/-45/-15/-35"
    target_directory = "/run/media/cilli/data190208/1_Work/1_OzFlux/Sites/ERA5/BRA/"
elif sys.argv[1].lower()=="nz":  #"Get ERAI data for New Zealand"
    erai_info["area"] = "-30/165/-50/180"
    target_directory = "/run/media/cilli/cillidata/cilli/1_Work/1_OzFlux/Sites/ERAI/NZ/"
elif sys.argv[1].lower()=="china":
    erai_info["area"] = "60/80/20/140"
    target_directory = "/home/peter/ChinaFlux/ERAI/"
elif sys.argv[1].lower()=="borneo":  #"Get ERAI data for Borneo/Indonesia"
    erai_info["area"] = "8.25/108/3.75/120"
    target_directory = "/run/media/cilli/cillidata/cilli/1_Work/1_OzFlux/Sites/ERAI/BORNEO/"
else:
    print("Unrecognised country option entered on command line")
    print("Valid country options are:")
    print(" australia")
    print(" usa")
    sys.exit()

# Specify product type ------------------
era5_info["product_type"] = "reanalysis"
# Specify format (grib/netcdf)-----------
era5_info["format"] = "netcdf"
# Specify variables ---------------------
era5_info["variable"] = [
            '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
            '2m_temperature','boundary_layer_height','soil_temperature_level_1',
            'soil_temperature_level_2','soil_temperature_level_3','soil_temperature_level_4',
            'surface_latent_heat_flux','surface_net_solar_radiation','surface_net_thermal_radiation',
            'surface_pressure','surface_sensible_heat_flux','surface_solar_radiation_downwards',
            'surface_thermal_radiation_downwards','total_precipitation','volumetric_soil_water_layer_1',
            'volumetric_soil_water_layer_2','volumetric_soil_water_layer_3','volumetric_soil_water_layer_4'
        ]
# Specify resolution------------------
era5_info["grid"] = [0.25,0.25] 
# Specify hours of the day------------
era5_info["time"] = [
            '00:00','01:00','02:00','03:00','04:00','05:00',
            '06:00','07:00','08:00','09:00','10:00','11:00',
            '12:00','13:00','14:00','15:00','16:00','17:00',
            '18:00','19:00','20:00','21:00','22:00','23:00'
        ]

# Call the CDS server --
server = cdsapi.Client()

#print " Processing year, month, day: ",str(sds[counter]),"to",str(eds[counter])
#era5_info["date"] = str(sds[counter]) + "/" + str(eds[counter])
if individual == "yes":
    era5_info["date"] = start_date + "/" + end_date
    print(era5_info["date"])
    era5_info["target"] = target_directory+"ERA5_"+start_date+"_to_"+end_date+".nc"
    print(era5_info["target"])
    server.retrieve(
        'reanalysis-era5-single-levels',
        {
        'product_type' : era5_info["product_type"],
        'format'       : era5_info["format"],
        'variable'     : era5_info["variable"],
        'grid'         : era5_info["grid"],
        'time'         : era5_info["time"],
        # === use "date" for multiple days
        'date'         : era5_info["date"],
        # use year, month, day for single day
        #'year': '2019',
        #'month': '09',
        #'day': '30',
        'area'         : era5_info["area"],
        },
        era5_info["target"])
else:
    for counter, value in enumerate(sds):
        era5_info["date"] = sds[counter] + "/" + eds[counter]
        print(era5_info["date"])
        era5_info["target"] = target_directory+"ERA5_"+str(sds[counter])+"_to_"+str(eds[counter]).zfill(2)+".nc"
        print(era5_info["target"])
        server.retrieve(
            'reanalysis-era5-single-levels',
            {
            'product_type' : era5_info["product_type"],
            'format'       : era5_info["format"],
            'variable'     : era5_info["variable"],
            'grid'         : era5_info["grid"],
            'time'         : era5_info["time"],
            # === use "date" for multiple days
            'date'         : era5_info["date"],
            # use year, month, day for single day
            #'year': '2019',
            #'month': '09',
            #'day': '30',
            'area'         : era5_info["area"],
            },
            era5_info["target"])
    
