"""
This script reads hourly ERA5 data and picks out the tower site data as a timeline
it then converts these to tower time step using an Akima univariate spline interpolation.
Akima univariante spline interpolation insures the spline is going through the data point.
"""
# standard modules
import copy
import datetime
import faulthandler
import logging
import os
import shutil
import sys
import _pickle as cPickle
# 3rd party
from configobj import ConfigObj
import cftime
import dateutil
import netCDF4
import numpy
import xlrd
# PFP modules
# cwd = os.getcwd()
cwd = os.path.dirname(__file__)
sys.path.insert(0, cwd[:cwd.index("external_mysite")-1])
import scripts.constants as c
import scripts.meteorologicalfunctions as pfp_mf
import scripts.pfp_io as pfp_io
import scripts.pfp_log as pfp_log
import scripts.pfp_ts as pfp_ts
import scripts.pfp_utils as pfp_utils

faulthandler.enable()

now = datetime.datetime.now()
log_file_name = "era52nc_" + now.strftime("%Y%m%d%H%M") + ".log"
log_file_path = "logfiles"
if not os.path.isdir(log_file_path):
    os.mkdir(log_file_path)
log_file_uri = os.path.join(log_file_path, log_file_name)
logger = pfp_log.CreateLogger("pfp_log", log_file_name=log_file_uri, to_screen=True)

def read_site_master(site_master_file_path, sheet_name):
    """
    """
    xl_book = xlrd.open_workbook(site_master_file_path)
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
    site_info = {}
    for n in range(first_data_row,last_row):
        site_name = xl_sheet.cell(n,0).value
        site_name = site_name.replace(" ","")
        site_info[site_name] = {}
        for item in header_row_values[1:]:
            i = header_row_values.index(item)
            site_info[site_name][item] = xl_sheet.cell(n,i).value
    return site_info
def build_concatenation_dictionary(site_info, cfg):
    existing_era5_base_path = cfg["Files"]["existing_era5_base_path"]
    concatenation_info = {}
    sites = list(site_info.keys())
    for site in sites:
        era5_file_path = os.path.join(existing_era5_base_path, site, "Data", "ERA5")
        era5_file_name = site + "_ERA5" + ".nc"
        out_file_name = os.path.join(era5_file_path, era5_file_name)
        #if cfg["Files"].as_bool("from_cleaned"):
        #    era5_file_path = os.path.join(era5_file_path, "cleaned")
        in_file_name = os.path.join(era5_file_path, era5_file_name)
        in_file_names = [in_file_name]
        concatenation_info[site] = {"NetCDFConcatenate": {}, "RemoveIntermediateSeries": {}}
        cisn = concatenation_info[site]["NetCDFConcatenate"]
        cisn["OK"] = True
        cisn["in_file_names"] = in_file_names
        cisn["out_file_name"] = out_file_name
        cisn["NumberOfDimensions"] = 3
        cisn["MaxGapInterpolate"] = 0
        cisn["FixTimeStepMethod"] = "round"
        cisn["Truncate"] = "No"
        cisn["TruncateThreshold"] = 50
        cisn["SeriesToCheck"] = ["all"]
        cisn["time_coverage_start"] = []
        cisn["time_coverage_end"] = []
        cisn["chrono_files"] = []
        cisn["labels"] = []
        cisn["attributes"] = ['long_name', 'standard_name', 'statistic_type', 'units']
        cisr = concatenation_info[site]["RemoveIntermediateSeries"]
        cisr["KeepIntermediateSeries"] = "No"
        cisr["not_output"] = []
    return concatenation_info
def convert_utc_to_local_standard(utc_notz, site_timezone):
    # get the UTC timezone object
    from_zone = dateutil.tz.gettz('UTC')
    # get the site timezone object
    to_zone = dateutil.tz.gettz(site_timezone)
    # add the UTC timezone information to the ERA5 datetime
    utc_tz = [dt.replace(tzinfo=from_zone) for dt in utc_notz]
    # get the local daylight savings time at the site from UTC
    loc_dst = [dt.astimezone(to_zone) for dt in utc_tz]
    # convert from daylight savings time to standard time at the site
    loc_ast = [dt-dt.dst() for dt in loc_dst]
    # remove the timezone information from the local standard time at the site
    loc_ast_notz = [dt.replace(tzinfo=None) for dt in loc_ast]
    return loc_ast_notz
def delete_intermediate_variables(ds, intermediate_variables):
    labels = list(ds.root["Variables"].keys())
    for iv_label in intermediate_variables:
        label = iv_label
        if label in labels:
            ds.root["Variables"].pop(label)
    return
def init_data(cfg, site_info):
    labels = list(cfg["Variables"].keys())
    sites = list(site_info.keys())
    # initialise the data dictionary
    data = {}
    for site in sites:
        data[site] = {}
        data[site]["dt_utc"] = []
        data[site]["DateTime"] = {"Data": [],
                                  "Flag": [],
                                  "Attr": {"long_name": "Datetime in local timezone",
                                           "units": ""}}
        for label in labels:
            attr = copy.deepcopy(cfg["Variables"][label]["Attr"])
                    # we use lists because they are fast
            data[site][label] = {"Data": [], "Flag": [], "Attr": attr}
    return data
# def calculate_upwelling_shortwave_radiation(ds):
#     nrecs = int(ds.root["Attributes"]["nc_nrecs"])
#     zeros = numpy.zeros(nrecs)
#     ones = numpy.ones(nrecs)
#     Fsd = pfp_utils.GetVariable(ds, "Fsd")
#     Fnsw = pfp_utils.GetVariable(ds, "Fnsw")
#     Fsu = pfp_utils.CreateEmptyVariable("Fsu", nrecs,
#                                          datetime=Fsd["DateTime"])
#     Fsu["Data"] = Fsd["Data"] - Fnsw["Data"]
#     Fsu["Flag"] = numpy.where(numpy.ma.getmaskarray(Fsd["Data"]) == True, ones, zeros)
#     Fsu["Attr"] = {"standard_name": "surface_upwelling_shortwave_flux_in_air",
#                     "long_name": "Up-welling shortwave radiation",
#                     "units": "W/m^2", "statistic_type": "average"}
#     pfp_utils.CreateVariable(ds, Fsu)
#     return
# def calculate_upwelling_longwave_radiation(ds):
#     nrecs = int(ds.root["Attributes"]["nc_nrecs"])
#     zeros = numpy.zeros(nrecs)
#     ones = numpy.ones(nrecs)
#     Fld = pfp_utils.GetVariable(ds, "Fld")
#     Fnlw = pfp_utils.GetVariable(ds, "Fnlw")
#     Flu = pfp_utils.CreateEmptyVariable("Flu", nrecs,
#                                          datetime=Fld["DateTime"])
#     Flu["Data"] = Fld["Data"] - Fnlw["Data"]
#     Flu["Flag"] = numpy.where(numpy.ma.getmaskarray(Fld["Data"]) == True, ones, zeros)
#     Flu["Attr"] = {"standard_name": "surface_upwelling_longwave_flux_in_air",
#                     "long_name": "Up-welling shortwave radiation",
#                     "units": "W/m^2", "statistic_type": "average"}
#     pfp_utils.CreateVariable(ds, Flu)
#     return
def calculate_net_radiation(ds):
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    Fnsw = pfp_utils.GetVariable(ds, "Fnsw")
    Fnlw = pfp_utils.GetVariable(ds, "Fnlw")
    Fn = pfp_utils.CreateEmptyVariable("Fn", nrecs,
                                        datetime=Fnsw["DateTime"])
    Fn["Data"] = Fnsw["Data"] + Fnlw["Data"]
    Fn["Flag"] = numpy.where(numpy.ma.getmaskarray(Fn["Data"]) == True, ones, zeros)
    Fn["Attr"] = {"standard_name": "surface_net_downwawrd_radiative_flux",
                   "long_name": "Net radiation",
                   "units": "W/m^2", "statistic_type": "average"}
    pfp_utils.CreateVariable(ds, Fn)
    return
# def calculate_available_energy(ds):
#     nrecs = int(ds.root["Attributes"]["nc_nrecs"])
#     zeros = numpy.zeros(nrecs)
#     ones = numpy.ones(nrecs)
#     Fh = pfp_utils.GetVariable(ds, "Fh")
#     Fe = pfp_utils.GetVariable(ds, "Fe")
#     Fa = pfp_utils.CreateEmptyVariable("Fa", nrecs,
#                                         datetime=Fh["DateTime"])
#     Fa["Data"] = Fh["Data"] + Fe["Data"]
#     Fa["Flag"] = numpy.where(numpy.ma.getmaskarray(Fa["Data"]) == True, ones, zeros)
#     Fa["Attr"] = {"long_name": "Available energy",
#                    "units": "W/m^2", "statistic_type": "average"}
#     pfp_utils.CreateVariable(ds, Fa)
#     return
# def calculate_ground_heat_flux(ds):
#     # as residual from net rad - avail energy
#     nrecs = int(ds.root["Attributes"]["nc_nrecs"])
#     zeros = numpy.zeros(nrecs)
#     ones = numpy.ones(nrecs)
#     Fn = pfp_utils.GetVariable(ds, "Fn")
#     Fa = pfp_utils.GetVariable(ds, "Fa")
#     Fg = pfp_utils.CreateEmptyVariable("Fg", nrecs,
#                                         datetime=Fn["DateTime"])
#     Fg["Data"] = Fn["Data"] - Fa["Data"]
#     Fg["Flag"] = numpy.where(numpy.ma.getmaskarray(Fg["Data"]) == True, ones, zeros)
#     Fg["Attr"] = {"standard_name": "downward_heat_flux_in_soil",
#                    "long_name": "Ground heat flux",
#                    "units": "W/m^2", "statistic_type": "average"}
#     pfp_utils.CreateVariable(ds, Fg)
#     return
def calculate_relative_humidity(ds):
    # from dew point temperature
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    Td = pfp_utils.GetVariable(ds, "Td")
    Ta = pfp_utils.GetVariable(ds, "Ta")
    RH = pfp_utils.CreateEmptyVariable("RH", nrecs,
                                        datetime=Td["DateTime"])
    RH["Data"] = pfp_mf.relativehumidityfromdewpoint(Td["Data"], Ta["Data"])
    RH["Flag"] = numpy.where(numpy.ma.getmaskarray(RH["Data"]) == True, ones, zeros)
    RH["Attr"] = {"standard_name": "relative_humidity",
                   "long_name": "Relative humidity",
                   "units": "percent", "statistic_type": "average"}
    pfp_utils.CreateVariable(ds, RH)
    return
# def calculate_absolute_humidity(ds):
#     # from relative humidity
#     nrecs = int(ds.root["Attributes"]["nc_nrecs"])
#     zeros = numpy.zeros(nrecs)
#     ones = numpy.ones(nrecs)
#     RH = pfp_utils.GetVariable(ds, "RH")
#     Ta = pfp_utils.GetVariable(ds, "Ta")
#     AH = pfp_utils.CreateEmptyVariable("AH", nrecs,
#                                         datetime=RH["DateTime"])
#     AH["Data"] = pfp_mf.absolutehumidityfromrelativehumidity(Ta["Data"], RH["Data"])
#     AH["Flag"] = numpy.where(numpy.ma.getmaskarray(AH["Data"]) == True, ones, zeros)
#     AH["Attr"] = {"standard_name": "mass_concentration_of_water_vapor_in_air",
#                    "long_name": "Absolute humidity",
#                    "units": "g/m^3", "statistic_type": "average"}
#     pfp_utils.CreateVariable(ds, AH)
#     return
# def calculate_specific_humidity(ds):
#     # from relative humidity
#     nrecs = int(ds.root["Attributes"]["nc_nrecs"])
#     zeros = numpy.zeros(nrecs)
#     ones = numpy.ones(nrecs)
#     Ta = pfp_utils.GetVariable(ds,"Ta")
#     ps = pfp_utils.GetVariable(ds,"ps")
#     RH = pfp_utils.GetVariable(ds,"RH")
#     SH = pfp_utils.CreateEmptyVariable("SH", nrecs,
#                                         datetime=RH["DateTime"])
#     SH["Data"] = pfp_mf.specifichumidityfromrelativehumidity(RH["Data"], Ta["Data"], ps["Data"])
#     SH["Flag"] = numpy.where(numpy.ma.getmaskarray(SH["Data"]) == True, ones, zeros)
#     SH["Attr"] = {"standard_name": "specific_humidity",
#                    "long_name": "Specific humidity",
#                    "units": "kg/kg", "statistic_type": "average"}
#     pfp_utils.CreateVariable(ds, SH)
def calculate_vapour_pressure_deficit(ds):
    # from dew point temperature
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    Ta = pfp_utils.GetVariable(ds, "Ta")
    RH = pfp_utils.GetVariable(ds, "RH")
    VPD = pfp_utils.CreateEmptyVariable("VPD", nrecs,
                                        datetime=Ta["DateTime"])
    VPD["Data"] = pfp_mf.vpd_from_ta_rh(Ta["Data"], RH["Data"])
    VPD["Flag"] = numpy.where(numpy.ma.getmaskarray(VPD["Data"]) == True, ones, zeros)
    VPD["Attr"] = {"standard_name": "water_vapor_saturation_deficit_in_air",
                   "long_name": "Vapour pressure deficit",
                   "units": "hPa", "statistic_type": "average"}
    pfp_utils.CreateVariable(ds, VPD)
    return
def calculate_ws_and_wd_from_u_and_v(ds):
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    U = pfp_utils.GetVariable(ds, "U")
    V = pfp_utils.GetVariable(ds, "V")
    Ws = pfp_utils.CreateEmptyVariable("Ws", nrecs,
                                        datetime=U["DateTime"])
    Wd = pfp_utils.CreateEmptyVariable("Wd", nrecs,
                                        datetime=U["DateTime"])
    # get the wind speed and direction from the components
    Wd["Data"] = float(270) - (numpy.degrees(numpy.ma.arctan2(V["Data"], U["Data"])))
    Wd["Data"] = numpy.ma.mod(Wd["Data"], 360)
    Ws["Data"] = numpy.ma.sqrt(U["Data"]*U["Data"] + V["Data"]*V["Data"])
    # mask wind direction when the wind speed is less than 0.01
    Wd["Data"] = numpy.ma.masked_where(Ws["Data"] < 0.01, Wd["Data"])
    # now set the QC flag
    Ws["Flag"] = numpy.where(numpy.ma.getmaskarray(Ws["Data"]) == True, ones, zeros)
    Wd["Flag"] = numpy.where(numpy.ma.getmaskarray(Wd["Data"]) == True, ones, zeros)
    # update the variable attributes
    Ws["Attr"] = {"standard_name": "wind_speed", "long_name": "Wind speed",
                  "units": "m/s", "statistic_type": "average"}
    Wd["Attr"] = {"standard_name": "wind_from_direction", "long_name": "Wind direction",
                  "units": "degrees", "statistic_type": "average"}
    pfp_utils.CreateVariable(ds, Ws)
    pfp_utils.CreateVariable(ds, Wd)
    return
# ======================================= MAIN ========================================
logger = logging.getLogger("pfp_log")
# read the control file
cfg_file_path = "process_era52nc.txt"
msg = " Loading the control file"
logger.info(msg)
cfg = pfp_io.get_controlfilecontents(cfg_file_path)
cfg_labels = [l for l in list(cfg["Variables"].keys()) if "nc" in list(cfg["Variables"][l].keys())]
# read the site master workbook
site_info = read_site_master(cfg["Files"]["site_master_file_path"], cfg["Files"]["xl_sheet_name"])
sites = list(site_info.keys())
# build the information dictionary for concatenation
concatenation_info = build_concatenation_dictionary(site_info, cfg)

dt_utc = []
data = init_data(cfg, site_info)
new_era5_base_path = cfg["Files"]["new_era5_base_path"]
new_era5_files = sorted([f.path for f in os.scandir(new_era5_base_path) if f.is_file()])
n_new_era5_files = len(new_era5_files)
era5_timestep = 60
if n_new_era5_files == 0:
    msg = " No files found in " + new_era5_base_path
    logger.info(msg)
    sys.exit()

if cfg["Files"].as_bool("from_file"):
    # loop over the ERA5 files and read the ERA5 files,
    # and put the data into a dictionary ordered by site
    msg = " Processing ERA5 files from " + new_era5_files[0].split(os.sep)[-1]
    msg += " to " + new_era5_files[-1].split(os.sep)[-1]
    logger.info(msg)
    warnings = {"all": [], "cant_read": [], "netcdf_error": []}
    for new_era5_file in new_era5_files:
        #new_era5_file_path = new_era5_file
        try:
            era5_file = netCDF4.Dataset(new_era5_file)
        except:
            base_name = os.path.basename(new_era5_file)
            msg = " Can't read " + base_name + ", skipping file ..."
            #logger.warning(msg)
            warnings["cant_read"].append(msg)
            #era5_file.close()
            continue
        # get the Python datetime (UTC)
        nc_time_data = numpy.ma.filled(era5_file.variables["time"])
        nc_time_units = getattr(era5_file.variables["time"], "units")
        ntime = len(nc_time_data) #era5_file.variables["time"])
        dt_utc = cftime.num2pydate(nc_time_data, nc_time_units)
        # get the latitude and longitude
        latitude = era5_file.variables["latitude"]
        longitude = era5_file.variables["longitude"]
        # get the latitude and longitude resolution
        delta_latitude = (max(latitude) - min(latitude))/(len(latitude)-1)
        delta_longitude = (max(longitude) - min(longitude))/(len(longitude)-1)
        # loop over sites
        for site in sites:
            # get the site latitude and longitude
            site_latitude = site_info[site]["Latitude"]
            site_longitude = site_info[site]["Longitude"]
            #
            # we need a check here to make sure the selected grid square is land, not water
            #
            y = int(((latitude[0]-site_latitude)/delta_latitude)+0.5)
            x = int(((site_longitude-longitude[0])/delta_longitude)+0.5)
            # loop over variables to read from the ERA5 file
            data[site]["dt_utc"].append(dt_utc[:])
            n_cfg_labels = len(cfg_labels)
            for label in cfg_labels:
                d = numpy.full([ntime], -9999, dtype=numpy.float32)
                era5_label = cfg["Variables"][label]["nc"]["name"]
                era5_variable = era5_file.variables[era5_label]
                # check the variable dimensions
                if len(era5_variable.shape) == 3:
                    # ERA5 variables have 3 dimensions if all final or all temp data
                    try:
                        d = era5_variable[:, y, x].data
                    except:
                        base_name = os.path.basename(new_era5_file)
                        msg = " netCDF error: " + base_name + " (" + site + ", "
                        msg += label + ")"
                        #logger.warning(msg)
                        warnings["netcdf_error"].append(msg)
                        # create dummy variable values
                        d = numpy.full((ntime), float(c.missing_value), dtype=numpy.float32)
                elif len(era5_variable.shape) == 4:
                    # but if mix of final and temp data then it is 4-D ...
                    try:
                        d0 = era5_variable[:, 0, y, x]#.data
                        d1 = era5_variable[:, 1, y, x]#.data
                        d = numpy.ma.where(d0.mask == False,d0,d1).data
                    except:
                        base_name = os.path.basename(new_era5_file)
                        msg = " netCDF error: " + base_name + " (" + site + ", "
                        msg += label + ")"
                        #logger.warning(msg)
                        warnings["netcdf_error"].append(msg)
                        # create dummy variable values
                        d = numpy.full((ntime), float(c.missing_value), dtype=numpy.float32)
                else:
                    msg = " Unexpected dimensions in netCDF file"
                    logger.error(msg)
                    raise Exception
                # loop over the 1x1 cut out around the site
                data[site][label]["Data"].append(d[:])
        era5_file.close()
    # move the processed files from "download" to "processed"
    if cfg["Options"].as_bool("move_to_processed"):
        msg = " Moving " + new_era5_file.split(os.sep)[-1] + " to processed"
        logger.info(msg)
    #processed_era5_directory = new_era5_directory.replace("downloaded", "processed")
    #shutil.move(new_era5_directory, processed_era5_directory)
    wlist = [warnings[k] for k in warnings.keys()]
    warnings["all"] = [item for sublist in wlist for item in sublist]
    if len(warnings["all"]) != 0:
        n_warnings = len(list(set(warnings["all"])))
        msg = str(n_warnings) + " warnings occurred"
        logger.warning("")
        logger.warning(msg)
        logger.warning("")
    # pickle the data in case we want to use it again
    with open(r"era5.data.pickle", "wb") as output_file:
        cPickle.dump(data, output_file)
else:
    # read the pickled data
    msg = " Reading from pickle file"
    logger.info(msg)
    import _pickle as cPickle
    with open(r"era5.data.pickle", "rb") as input_file:
        data = cPickle.load(input_file)

# dictionary to hold data structures for each site
dss_ats = {}
msg = " Creating data structure for each site"
logger.info(msg)
for site in sites:
    # create a data structure for this site
    dss_ats[site] = pfp_io.DataStructure()
    # convert UTC datetime to local standard time
    dt_utc = numpy.hstack(data[site]["dt_utc"])
    dt_loc = numpy.array(convert_utc_to_local_standard(dt_utc, site_info[site]["Time zone"]))
    nrecs = len(dt_loc)
    # add the global attributes
    dss_ats[site].root["Attributes"]["site_name"] = site.replace(" ", "")
    dss_ats[site].root["Attributes"]["time_zone"] = site_info[site]["Time zone"]
    dss_ats[site].root["Attributes"]["latitude"] = site_info[site]["Latitude"]
    dss_ats[site].root["Attributes"]["longitude"] = site_info[site]["Longitude"]
    dss_ats[site].root["Attributes"]["time_coverage_start"] = str(dt_loc[0])
    dss_ats[site].root["Attributes"]["time_coverage_end"] = str(dt_loc[-1])
    dss_ats[site].root["Attributes"]["processing_level"] = "L1"
    # 
    dss_ats[site].root["Attributes"]["nc_nrecs"] = nrecs
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    # put the Python datetime in the data structure
    attr = {"long_name": "Datetime in local timezone", "units": ""}
    variable = {"Label": "DateTime", "Data": dt_loc, "Flag": zeros, "Attr": attr}
    pfp_utils.CreateVariable(dss_ats[site], variable)
    # get the time step
    ts = numpy.mean(pfp_utils.get_timestep(dss_ats[site])//60)
    ts = pfp_utils.roundtobase(ts, base=30)
    dss_ats[site].root["Attributes"]["time_step"] = str(int(ts))
    # put the UTC datetime in the data structure
    attr = {"long_name": "Datetime in UTC", "units": ""}
    variable = {"Label": "DateTime_UTC", "Data": dt_utc, "Flag": zeros, "Attr": attr}
    pfp_utils.CreateVariable(dss_ats[site], variable)
    # loop over variables to read from the ERA5 file
    these_labels = cfg_labels #[l for l in cfg_labels if "Precip" not in l]
    for label in these_labels:
        #d = data[site][label]["Data"]
        d = numpy.hstack(data[site][label]["Data"])
        # create the variable
        variable = {"Label": label,
                    "Data": d,
                    "Flag": numpy.where(numpy.ma.getmaskarray(d) == True, ones, zeros),
                    "Attr": data[site][label]["Attr"]}
        # put it into the data structure
        pfp_utils.CreateVariable(dss_ats[site], variable)
# conversion from cumulative values to 1 hour values
msg = " Converting from cumulative to 1 hour for Fsd, Fsu, Fld, Flu, Fh, Fe"
logger.info(msg)
for site in sites:
    # shortwave downwelling
    Fsd = pfp_utils.GetVariable(dss_ats[site], "Fsd")
    Fsd["Data"] = Fsd["Data"]/(era5_timestep*60.0)
    pfp_utils.CreateVariable(dss_ats[site], Fsd)
    # net shortwave
    Fnsw = pfp_utils.GetVariable(dss_ats[site], "Fnsw")
    Fnsw["Data"] = Fnsw["Data"]/(era5_timestep*60.0)
    pfp_utils.CreateVariable(dss_ats[site], Fnsw)
    # longwave downwelling
    # Fld = pfp_utils.GetVariable(dss_ats[site], "Fld")
    # Fld["Data"] = Fld["Data"]/(era5_timestep*60.0)
    # pfp_utils.CreateVariable(dss_ats[site], Fld)
    # net longwave
    Fnlw = pfp_utils.GetVariable(dss_ats[site], "Fnlw")
    Fnlw["Data"] = Fnlw["Data"]/(era5_timestep*60.0)
    pfp_utils.CreateVariable(dss_ats[site], Fnlw)
    # sensible heat flux
    # Fh = pfp_utils.GetVariable(dss_ats[site], "Fh")
    # Fh["Data"] = float(-1)*Fh["Data"]/(era5_timestep*60.0)
    # pfp_utils.CreateVariable(dss_ats[site], Fh)
    # latent heat flux
    # Fe = pfp_utils.GetVariable(dss_ats[site], "Fe")
    # Fe["Data"] = float(-1)*Fe["Data"]/(era5_timestep*60.0)
    # pfp_utils.CreateVariable(dss_ats[site], Fe)
    # Precipitation is in m, convert to mm
    # Precip = pfp_utils.GetVariable(dss_ats[site], "Precip")
    # Precip["Data"] = Precip["Data"]*float(1000)
    # pfp_utils.CreateVariable(dss_ats[site], Precip)
    # leaf area index
    LAI = pfp_utils.GetVariable(dss_ats[site], "LAI")
    pfp_utils.CreateVariable(dss_ats[site], LAI)
# do the units conversion
msg = " Converting units for ps, Ta, Td, Ts and Sws"
logger.info(msg)
for site in sites:
    # surface pressure units Pa to kPa
    # ps = pfp_utils.GetVariable(dss_ats[site], "ps")
    # ps = pfp_utils.convert_units_func(dss_ats[site], ps, "kPa")
    # pfp_utils.CreateVariable(dss_ats[site], ps)
    # air temperature
    Ta = pfp_utils.GetVariable(dss_ats[site], "Ta")
    Ta = pfp_utils.convert_units_func(dss_ats[site], Ta, "degC")
    pfp_utils.CreateVariable(dss_ats[site], Ta)
    # dew point temperature
    Td = pfp_utils.GetVariable(dss_ats[site], "Td")
    Td = pfp_utils.convert_units_func(dss_ats[site], Td, "degC")
    pfp_utils.CreateVariable(dss_ats[site], Td)
    # soil temperature level 1
    Ts = pfp_utils.GetVariable(dss_ats[site], "Ts")
    Ts = pfp_utils.convert_units_func(dss_ats[site], Ts, "degC")
    pfp_utils.CreateVariable(dss_ats[site], Ts)
    # soil temperature level 2
    # Ts2 = pfp_utils.GetVariable(dss_ats[site], "Ts2")
    # Ts2 = pfp_utils.convert_units_func(dss_ats[site], Ts2, "degC")
    # pfp_utils.CreateVariable(dss_ats[site], Ts2)
    # soil temperature level 3
    # Ts3 = pfp_utils.GetVariable(dss_ats[site], "Ts3")
    # Ts3 = pfp_utils.convert_units_func(dss_ats[site], Ts3, "degC")
    # pfp_utils.CreateVariable(dss_ats[site], Ts3)
    # soil temperature level 4
    # Ts4 = pfp_utils.GetVariable(dss_ats[site], "Ts4")
    # Ts4 = pfp_utils.convert_units_func(dss_ats[site], Ts4, "degC")
    # pfp_utils.CreateVariable(dss_ats[site], Ts4)
    # soil moisture level 1
    Sws = pfp_utils.GetVariable(dss_ats[site], "Sws")
    Sws = pfp_utils.convert_units_func(dss_ats[site], Sws, "m^3/m^3")
    pfp_utils.CreateVariable(dss_ats[site], Sws)
    # soil moisture level 2
    # Sws2 = pfp_utils.GetVariable(dss_ats[site], "Sws2")
    # #Sws2 = pfp_utils.convert_units_func(dss_ats[site], Sws2, "m^3/m^3")
    # pfp_utils.CreateVariable(dss_ats[site], Sws2)
    # soil moisture level 3
    # Sws3 = pfp_utils.GetVariable(dss_ats[site], "Sws3")
    # #Sws3 = pfp_utils.convert_units_func(dss_ats[site], Sws3, "m^3/m^3")
    # pfp_utils.CreateVariable(dss_ats[site], Sws3)
    # soil moisture level 4
    # Sws4 = pfp_utils.GetVariable(dss_ats[site], "Sws4")
    # #Sws4 = pfp_utils.convert_units_func(dss_ats[site], Sws4, "m^3/m^3")
    # pfp_utils.CreateVariable(dss_ats[site], Sws4)

# get the derived quantities
msg = " Deriving radiation, ground heat flux, absolute humidity, wind speed and direction"
logger.info(msg)
for site in sites:
    # calculate_upwelling_shortwave_radiation(dss_ats[site])
    # calculate_upwelling_longwave_radiation(dss_ats[site])
    calculate_net_radiation(dss_ats[site])
    # calculate_available_energy(dss_ats[site])
    # calculate_ground_heat_flux(dss_ats[site])
    calculate_relative_humidity(dss_ats[site])
    # calculate_absolute_humidity(dss_ats[site])
    # calculate_specific_humidity(dss_ats[site])
    calculate_ws_and_wd_from_u_and_v(dss_ats[site])
    calculate_vapour_pressure_deficit(dss_ats[site])
    delete_intermediate_variables(dss_ats[site], ["Fnlw", "Fnsw"])
    # discard duplicates or overlapping times found
    msg = " FixTimeStep: duplicate or overlapping times found, removing ..."
    logger.info(msg)
    pfp_utils.RemoveDuplicateRecords(dss_ats[site])

# interpolate from the ERA5 time step (60 minutes) to the tower time step as required
dss_tts = {}
msg = " Interpolating data for each site "
logger.info(msg)
for site in sites:
    era5_ts = int(float(dss_ats[site].root["Attributes"]["time_step"]))
    tower_ts = int(float(site_info[site]["Time step"]))
    labels = [l for l in list(dss_ats[site].root["Variables"].keys()) if "DateTime" not in l]
    # interpolate from the ERA5 time step (60 minutes) to the tower time step
    dss_tts[site] = pfp_ts.InterpolateDataStructure(dss_ats[site], labels=labels,
                                                    new_time_step=tower_ts,
                                                    sums="interpolate",
                                                    mode="quiet")

# write out the files of ERA5 data
existing_era5_base_path = cfg["Files"]["existing_era5_base_path"]
msg = " Writing ERA5 netCDF files"
logger.info(msg)
for site in sites:
    cis = concatenation_info[site]

    # What if in_file_name does not exist 
    main_in_file_name = ''.join(cis["NetCDFConcatenate"]["in_file_names"])
    if not os.path.isfile(main_in_file_name):
        cis["NetCDFConcatenate"]["in_file_names"].remove(main_in_file_name)
    # get the ERA5 data path
    era5_file_path = os.path.join(existing_era5_base_path, site, "Data", "ERA5", "processed")
    # check the directory exists and create if it doesn't
    if not os.path.isdir(era5_file_path):
        os.makedirs(era5_file_path)
    # build the ERA5 file name
    dt = pfp_utils.GetVariable(dss_tts[site], "DateTime")
    start_date = dt["Data"][0].strftime("%Y%m%d%H%M")
    end_date = dt["Data"][-1].strftime("%Y%m%d%H%M")
    era5_file_name = site + "_ERA5" + "_" + start_date + "_" + end_date + ".nc"
    # get the full path including the file name
    era5_file_uri = os.path.join(era5_file_path, era5_file_name)
    # and write the ERA5 data to a netCDF file
    pfp_io.NetCDFWrite(era5_file_uri, dss_tts[site])
    cis["NetCDFConcatenate"]["in_file_names"].append(era5_file_uri)
    # concatenate with the existing ERA5 data
    pfp_io.NetCDFConcatenate(cis)
    # save the concatenation_info dictionary
    cis_file_name = "cis_" + site + ".txt"
    cis_file_path = os.path.join(existing_era5_base_path, site, "Data", "ERA5", "cis")
    if not os.path.isdir(cis_file_path):
        os.mkdir(cis_file_path)
    cfg_cis = ConfigObj(cis, indent_type="    ", list_values=False)
    cfg_cis.filename = os.path.join(cis_file_path, cis_file_name)
    cfg_cis.write()

logger.info("All finished")
