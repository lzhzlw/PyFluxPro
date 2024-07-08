import netCDF4
import os
import pandas

def read_netcdf(nc_full_name,variable_list=[]):
    """
    Purpose:
     Read an OzFlux netCDF file and return the data in an Pandas data frame.
    Usage:
     df = qcio.nc_read_todf(nc_full_name,variable_list=variable_list)
      where nc_full_name is the full name of the netCDF file.
            variable_list (optional) is an list of variables to be read
     If variable_list is not passed, all variables in the netCDF are returned.
    Side effects:
     Returns a Pandas data frame containing the data indexed by datetime and
     a dictionary containing the global and variable attributes.
    Author: PRI using code originally written by Ian McHugh
    Date: June 2015
    """
    # check to see if the file exists
    if "http" not in nc_full_name.lower():
        if not os.path.exists(nc_full_name):
            raise Exception("read_netcdf: input file "+nc_full_name+" not found")
    # read the netCDF file
    nc_file = netCDF4.Dataset(nc_full_name,"r")
    # create a dictionary to hold the global and variable attributes
    attr = {}
    attr["global"] = {}
    attr["variable"] = {}
    # now deal with the global attributes
    gattrlist = nc_file.ncattrs()
    if len(gattrlist)!=0:
        for item in gattrlist:
            attr["global"][item] = getattr(nc_file,item)
    # get a list of Python datetimes from the xlDatetime
    time = time = nc_file.variables["time"][:]
    time_units = getattr(nc_file.variables["time"],"units")
    dates_list = list(netCDF4.num2date(time,time_units))
    # get a list of variables to read from the netCDF file
    # was a variable list passed in as variable_list?
    if len(variable_list)==0:
        # if not, get the variable list from the netCDF file contents
        variable_list = list(nc_file.variables.keys())
    else:
        # if so, add the QC flags to the list entered as an argument
        flag_list = []
        for item in variable_list: flag_list.append(item+"_QCFlag")
        variable_list = variable_list+flag_list
    # read the variables and attributes from the netCDF file
    # create a dictionary to hold the data
    data = {}
    # loop over the variables to be read
    for item in variable_list:
        # get the number of dimensions
        # variables in OzFlux netCDF files can have 1 (time) or 3 dimensions (time,latitude,longitude)
        ndims = len(nc_file.variables[item].shape)
        if ndims==1:
            data[item] = nc_File.variables[item][:]
        elif ndims==3:
            # drop the degenerate dimensions (latitude and longitude)
            data[item] = nc_file.variables[item][:,0,0]
            print((item,ndims,data[item]))
        else:
            raise Exception("unrecognised number of dimensions for variable"+str(item))
        # get the variable attributes
        vattrlist = nc_file.variables[item].ncattrs()
        if len(vattrlist)!=0:
            attr["variable"][item] = {}
            for vattr in vattrlist:
                attr["variable"][item][vattr] = getattr(nc_file.variables[item],vattr)
    nc_file.close()
    # convert the dictionary to a Pandas data frame
    df = pandas.DataFrame(data,index=dates_list)
    return df,attr

# ==== MAIN ==== #
in_list = ['Fe','Fh','Fc']
#['Fe','Fh','Fc','Ah','Ah_15m','Ah_2m','CO2','VbatCR3000','VbatCR1000','TpanelCR1000','TpanelCR3000','Fsd',\
#'Fn_NR','ps','Precip','RH_15m','RH_2m','Sws_10cm','Sws_10cm','Sws_20cm','Sws_20cm','Sws_160cm',\
#'Ta_15m','Ta_2m','Ts_200cm','Ts_10cm','Ts_4cm','IR_TargT','Wd_WS4_16m','Ws_WS4_16m','Ws_WS4_2m']

#out_list = \
#['Fe_raw_Avg','Fh_Avg','Fc_raw_Avg','Ah_7500_Avg','Ah_HMP_15m_Avg','Ah_HMP_2m_Avg','Cc_7500_Avg','CR3000_batt_volt_Avg','CR1000_batt_volt_Avg','CR1000_PTemp_Avg','CR3000_PTemp_Avg','Fsd_CNR1_Avg',
#'NRLITE_Wm2_Avg','ps_7500_Avg','Rain_Tot','RH_HMP_15m_Avg','RH_HMP_2m_Avg','Sws_616_NW_10cm_a_Avg','Sws_616_W_10cm_a_Avg','Sws_616_W_20cm_Avg','Sws_616_NW_20cm_Avg','Sws_616_NW_160cm_a_Avg',
#'Ta_HMP_15m_Avg','Ta_HMP_2m_Avg','Sws_CS650_mux1_2m_Temp_Avg','Ts_T108_mux1_10cm_Avg','Ts_TCAVa_4cm_Avg','TargTempC_Avg','WD_WS4_16m_Avg_Avg','WS_WS4_16m_Avg_Avg','WS_WS4_2m_Avg_Avg']

# read the variables from the local netCDF file
nc_full_name = "../../Sites/Loxton/Data/Processed/all/Loxton_2008_to_2009_L3.nc"
variable_list = in_list
print ("reading local netCDF file")
df,attr = read_netcdf(nc_full_name,variable_list=variable_list)

