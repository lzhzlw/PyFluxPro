#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 15:16:21 2019

@author: ian
"""

#------------------------------------------------------------------------------
### MODULES (STANDARD) ###
#------------------------------------------------------------------------------

import datetime as dt
import glob
import numpy as np
import os
import pandas as pd
from pytz import timezone
import xarray as xr

#------------------------------------------------------------------------------
### MODULES (CUSTOM) ###
#------------------------------------------------------------------------------

import utils
import met_funcs

#------------------------------------------------------------------------------
### CONFIGURATIONS ###
#------------------------------------------------------------------------------

configs = utils.get_configs()
conventions = utils.get_conventions()
raw_file_path = configs['raw_data_write_paths']['access']
#raw_file_path_prev = configs['raw_data_write_paths']['access_previous']
raw_file_path_prev = '/rdsi/market/access_old_site_files/monthly'
nc_write_path = configs['nc_data_write_paths']['access']

#------------------------------------------------------------------------------
### CLASSES ###
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
class access_data_converter():

    """Conversion class for converting raw ACCESS data into PFP format
       (note that previously, ACCESS data was collected and partially 
       converted before dumping the raw data... tsk tsk)
       Args:
           * site_details (pandas dataframe): df containing the details of the
             site as documented in the site master file
           * include_prior_data (boolean): if True, concatenates the older 
             semi-formatted (see above) data as well as the new data (which
             preserves the format and attributes of th original ACCESS nc file
    """
    
    def __init__(self, site_details, include_prior_data=False):

        self.site_details = site_details
        self.include_prior_data = include_prior_data

    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def create_dataset(self):

        """Main method for concatenation and conversion of raw data to PFP
           format"""
        
        # Do formatting, conversions and attributes
        ds = self.get_raw_file()
        ds = ds.compute() # This converts from dask to numpy :)
        ds = ds[list(vars_dict.keys())]
        _apply_range_limits(ds)
        ds = ds.rename(vars_dict)
        ds = _reindex_time(ds)
        if self.site_details['Time step'] == 30: ds = _resample_dataset(ds)
        do_conversions(ds)
        get_energy_components(ds)
        ds = ds.fillna(-9999.0)
        offset = self.get_utc_offset()
        ds.time.data = (pd.to_datetime(ds.time.data) + dt.timedelta(hours=offset))
        _set_var_attrs(ds)

        # Rebuild dataset with separated pixels
        ds_list = []
        for i, this_lat in enumerate(ds.lat):
            for j, this_lon in enumerate(ds.lon):

                # Deal with dimensions
                sub_ds = ds.sel(lat=ds.lat[i], lon=ds.lon[j])
                for x in ['Ts', 'Sws']:
                    sub_ds[x] = sub_ds[x].sel(soil_lvl=sub_ds.soil_lvl[0])

                # Dump extraneous dims and coords
                sub_ds = sub_ds.drop_dims(['soil_lvl'])
                sub_ds = sub_ds.reset_coords(['lat', 'lon'], drop=True)

                # Add lat long to variable attributes
                for var in sub_ds:
                    sub_ds[var].attrs['latitude'] = round(this_lat.item(), 4)
                    sub_ds[var].attrs['longitude'] = round(this_lon.item(), 4)

                # Rename with variable numbering
                var_suffix = '_{}{}'.format(str(i), str(j))
                new_dict = {x: x + var_suffix for x in list(sub_ds.variables)
                            if not x == 'time'}
                sub_ds = sub_ds.rename(new_dict)

                # Append
                ds_list.append(sub_ds)

        # Merge qc flags
        final_ds = xr.merge(ds_list)
        final_ds = xr.merge([final_ds, make_qc_flags(final_ds)])
        
        # Add old data (or not), set global attrs and return
        if self.include_prior_data:
            try: final_ds = _combine_datasets(final_ds, self.site_details.name)
            except OSError: pass
        
        # Insert old variable names temporarily
        final_ds = xr.merge([final_ds, _make_temp_vars(final_ds)])
        _set_global_attrs(final_ds, self.site_details)
        
        return final_ds

    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def get_file_list(self):

        """Return list of files found in the path for the requested site"""
        
        search_str = self.site_details.name.replace(' ', '')
        return sorted(glob.glob(raw_file_path +
                                '/Monthly_files/**/{}*'.format(search_str)))
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def get_raw_file(self):

        """Concatenate all files found in the path, dropping any dupe indices
        """
        
        def preproc(ds):
            idx = np.unique(ds.time.data, return_index=True)[1]
            return ds.isel(time=idx)

        return xr.open_mfdataset(self.get_file_list(), combine='by_coords',
                                 preprocess=preproc)
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def get_utc_offset(self):

        """Get the UTC offset of the site"""
        
        tz_obj = timezone(self.site_details['Time zone'])
        now_time = dt.datetime.now()
        return (tz_obj.utcoffset(now_time) - tz_obj.dst(now_time)).seconds / 3600
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def write_to_netcdf(self, write_path):

        print('Writing netCDF file for site {}'.format(self.site_details.name))
        dataset = self.create_dataset()
        fname = '{}_ACCESS.nc'.format(self.site_details.name.replace(' ', ''))
        target = os.path.join(write_path, fname)
        dataset.to_netcdf(target, format='NETCDF4')
    #--------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
### FUNCTIONS ###
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _apply_range_limits(ds):

    """Apply variable range limits as documented in range_dict"""
    
    for var in ds.variables:
        if var in ds.dims: continue
        lims = range_dict[var]
        ds[var] = ds[var].where(cond=(ds[var] >= lims[0]) & (ds[var] <= lims[1]))
    ds['inst_prcp'] = ds.inst_prcp.where(cond=((ds.inst_prcp < -1) |
                                               (ds.inst_prcp > 0.001)), other=0)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _collate_prior_data(site_str):
    
    """Stack the individual data files from the old ACCESS collection together;
       note that we suspect the rainfall calculation for the mod6 = 0 hours is
       wrong in this data, so we cut to the new collection in November of '19,
       the first month for which we collected complete data"""
    
    def preproc(sub_ds):
        return sub_ds.drop_sel(time=sub_ds.time[-1].data)

    cutoff_month = '201910'
    f_list = sorted(glob.glob(raw_file_path_prev + '/**/{}*'.format(site_str)))
    if len(f_list) == 0: raise OSError
    months = [os.path.splitext(x)[0].split('_')[-1] for x in f_list]
    for i, month in enumerate(months): 
        if month > cutoff_month: break
    f_list = f_list[:i]
    prior_ds = xr.open_mfdataset(f_list, combine='by_coords', preprocess=preproc).compute()
    for var in list(prior_ds.variables):
        prior_ds[var].encoding['_FillValue'] = -9999.0
        try: prior_ds[var].encoding.pop('missing_value')
        except KeyError: pass
        if not 'Precip' in var: continue
        prior_ds[var] = prior_ds[var].where(cond=prior_ds[var]<9999, other=-9999.0)
    return prior_ds
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _combine_datasets(current_ds, site_name):
    
    """Collate the newer data with the older data"""
    
    site_name = site_name.replace(' ','')
    prior_ds = _collate_prior_data(site_name)
    prior_ds = prior_ds.drop(labels=[x for x in prior_ds.variables 
                                     if not x in current_ds.variables])
    for var in current_ds.variables:
        if var in current_ds.dims: continue
        prior_ds[var].attrs = current_ds[var].attrs
    current_ds = xr.concat([prior_ds, current_ds], dim='time')
    idx = np.unique(current_ds.time.data, return_index=True)[1]
    current_ds = current_ds.isel(time=idx)
    return current_ds
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------
def do_conversions(ds):

    """Convert to desired units for output to PFP format files"""
    
    ds['Ws'] = met_funcs.get_ws_from_uv(ds.u, ds.v)
    ds['Wd'] = met_funcs.get_wd_from_uv(ds.u, ds.v)
    ds['ps'] = met_funcs.convert_Pa_to_kPa(ds.ps)
    ds['RH'] = (met_funcs.get_e_from_q(ds.SH, ds.ps) / 
                met_funcs.get_es(ds.Ta)) * 100
    ds['AH'] = met_funcs.get_Ah(ds.Ta, ds.SH, ds.ps)
    ds['Ta'] = met_funcs.convert_Kelvin_to_celsius(ds.Ta)
    ds['Ts'] = met_funcs.convert_Kelvin_to_celsius(ds.Ts) 
    ds['Sws'] = ds['Sws'] / 100
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_energy_components(ds):

    """Calculate energy balance components"""
    
    ds['Fsu'] = ds.Fsd - ds.Fn_sw
    ds['Flu'] = ds.Fld - ds.Fn_lw
    ds['Fn'] = (ds.Fsd - ds.Fsu) + (ds.Fld - ds.Flu)
    ds['Fa'] = ds.Fh + ds.Fe
    ds['Fg'] = ds.Fn - ds.Fa
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def make_qc_flags(ds):

    """Generate QC flags for all variables in the ds"""

    da_list = []
    for var in ds.variables:
        if var in ds.dims: continue
        da = xr.where(~np.isnan(ds[var]), 0, 10)
        da.name = var + '_QCFlag'
        da_list.append(da)
    return xr.merge(da_list)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _make_temp_vars(ds):

    """Temporary addition of old variable names"""
    
    ds_list = []
    rename_dict = {'AH': 'Ah', 'SH': 'q'}
    for this_str in rename_dict.keys():
        vars_list = [x for x in ds.variables if this_str in x]
        vars_mapper = {var: var.replace(this_str, rename_dict[this_str]) 
                       for var in vars_list}
        sub_ds = ds[vars_list].rename(vars_mapper)
        ds_list.append(sub_ds)
    return xr.merge(ds_list)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _reindex_time(ds):

    """Reindex dataset to include missing cases"""
    
    new_index = pd.date_range(ds.time[0].item(), ds.time[-1].item(), freq='60T')
    return ds.reindex(time=new_index)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _resample_dataset(ds):

    """Resample to half-hourly and interpolate only gaps created by
       resampling; note that rainfall must be converted to a cumulative sum
       to be resampled, then redifferenced to recover hourly total rainfall"""

    new_ds = ds.copy()
    new_ds['cml_precip'] = new_ds.Precip.cumsum(dim='time')
    new_dates = pd.date_range(start=ds.time[0].item(), end=ds.time[-1].item(),
                              freq='30T')
    new_ds = (new_ds.reindex(time=new_dates)
              .interpolate_na(dim='time', max_gap=pd.Timedelta(hours=1)))
    new_ds['cml_precip'] = new_ds.cml_precip.where(~np.isnan(new_ds.Precip))
    new_ds['Precip'] = new_ds.cml_precip - new_ds.cml_precip.shift(time=1)
    new_ds = new_ds.drop('cml_precip')
    return new_ds
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _set_var_attrs(ds):

    for this_var in list(ds.variables):
        if this_var in ds.dims: continue
        ds[this_var].attrs = _get_var_attrs(this_var)
        ds[this_var].encoding = {'_FillValue': -9999.0}
#------------------------------------------------------------------------------
        
#------------------------------------------------------------------------------        
def _get_var_attrs(var):
    
    unit_list_level_dict = {'RH': 1, 'Ta': 0, 'ps': 0, 'Ts': 0}
    if var in ['crs', 'time']:
        base_dict = {}
    else:
        base_dict = generic_dict.copy()
    base_dict.update(attrs_dict[var])
    try:
        standard_attrs = conventions['variable_attributes'][var]
        if isinstance(standard_attrs['units'], list):
            idx = unit_list_level_dict[var]
            standard_attrs['units'] = standard_attrs['units'][idx]
        base_dict.update(standard_attrs)
    except KeyError:
        pass
    return base_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _set_global_attrs(ds, site_details):

    ds.attrs = {'nc_nrecs': len(ds.time),
                'nc_level': 'L1',
                'nc_rundatetime': (dt.datetime.now().
                                   strftime('%Y-%m-%d %H:%M:%S')),
                'start_date': (dt.datetime.strftime
                               (pd.Timestamp(ds.time[0].item()),
                                '%Y-%m-%d %H:%M:%S')),
                'end_date': (dt.datetime.strftime
                             (pd.Timestamp(ds.time[-1].item()),
                              '%Y-%m-%d %H:%M:%S')),
                'latitude': site_details.Latitude,
                'longitude': site_details.Longitude,
                'site_name': site_details.name,
                'time_step': str(int(site_details['Time step'])),
                'time_zone': site_details['Time zone'],
                'xl_datemode': '0'}
    ds.time.encoding = {'units': 'days since 1800-01-01',
                        '_FillValue': None}
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------
### GLOBALS ###
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
attrs_dict = {'AH': {'source': 'Calculated'}, 
              'Fa': {'source': 'Calculated'},
              'Fe': {'level_type': 'single', 'accum_type': 'instantaneous', 
                     'stash_code': '3234', 'grid_type': 'spatial', 
                     'source': 'access'}, 
              'Fg': {'source': 'Calculated'}, 
              'Fh': {'level_type': 'single', 'accum_type': 'instantaneous', 
                     'stash_code': '3217', 'grid_type': 'spatial',
                     'source': 'ACCESS'},
              'Fld': {'level_type': 'single', 'accum_type': 'mean',
                      'stash_code': '2207', 'accum_units': 'hrs',
                      'accum_value': '1', 'grid_type': 'spatial',
                      'source': 'ACCESS'},
              'Flu': {'source': 'Calculated'},
              'Fn': {'source': 'Calculated'},
              'Fn_lw': {'long_name': 'Net longwave radiation', 
                        'group_name': 'radiation', 'units': 'W/m^2',
                        'level_type': 'single', 'accum_type': 'mean', 
                        'stash_code': '2201', 'accum_units': 'hrs', 
                        'accum_value': '1', 'grid_type': 'spatial', 
                        'source': 'ACCESS'},
              'Fn_sw': {'long_name': 'Net shortwave radiation', 
                        'group_name': 'radiation', 'units': 'W/m^2',
                        'level_type': 'single', 'coverage_L1': '66',
                        'accum_type': 'mean', 'stash_code': '1202',
                        'accum_units': 'hrs', 'accum_value': '1',
                        'grid_type': 'spatial', 'source': 'ACCESS'},
              'Fsd': {'level_type': 'single', 'accum_type': 'mean', 
                      'stash_code': '1235', 'accum_units': 'hrs', 
                      'accum_value': '1', 'grid_type': 'spatial',
                      'source': 'ACCESS'},
              'Fsu': {'source': 'Calculated'},
              'Habl': {'level_type': 'single', 
                       'long_name': 'Planetary boundary layer height',
                       'accum_type': 'instantaneous', 'stash_code': '25', 
                       'units': 'm', 'grid_type': 'spatial', 
                       'source': 'ACCESS'},
              'Precip': {'level_type': 'single', 'accum_type': 'accumulative', 
                         'stash_code': '5226', 'accum_units': 'hrs', 
                         'accum_value': '3', 'grid_type': 'spatial',
                         'source': 'ACCESS'},
              'RH': {'source': 'Calculated'},
              'Sws': {'accum_type': 'instantaneous', 'stash_code': '8223', 
                      'level_type': 'multi', 'grid_type': 'spatial',
                      'source': 'ACCESS'},
              'Ta': {'accum_type': 'instantaneous', 'stash_code': '3236', 
                     'level_type': 'single', 'grid_type': 'spatial',
                     'source': 'ACCESS'},
              'Ts': {'accum_type': 'instantaneous', 'stash_code': '8225', 
                     'level_type': 'multi', 'grid_type': 'spatial',
                     'source': 'ACCESS'},
              'Wd': {'height': '10m', 'source': 'Calculated'},
              'Ws': {'height': '10m', 'source': 'Calculated'},
              'ps': {'accum_type': 'instantaneous', 'stash_code': '409', 
                     'level_type': 'single', 'grid_type': 'spatial',
                     'source': 'ACCESS'},
              'crs': {'grid_mapping_name': 'latitude_longitude',
                      'long_name': 'WGS 1984 datum',
                      'longitude_of_prime_meridian': '0.0',
                      'semi_major_axis': '6378137.0',
                      'inverse_flattening': '298.257223563',
                      'source': ''},
              'SH': {'level_type': 'single', 'accum_type': 'instantaneous', 
                     'stash_code': '3237', 'grid_type': 'spatial',
                     'source': 'ACCESS'},
              'time': {'long_name': 'time', 'standard_name': 'time',
                       'source': ''},
              'u': {'level_type': 'single', 'long_name': 'Wind u component',
                    'accum_type': 'instantaneous', 'stash_code': '3209', 
                    'units': 'm/s', 'height': '10m', 'grid_type': 'spatial',
                    'source': 'ACCESS'},
              'ustar': {'level_type': 'single', 'accum_type': 'instantaneous', 
                        'stash_code': '3465', 'grid_type': 'spatial', 
                        'source': 'ACCESS'},
              'v': {'level_type': 'single', 'long_name': 'Wind v component',
                    'accum_type': 'instantaneous', 'stash_code': '3210', 
                    'units': 'm/s', 'height': '10m', 'grid_type': 'spatial',
                    'source': 'ACCESS'}
              }
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
range_dict = {'av_swsfcdown': [0, 1400],
              'av_netswsfc': [0, 1400],
              'av_lwsfcdown': [200, 600],
              'av_netlwsfc': [-300, 300],
              'temp_scrn': [230, 330],
              'qsair_scrn': [0, 1],
              'soil_mois': [0, 100],
              'soil_temp': [210, 350],
              'u10': [-50, 50],
              'v10': [-50, 50],
              'sfc_pres': [75000, 110000],
              'inst_prcp': [-1, 100],
              'sens_hflx': [-200, 1000],
              'lat_hflx': [-200, 1000],
              'abl_ht': [0, 5000]}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
vars_dict = {'av_swsfcdown': 'Fsd',
             'av_netswsfc': 'Fn_sw',
             'av_lwsfcdown': 'Fld',
             'av_netlwsfc': 'Fn_lw',
             'temp_scrn': 'Ta',
             'qsair_scrn': 'SH',
             'soil_mois': 'Sws',
             'soil_temp': 'Ts',
             'u10': 'u',
             'v10': 'v',
             'sfc_pres': 'ps',
             'inst_prcp': 'Precip',
             'sens_hflx': 'Fh',
             'lat_hflx': 'Fe',
             'abl_ht': 'Habl'}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
generic_dict = {'height': 'not defined',
                'standard_name': 'not defined',
                'valid_range': '-1e+35,1e+35'}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
### MAIN PROGRAM ###
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
if __name__ == "__main__":

    sites = utils.get_ozflux_site_list()#(active_sites_only=True)
    for site in sites.index[1:]:
        specific_file_path = nc_write_path.format(site.replace(' ', ''))
        converter = access_data_converter(sites.loc[site], include_prior_data=True)
        converter.write_to_netcdf(specific_file_path)
#------------------------------------------------------------------------------