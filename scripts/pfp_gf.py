# standard modules
import copy
import datetime
import os
import logging
# 3rd party modules
import dateutil
import numpy
import xlrd
# PFP modules
from scripts import constants as c
from scripts import pfp_gui
from scripts import pfp_io
from scripts import pfp_ts
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def CheckL5Drivers(ds, l5_info):
    """
    Purpose:
     Check the drvers specified for gap filling for missing data.
    Usage:
    Side effects:
    Author: PRI
    Date: October 2020
    """
    msg = " Checking drivers for missing data"
    logger.info(msg)
    drivers_with_missing = []
    for label in list(l5_info["CheckL5Drivers"]["drivers"]):
        if label not in ds.root["Variables"].keys():
            msg = "  Requested driver (" + label + ") not found in data structure"
            logger.error(msg)
            ds.info["returncodes"]["message"] = msg
            ds.info["returncodes"]["value"] = 1
            return
        var = pfp_utils.GetVariable(ds, label)
        if numpy.any(numpy.ma.getmaskarray(var["Data"])):
            drivers_with_missing.append(label)
    if len(drivers_with_missing) == 0:
        msg = "  No missing data found in drivers"
        logger.info(msg)
        ds.info["returncodes"] = {"value": 0, "message": msg}
    else:
        dwm = ",".join(drivers_with_missing)
        msg = " Drivers " + dwm + " have missing data, aborting L5 ..."
        logger.error("!!!!!")
        logger.error(msg)
        logger.error("!!!!!")
        ds.info["returncodes"] = {"value": 1, "message": msg}
    return

def CheckGapLengths(cfg, ds, l5_info):
    """
    Purpose:
     Check to see if any of the series being gap filled have long gaps and
     tell the user if they are present.  The user can then opt to continue
     or to add a long-gap filling method (e.g. GapFillLongSOLO) to the
     control file.
    Usage:
    Side effects:
    Author: PRI
    Date: June 2019
    """
    ds.info["returncodes"]["value"] = 0
    l5_info["CheckGapLengths"] = {}
    # get the maximum length for "short" gaps in days
    opt = pfp_utils.get_keyvaluefromcf(cfg, ["Options"], "MaxShortGapDays", default=30)
    max_short_gap_days = int(opt)
    # maximum length in records
    ts = int(float(ds.root["Attributes"]["time_step"]))
    nperday = 24 * 60//ts
    max_short_gap_records = max_short_gap_days * nperday
    targets_with_long_gaps = []
    # loop over the targets, get the duration and check to see if any exceed the maximum
    targets = l5_info["CheckL5Targets"]["targets"]
    for target in targets:
        # initialise dictionary entry
        l5_info["CheckGapLengths"][target] = {"got_long_gaps": False,
                                              "got_long_gap_method": False}
        # loop over possible long gap filling methods
        for long_gap_method in ["GapFillLongSOLO"]:
            if long_gap_method in list(cfg["Fluxes"][target].keys()):
                # set logical true if long gap filling method present
                l5_info["CheckGapLengths"][target]["got_long_gap_method"] = True
        # get the data
        variable = pfp_utils.GetVariable(ds, target)
        # get the mask
        mask = numpy.ma.getmaskarray(variable["Data"])
        # get the gaps
        gap_start_end = pfp_utils.contiguous_regions(mask)
        # loop over the gaps
        for start, end in gap_start_end:
            gap_length = end - start
            # check to see if any gaps are longer than the max
            if gap_length > max_short_gap_records:
                # set logical if long gaps present
                l5_info["CheckGapLengths"][target]["got_long_gaps"] = True
                targets_with_long_gaps.append(target)
                break
    # write an info message to the log window
    if len(targets_with_long_gaps) != 0:
        msg = " Series " + ",".join(targets_with_long_gaps) + " have gaps longer than "
        msg = msg + str(max_short_gap_days) + " days"
        logger.warning(msg)
    # check for targets with long gaps but no long gap fill method
    targets_without = []
    for target in targets:
        if (l5_info["CheckGapLengths"][target]["got_long_gaps"] and
            not l5_info["CheckGapLengths"][target]["got_long_gap_method"]):
            targets_without.append(target)
    # if we have any, put up a warning message and let the user decide
    if len(targets_without) != 0:
        if cfg["Options"]["call_mode"].lower() == "interactive":
            # put up a message box, continue or quit
            msg = "The following series have long gaps but no long gap filling method\n"
            msg = msg + "is specified in the control file.\n"
            msg = msg + "    " + ",".join(targets_without) + "\n"
            msg = msg + "To add a long gap fill method, press 'Quit' and edit the control file\n"
            msg = msg + "or press 'Continue' to ignore this warning."
            result = pfp_gui.MsgBox_ContinueOrQuit(msg, title="Warning: Long gaps")
            if result.clickedButton().text() == "Quit":
                # user wants to edit the control file
                msg = " Quitting L5 to edit control file"
                logger.warning(msg)
                ds.info["returncodes"]["message"] = msg
                ds.info["returncodes"]["value"] = 1
            else:
                # user wants to continue, turn on auto-complete for SOLO ...
                if "GapFillUsingSOLO" in l5_info:
                    l5_info["GapFillUsingSOLO"]["gui"]["auto_complete"] = True
                # ... and disable masking of long gaps with MDS
                if "GapFillUsingMDS" in l5_info:
                    l5_info["GapFillUsingMDS"]["info"].pop("MaxShortGapRecords", None)
        else:
            # batch mode, turn on auto-complete for SOLO ...
            if "GapFillUsingSOLO" in l5_info:
                l5_info["GapFillUsingSOLO"]["gui"]["auto_complete"] = True
            # ... and disable masking of long gaps with MDS
            if "GapFillUsingMDS" in l5_info:
                l5_info["GapFillUsingMDS"]["info"].pop("MaxShortGapRecords", None)
    return

def CheckL5Targets(ds, l5_info):
    """
    Purpose:
     Check the targets specified for gap filling at L5 to see if any of them
     still contain missing data.
    Usage:
    Side effects:
    Author: PRI
    Date: October 2020
    """
    msg = " Checking targets for missing data"
    logger.info(msg)
    # get a list of target variables
    targets = l5_info["CheckL5Targets"]["targets"]
    series_with_missing_data = []
    for target in targets:
        var = pfp_utils.GetVariable(ds, target)
        if numpy.any(numpy.ma.getmaskarray(var["Data"])):
            series_with_missing_data.append(target)
    if len(series_with_missing_data) == 0:
        msg = "  No missing data found in targets"
        logger.info(msg)
        ds.info["returncodes"] = {"value": 0, "message": msg}
    else:
        s = ",".join(series_with_missing_data)
        msg = " Targets " + s + " contain missing data, aborting L5 ..."
        logger.error("!!!!!")
        logger.error(msg)
        logger.error("!!!!!")
        ds.info["returncodes"] = {"value": 1, "message": msg}
    return

def ParseL4ControlFile(cfg, ds):
    """
    Purpose:
     Create the L4 information and setting dictionary.
    Usage:
    Side effects:
    Author: PRI
    Date: Back in the day
    """
    ds.info["returncodes"]["message"] = "OK"
    ds.info["returncodes"]["value"] = 0
    l4_info = {}
    l4_info["cfg"] = copy.deepcopy(cfg)
    # add key for suppressing output of intermediate variables e.g. Ta_aws
    opt = pfp_utils.get_keyvaluefromcf(cfg, ["Options"], "KeepIntermediateSeries", default="No")
    l4_info["RemoveIntermediateSeries"] = {"KeepIntermediateSeries": opt, "not_output": []}
    # add key for interpolation
    interpolate_type = pfp_utils.get_keyvaluefromcf(cfg, ["Options"], "InterpolateType", default="Akima")
    max_gap_interpolate = pfp_utils.get_keyvaluefromcf(cfg, ["Options"], "MaxGapInterpolate", default=3)
    l4_info["GapFillUsingInterpolation"] = {"InterpolateType": str(interpolate_type),
                                            "MaxGapInterpolate": int(max_gap_interpolate)}
    # get a list of keys in the control file
    labels = list(cfg["Drivers"].keys())
    # get a list of targets being gap filled
    targets = sorted(list(cfg["Drivers"].keys()))
    l4_info["GapFillUsingInterpolation"]["targets"] = targets.copy()
    # loop over target variables
    for label in labels:
        if "GapFillFromAlternate" in list(cfg["Drivers"][label].keys()):
            gfalternate_createdict(cfg, ds, l4_info, label, "GapFillFromAlternate")
            # check to see if something went wrong
            if ds.info["returncodes"]["value"] != 0:
                # if it has, return to calling routine
                return l4_info
        if "GapFillFromClimatology" in list(cfg["Drivers"][label].keys()):
            gfClimatology_createdict(cfg, ds, l4_info, label, "GapFillFromClimatology")
            if ds.info["returncodes"]["value"] != 0:
                return l4_info
        if "MergeSeries" in list(cfg["Drivers"][label].keys()):
            gfMergeSeries_createdict(cfg, ds, l4_info, label, "MergeSeries")
    # check to make sure at least 1 output is defined
    outputs = []
    for method in ["GapFillFromAlternate", "GapFillFromClimatology"]:
        if method in list(l4_info.keys()):
            outputs = outputs + list(l4_info[method]["outputs"].keys())
    if len(outputs) == 0:
        msg = " No output variables defined, quitting L4 ..."
        logger.error(msg)
        ds.info["returncodes"]["message"] = msg
        ds.info["returncodes"]["value"] = 1
    return l4_info

def ParseL5ControlFile(cfg, ds):
    """
    Purpose:
     Create the L5 information and setting dictionary.
    Usage:
    Side effects:
    Author: PRI
    Date: Back in the day
    """
    ds.info["returncodes"]["message"] = "OK"
    ds.info["returncodes"]["value"] = 0
    l5_info = {}
    l5_info["cfg"] = copy.deepcopy(cfg)
    # add key for suppressing output of intermediate variables e.g. Ta_aws
    opt = pfp_utils.get_keyvaluefromcf(cfg, ["Options"], "KeepIntermediateSeries", default="No")
    l5_info["RemoveIntermediateSeries"] = {"KeepIntermediateSeries": opt, "not_output": []}
    # add key for interpolation
    interpolate_type = pfp_utils.get_keyvaluefromcf(cfg, ["Options"], "InterpolateType", default="Akima")
    max_gap_interpolate = pfp_utils.get_keyvaluefromcf(cfg, ["Options"], "MaxGapInterpolate", default=3)
    l5_info["GapFillUsingInterpolation"] = {"InterpolateType": str(interpolate_type),
                                            "MaxGapInterpolate": int(max_gap_interpolate)}
    # get a list of targets being gap filled, this may not simply be a list of keys
    # in the control file if 'target' is specified in the gap filling method
    targets = list(cfg["Fluxes"].keys())
    gf_methods = []
    for target in list(targets):
        gf_methods = [m for m in cfg["Fluxes"][target].keys() if m != "MergeSeries"]
        gf_method_targets = []
        for gf_method in list(gf_methods):
            gf_method_labels = list(cfg["Fluxes"][target][gf_method].keys())
            for gf_method_label in list(gf_method_labels):
                if "target" in cfg["Fluxes"][target][gf_method][gf_method_label]:
                    real_target = cfg["Fluxes"][target][gf_method][gf_method_label]["target"]
                    gf_method_targets.append(real_target)
                    targets[targets.index(target)] = real_target
    targets = targets + gf_method_targets
    targets = sorted(list(set(targets)))
    l5_info["GapFillUsingInterpolation"]["targets"] = targets.copy()
    l5_info["CheckL5Targets"] = {"targets": targets.copy()}
    l5_info["CheckL5Drivers"] = {"drivers": []}
    # get a list of keys in the control file
    labels = sorted(list(cfg["Fluxes"].keys()))
    for label in labels:
        if "GapFillUsingSOLO" in list(cfg["Fluxes"][label].keys()):
            gfSOLO_createdict(cfg, ds, l5_info, label, "GapFillUsingSOLO", 510)
            # check to see if something went wrong
            if ds.info["returncodes"]["value"] != 0:
                # if it has, return to calling routine
                return l5_info
        if "GapFillLongSOLO" in list(cfg["Fluxes"][label].keys()):
            gfSOLO_createdict(cfg, ds, l5_info, label, "GapFillLongSOLO", 520)
            if ds.info["returncodes"]["value"] != 0:
                return l5_info
        if "GapFillUsingMDS" in list(cfg["Fluxes"][label].keys()):
            gfMDS_createdict(cfg, ds, l5_info, label, "GapFillUsingMDS", 530)
            if ds.info["returncodes"]["value"] != 0:
                return l5_info
        if "MergeSeries" in list(cfg["Fluxes"][label].keys()):
            gfMergeSeries_createdict(cfg, ds, l5_info, label, "MergeSeries")
    l5_info["CheckL5Drivers"]["drivers"] = list(set(l5_info["CheckL5Drivers"]["drivers"]))
    return l5_info

def ReadAlternateFiles(ds, l4_info):
    ds_alt = {}
    l4ao = l4_info["GapFillFromAlternate"]["outputs"]
    # get a list of file names
    files = [l4ao[output]["file_name"] for output in list(l4ao.keys())]
    # read the alternate files
    for f in files:
        # if the file has not already been read, do it now
        if f not in ds_alt:
            ds_alternate = pfp_io.NetCDFRead(f, fixtimestepmethod="round")
            if ds_alternate.info["returncodes"]["value"] != 0: return ds_alt
            ds_alt[f] = gfalternate_matchstartendtimes(ds, ds_alternate)
    return ds_alt

def gfalternate_createdict(cf, ds, l4_info, label, called_by):
    """
    Purpose:
     Creates a dictionary in l4_info to hold information about the alternate data
     used to gap fill the tower data.
    Usage:
    Side effects:
    Author: PRI
    Date: August 2014
    """
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    # make the L4 "description" attrubute for the target variable
    descr_level = "description_" + ds.root["Attributes"]["processing_level"]
    # create the alternate data settings directory
    if called_by not in list(l4_info.keys()):
        # create the GapFillFromAlternate dictionary
        l4_info[called_by] = {"outputs": {}, "info": {}, "gui": {}}
        # only need to create the ["info"] dictionary on the first pass
        gfalternate_createdict_info(cf, ds, l4_info, called_by)
        if ds.info["returncodes"]["value"] != 0:
            return
        # only need to create the ["gui"] dictionary on the first pass
        gfalternate_createdict_gui(cf, ds, l4_info, called_by)
    # get the outputs section
    gfalternate_createdict_outputs(cf, l4_info, label, called_by)
    # create an empty series in ds if the alternate output series doesn't exist yet
    outputs = list(l4_info[called_by]["outputs"].keys())
    for output in outputs:
        if output not in list(ds.root["Variables"].keys()):
            l4_info["RemoveIntermediateSeries"]["not_output"].append(output)
            variable = pfp_utils.CreateEmptyVariable(output, nrecs)
            variable["Attr"][descr_level] = l4_info[called_by]["outputs"][output]["source"]
            pfp_utils.CreateVariable(ds, variable)
            variable = pfp_utils.CreateEmptyVariable(label + "_composite", nrecs)
            l4_info["RemoveIntermediateSeries"]["not_output"].append(label + "_composite")
            variable["Attr"][descr_level] = "Composite series of tower and alternate data"
            pfp_utils.CreateVariable(ds, variable)
    return

def gfalternate_createdict_info(cf, ds, l4_info, called_by):
    """
    Purpose:
     Create the "info" section of the l4_info[called_by] dictionary.
    Usage:
    Side effects:
    Author: PRI
    Date: Back in the day
          June 2019 - modified for new l4_info structure
    """
    # reset the return message and code
    ds.info["returncodes"]["message"] = "OK"
    ds.info["returncodes"]["value"] = 0
    # get a local pointer to the tower datetime series
    ldt = ds.root["Variables"]["DateTime"]["Data"]
    plot_path = pfp_utils.get_keyvaluefromcf(cf, ["Files"], "plot_path", default="./plots/")
    plot_path = os.path.join(plot_path, "L4", "")
    if not os.path.exists(plot_path):
        try:
            os.makedirs(plot_path)
        except OSError:
            msg = "Unable to create the plot path " + plot_path + "\n"
            msg = msg + "Press 'Quit' to edit the control file.\n"
            msg = msg + "Press 'Continue' to use the default path.\n"
            result = pfp_gui.MsgBox_ContinueOrQuit(msg, title="Warning: L4 plot path")
            if result.clickedButton().text() == "Quit":
                # user wants to edit the control file
                msg = " Quitting L4 to edit control file"
                logger.warning(msg)
                ds.info["returncodes"]["message"] = msg
                ds.info["returncodes"]["value"] = 1
            else:
                plot_path = "./plots/"
                cf["Files"]["plot_path"] = "./plots/"
    # check to see if this is a batch or an interactive run
    call_mode = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "call_mode", default="interactive")
    # create the alternate_info dictionary, this will hold much useful information
    l4a = l4_info[called_by]
    l4a["info"] = {"startdate": ldt[0].strftime("%Y-%m-%d %H:%M"),
                   "enddate": ldt[-1].strftime("%Y-%m-%d %H:%M"),
                   "called_by": called_by,
                   "plot_path": plot_path,
                   "call_mode": call_mode,
                   "time_step": int(float(ds.root["Attributes"]["time_step"])),
                   "site_name": ds.root["Attributes"]["site_name"]}
    out_filename = pfp_io.get_outfilenamefromcf(cf)
    xl_file_name = out_filename.replace('.nc', '_AlternateStats.xls')
    l4a["info"]["xl_file_name"] = xl_file_name
    return

def gfalternate_createdict_outputs(cf, l4_info, label, called_by):
    flag_codes = {"default": 400, "aws": 410, "access": 420, "erai": 430, "era5": 440}
    # name of alternate output series in ds
    outputs = list(cf["Drivers"][label][called_by].keys())
    # loop over the outputs listed in the control file
    l4ao = l4_info[called_by]["outputs"]
    cfalt = cf["Drivers"][label][called_by]
    for output in outputs:
        # disable output to netCDF file for this variable
        l4_info["RemoveIntermediateSeries"]["not_output"].append(output)
        # create the dictionary keys for this output
        l4ao[output] = {}
        # get the target
        sl = ["Drivers", label, called_by, output]
        l4ao[output]["target"] = pfp_utils.get_keyvaluefromcf(cf, sl, "target", default=label)
        # source name
        opt = pfp_utils.get_keyvaluefromcf(cf, sl, "source", default="")
        l4ao[output]["source"] = opt.lower()
        # output QC flag code
        l4ao[output]["flag_code"] = flag_codes["default"]
        if l4ao[output]["source"] in list(flag_codes.keys()):
            l4ao[output]["flag_code"] = flag_codes[l4ao[output]["source"]]
        # alternate data file name
        # first, look in the [Files] section for a generic file name
        file_list = list(cf["Files"].keys())
        lower_file_list = [item.lower() for item in file_list]
        if l4ao[output]["source"].lower() in lower_file_list:
            # found a generic file name
            i = lower_file_list.index(l4ao[output]["source"].lower())
            l4ao[output]["file_name"] = cf["Files"][file_list[i]]
        elif "file_name" in cfalt[output]:
            # no generic file name found, look for a file name in the variable section
            l4ao[output]["file_name"] = cfalt[output]["file_name"]
        else:
            # put up some warning messages
            msg = " Unable to resolve alternate source " + l4ao[output]["source"] + " for output " + output
            msg = msg + " of variable " + l4ao[output]["target"]
            logger.warning(msg)
            msg = " Output " + output + " for variable " + l4ao[output]["target"] + " will be skipped!"
            logger.warning(msg)
            # remove the entry in the l4ao dictionary
            l4ao.pop(output, None)
            continue
        # get the type of fit
        l4ao[output]["fit_type"] = "OLS"
        if "fit" in cfalt[output]:
            if cfalt[output]["fit"].lower() in ["ols", "ols_thru0", "mrev", "replace", "rma", "odr"]:
                l4ao[output]["fit_type"] = cfalt[output]["fit"]
            else:
                logger.info("gfAlternate: unrecognised fit option for series %s, used OLS", output)
        # correct for lag?
        if "lag" in cfalt[output]:
            if cfalt[output]["lag"].lower() in ["no", "false"]:
                l4ao[output]["lag"] = "no"
            elif cfalt[output]["lag"].lower() in ["yes", "true"]:
                l4ao[output]["lag"] = "yes"
            else:
                logger.info("gfAlternate: unrecognised lag option for series %s", output)
        else:
            l4ao[output]["lag"] = "yes"
        # choose specific alternate variable?
        if "usevars" in cfalt[output]:
            l4ao[output]["usevars"] = pfp_utils.string_to_list(cfalt[output]["usevars"])
        # alternate data variable name if different from name used in control file
        if "alternate_name" in cfalt[output]:
            l4ao[output]["alternate_name"] = cfalt[output]["alternate_name"]
        else:
            l4ao[output]["alternate_name"] = output
        # results of best fit for plotting later on
        l4ao[output]["results"] = {"startdate":[], "enddate":[], "No. points":[], "No. filled":[],
                                   "r":[], "Bias":[], "RMSE":[], "Frac Bias":[], "NMSE":[],
                                   "Avg (Tower)":[], "Avg (Alt)":[],
                                   "Var (Tower)":[], "Var (Alt)":[], "Var ratio":[]}

def gfalternate_createdict_gui(cf, ds, l4_info, called_by):
    """
    Purpose:
     Get settings for the GapFillFromAlternate routine from the [GUI]
     section of the L4 control file.
    Usage:
    Side effects:
    Author: PRI
    Date: July 2019
    """
    # local pointer to l4_info["GapFillFromAlternate"]
    l4a = l4_info[called_by]
    # local pointer to datetime
    ldt_tower = pfp_utils.GetVariable(ds, "DateTime")
    # populate the l4_info["gui"] dictionary with things that will be useful
    ts = int(float(ds.root["Attributes"]["time_step"]))
    l4a["gui"]["nperhr"] = int(float(60)/ts + 0.5)
    l4a["gui"]["nperday"] = int(float(24)*l4a["gui"]["nperhr"] + 0.5)
    l4a["gui"]["max_lags"] = int(float(12)*l4a["gui"]["nperhr"] + 0.5)
    # items from the [GUI] section of the control file
    sl = ["GUI", "GapFillFromAlternate"]
    # window length
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "period_option", default="manual")
    if opt == "manual":
        l4a["gui"]["period_option"] = 1
    elif opt == "monthly":
        l4a["gui"]["period_option"] = 2
        opt = pfp_utils.get_keyvaluefromcf(cf, sl, "number_months", default=3)
        l4a["gui"]["number_months"] = int(opt)
    elif opt == "days":
        l4a["gui"]["period_option"] = 3
        opt = pfp_utils.get_keyvaluefromcf(cf, sl, "number_days", default=90)
        l4a["gui"]["number_days"] = int(opt)
    # overwrite option
    l4a["gui"]["overwrite"] = False
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "overwrite", default="no")
    if opt.lower() == "yes": l4a["gui"]["overwrite"] = True
    # show plots option
    if cf["Options"]["call_mode"].lower() == "interactive":
        l4a["gui"]["show_plots"] = True
    else:
        opt = pfp_utils.get_keyvaluefromcf(cf, sl, "show_plots", default="no")
        if opt.lower() == "no":
            l4a["gui"]["show_plots"] = False
    # show all plots option
    l4a["gui"]["show_all"] = False
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "show_all", default="no")
    if opt.lower() == "yes": l4a["gui"]["show_all"] = True
    # auto-complete option
    l4a["gui"]["auto_complete"] = True
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "auto_complete", default="yes")
    if opt.lower() == "no": l4a["gui"]["auto_complete"] = False
    l4a["gui"]["autoforce"] = False
    # minimum percentage of good points required
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "min_percent", default=50)
    l4a["gui"]["min_percent"] = max(int(opt), 1)
    # start date of period to be gap filled
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "start_date", default="YYYY-MM-DD HH:mm")
    try:
        sd = dateutil.parser.parse(opt)
    except (ValueError, TypeError):
        sd = ldt_tower["Data"][0].strftime("%Y-%m-%d %H:%M")
    l4a["gui"]["startdate"] = sd
    # end date of period to be gap filled
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "end_date", default="YYYY-MM-DD HH:mm")
    try:
        ed = dateutil.parser.parse(opt)
    except (ValueError, TypeError):
        ed = ldt_tower["Data"][-1].strftime("%Y-%m-%d %H:%M")
    l4a["gui"]["enddate"] = ed
    return

def gfalternate_matchstartendtimes(ds, ds_alternate):
    """
    Purpose:
     Match the start and end times of the alternate and tower data.
     The logic is as follows:
      - if there is no overlap between the alternate and tower data then
        dummy series with missing data are created for the alternate data
        for the period of the tower data
      - if the alternate and tower data overlap then truncate or pad (with
        missing values) the alternate data series so that the periods of the
        tower data and alternate data match.
    Usage:
     gfalternate_matchstartendtimes(ds,ds_alternate)
     where ds is the data structure containing the tower data
           ds_alternate is the data structure containing the alternate data
    Author: PRI
    Date: July 2015
    Modifications:
     June 2022 - rewrote to use pfp_utils.GetVariable() and pfp_utils.CreateVariable()
                 and to return a new data structure instead of modify ds_alternate
                 in place.
    """
    # check the time steps are the same
    ts_tower = int(float(ds.root["Attributes"]["time_step"]))
    ts_alternate = int(float(ds_alternate.root["Attributes"]["time_step"]))
    if ts_tower != ts_alternate:
        msg = " GapFillFromAlternate: time step for tower and alternate data are different, returning ..."
        logger.error(msg)
        ds.info["returncodes"]["GapFillFromAlternate"] = "error"
        return
    # get a list of alternate series
    labels_alternate = [item for item in list(ds_alternate.root["Variables"].keys()) if "_QCFlag" not in item]
    for label in ["DateTime", "DateTime_UTC"]:
        if label in labels_alternate:
            labels_alternate.remove(label)
    # number of records in truncated or padded alternate data
    nRecs_tower = int(ds.root["Attributes"]["nc_nrecs"])
    # create new data strucure to hold alternate data spanning period of tower data
    gattrs = ds_alternate.root["Attributes"]
    ds_matched = pfp_io.DataStructure(global_attributes=gattrs)
    # force the matched datetime to be the tower datetime
    ds_matched.root["Variables"]["DateTime"] = copy.deepcopy(ds.root["Variables"]["DateTime"])
    # update the number of records in the file
    ds_matched.root["Attributes"]["nc_nrecs"] = nRecs_tower
    # get the start and end times of the tower and the alternate data and see if they overlap
    ldt_alternate = ds_alternate.root["Variables"]["DateTime"]["Data"]
    start_alternate = ldt_alternate[0]
    ldt_tower = ds.root["Variables"]["DateTime"]["Data"]
    end_tower = ldt_tower[-1]
    # since the datetime is monotonically increasing we need only check the start datetime
    overlap = start_alternate <= end_tower
    # do the alternate and tower data overlap?
    if overlap:
        # index of alternate datetimes that are also in tower datetimes
        tower_index, alternate_index = pfp_utils.FindMatchingIndices(ldt_tower, ldt_alternate)
        # check that the indices point to the same times
        ldta = [ldt_alternate[i] for i in alternate_index]
        ldtt = [ldt_tower[i] for i in tower_index]
        if ldta != ldtt:
            # and exit with a helpful message if they dont
            msg = " Something went badly wrong at L4 and I'm giving up"
            logger.error(msg)
            raise RuntimeError(msg)
        # loop over the alternate series and truncate or pad as required
        # truncation or padding is handled by the indices
        for label in labels_alternate:
            # get the alternate data
            var_alternate = pfp_utils.GetVariable(ds_alternate, label)
            # create an empty variable of the required length
            attr = copy.deepcopy(var_alternate["Attr"])
            var_overlap = pfp_utils.CreateEmptyVariable(label, nRecs_tower, attr=attr)
            # replace missing data with alternate data where times match
            var_overlap["Data"][tower_index] = var_alternate["Data"][alternate_index]
            var_overlap["Flag"][tower_index] = var_alternate["Flag"][alternate_index]
            # write the truncated or padded series back into the matched data structure
            pfp_utils.CreateVariable(ds_matched, var_overlap)
    else:
        # there is no overlap between the alternate and tower data, create dummy series
        for label in labels_alternate:
            # get the alternate data
            var_alternate = pfp_utils.GetVariable(ds_alternate, label)
            # create an empty variable of the required length
            attr = copy.deepcopy(var_alternate["Attr"])
            var_overlap = pfp_utils.CreateEmptyVariable(label, nRecs_tower, attr=attr)
            # write the truncated or padded series back into the matched data structure
            pfp_utils.CreateVariable(ds_matched, var_overlap)
    ds.info["returncodes"]["GapFillFromAlternate"] = "normal"
    return ds_matched

def gfClimatology_createdict(cf, ds, l4_info, label, called_by):
    """
    Purpose:
     Creates a dictionary in l4_info to hold information about the climatological data
     used to gap fill the tower data.
    Usage:
    Side effects:
    Author: PRI
    Date: August 2014
    """
    # make the L4 "description" attrubute for the target variable
    descr_level = "description_" + ds.root["Attributes"]["processing_level"]
    # flag codes
    flag_codes = {"interpolated daily": 450, "monthly": 460}
    # create the climatology directory in the data structure
    if called_by not in list(l4_info.keys()):
        l4_info[called_by] = {"outputs": {}}
    # name of alternate output series in ds
    outputs = list(cf["Drivers"][label][called_by].keys())
    # loop over the outputs listed in the control file
    l4co = l4_info[called_by]["outputs"]
    cfcli = cf["Drivers"][label][called_by]
    for output in outputs:
        # disable output to netCDF file for this variable
        l4_info["RemoveIntermediateSeries"]["not_output"].append(output)
        # create the dictionary keys for this output
        l4co[output] = {}
        # get the target
        sl = ["Drivers", label, called_by, output]
        l4co[output]["target"] = pfp_utils.get_keyvaluefromcf(cf, sl, "target", default=label)
        # get the source
        l4co[output]["source"] = pfp_utils.get_keyvaluefromcf(cf, sl, "source", default="climatology")
        # Climatology file name
        file_list = list(cf["Files"].keys())
        lower_file_list = [item.lower() for item in file_list]
        # first, look in the [Files] section for a generic file name
        if l4co[output]["source"] in lower_file_list:
            # found a generic file name
            i = lower_file_list.index(l4co[output]["source"].lower())
            l4co[output]["file_name"] = cf["Files"][file_list[i]]
        elif "file_name" in cfcli[output]:
            # no generic file name found, look for a file name in the variable section
            l4co[output]["file_name"] = cfcli[output]["file_name"]
        else:
            # put up some warning messages
            msg = " Unable to resolve climatology source " + l4co[output]["source"] + " for output " + output
            msg = msg + " of variable " + l4co[output]["target"]
            logger.warning(msg)
            msg = " Output " + output + " for variable " + l4co[output]["target"] + " will be skipped!"
            logger.warning(msg)
            # remove the entry in the l4ao dictionary
            l4co.pop(output, None)
            continue
        # climatology variable name if different from name used in control file
        opt = pfp_utils.get_keyvaluefromcf(cfcli, [output], "climatology_name", default=label)
        l4co[output]["climatology_name"] = str(opt)
        # climatology gap filling method
        opt = pfp_utils.get_keyvaluefromcf(cfcli, [output], "method", default="interpolated daily")
        if opt in ["interpolated daily", "monthly"]:
            l4co[output]["method"] = opt
        else:
            l4co[output]["method"] = "interpolated daily"
        # get the flag code
        l4co[output]["flag_code"] = flag_codes[l4co[output]["method"]]
        # create an empty series in ds if the climatology output series doesn't exist yet
        if output not in list(ds.root["Variables"].keys()):
            nrecs = int(float(ds.root["Attributes"]["nc_nrecs"]))
            var = pfp_utils.CreateEmptyVariable(output, nrecs)
            var["Attr"][descr_level] = "climatology"
            pfp_utils.CreateVariable(ds, var)
    return

def gfMDS_createdict(cf, ds, l5_info, label, called_by, flag_code):
    """
    Purpose:
     Create an information dictionary for MDS gap filling from the contents
     of the control file.
    Usage:
     gfMDS_createdict(cf, ds, l5_info, label, called_by)
    Author: PRI
    Date: May 2018
    """
    # create the solo settings directory
    if called_by not in l5_info:
        l5_info[called_by] = {"outputs": {}, "info": {}}
    # file path and input file name
    l5_info[called_by]["info"]["file_path"] = cf["Files"]["file_path"]
    l5_info[called_by]["info"]["in_filename"] = cf["Files"]["in_filename"]
    # get the executable suffix
    l5_info[called_by]["info"]["executable_suffix"] = pfp_utils.get_executable_suffix()
    # get the plot path
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Files"], "plot_path", default="./plots/")
    plot_path = os.path.join(opt, "L5", "")
    if not os.path.exists(plot_path):
        try:
            os.makedirs(plot_path)
        except OSError:
            msg = "Unable to create the plot path " + plot_path + "\n"
            msg = msg + "Press 'Quit' to edit the control file.\n"
            msg = msg + "Press 'Continue' to use the default path.\n"
            result = pfp_gui.MsgBox_ContinueOrQuit(msg, title="Warning: L5 plot path")
            if result.clickedButton().text() == "Quit":
                # user wants to edit the control file
                msg = " Quitting L5 to edit control file"
                logger.warning(msg)
                ds.info["returncodes"]["message"] = msg
                ds.info["returncodes"]["value"] = 1
            else:
                plot_path = "./plots/"
                cf["Files"]["plot_path"] = "./plots/"
    l5_info[called_by]["info"]["plot_path"] = plot_path
    # get the maximum length for "short" gaps in days
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "MaxShortGapDays", default=14)
    max_short_gap_days = int(opt)
    l5_info[called_by]["info"]["MaxShortGapDays"] = max_short_gap_days
    # maximum length in records
    ts = int(float(ds.root["Attributes"]["time_step"]))
    nperday = 24 * 60//ts
    l5_info[called_by]["info"]["MaxShortGapRecords"] = max_short_gap_days * nperday
    # name of MDS output series in ds
    outputs = list(cf["Fluxes"][label]["GapFillUsingMDS"].keys())
    # loop over the outputs listed in the control file
    l5mo = l5_info[called_by]["outputs"]
    for output in outputs:
        # disable output to netCDF file for this variable
        if output not in list(cf["Fluxes"].keys()):
            l5_info["RemoveIntermediateSeries"]["not_output"].append(output)
        # create the dictionary keys for this series
        l5mo[output] = {}
        # get the target
        if "target" in cf["Fluxes"][label]["GapFillUsingMDS"][output]:
            l5mo[output]["target"] = cf["Fluxes"][label]["GapFillUsingMDS"][output]["target"]
        else:
            l5mo[output]["target"] = label
        # list of MDS settings
        if "mds_settings" in cf["Fluxes"][label]["GapFillUsingMDS"][output]:
            mdss_string = cf["Fluxes"][label]["GapFillUsingMDS"][output]["mds_settings"]
            mdss_string = mdss_string.replace(" ","")
            if "," in mdss_string:
                l5mo[output]["mds_settings"] = mdss_string.split(",")
            else:
                l5mo[output]["mds_settings"] = [mdss_string]
        # list of drivers
        drivers_string = cf["Fluxes"][label]["GapFillUsingMDS"][output]["drivers"]
        drivers_string = drivers_string.replace(" ","")
        if "," in drivers_string:
            if len(drivers_string.split(",")) == 3:
                l5mo[output]["drivers"] = drivers_string.split(",")
                # append to list of drivers to be checked for gaps
                l5_info["CheckL5Drivers"]["drivers"] += drivers_string.split(",")
            else:
                msg = " MDS: incorrect number of drivers for " + label + ", skipping ..."
                logger.error(msg)
                continue
        else:
            msg = " MDS: incorrect number of drivers for " + label + ", skipping ..."
            logger.error(msg)
            continue
        # list of tolerances
        tolerances_string = cf["Fluxes"][label]["GapFillUsingMDS"][output]["tolerances"]
        tolerances_string = tolerances_string.replace(" ","")
        tolerances_string = tolerances_string.replace("(","").replace(")","")
        if "," in tolerances_string:
            if len(tolerances_string.split(",")) == 4:
                parts = tolerances_string.split(",")
                l5mo[output]["tolerances"] = [(parts[0], parts[1]), parts[2], parts[3]]
            else:
                msg = " MDS: incorrect format for tolerances for " + label + ", skipping ..."
                logger.error(msg)
                continue
        else:
            msg = " MDS: incorrect format for tolerances for " + label + ", skipping ..."
            logger.error(msg)
            continue
    # check that all requested targets and drivers have a mapping to
    # a FluxNet label, remove if they don't
    fluxnet_label_map = {"Fco2":"NEE", "Fe":"LE", "Fh":"H",
                         "Fsd":"SW_IN", "Ta":"TA", "VPD":"VPD"}
    for mds_label in list(l5_info[called_by]["outputs"]):
        l5_info[called_by]["outputs"][mds_label]["mds_label"] = mds_label
        pfp_target = l5_info[called_by]["outputs"][mds_label]["target"]
        if pfp_target not in fluxnet_label_map:
            msg = " Target ("+pfp_target+") not supported for MDS gap filling"
            logger.warning(msg)
            del l5_info[called_by]["outputs"][mds_label]
            continue
        else:
            l5_info[called_by]["outputs"][mds_label]["target_mds"] = fluxnet_label_map[pfp_target]
        pfp_drivers = l5_info[called_by]["outputs"][mds_label]["drivers"]
        for pfp_driver in pfp_drivers:
            if pfp_driver not in fluxnet_label_map:
                msg = "Driver ("+pfp_driver+") not supported for MDS gap filling"
                logger.warning(msg)
                l5_info[called_by]["outputs"][mds_label]["drivers"].remove(pfp_driver)
            else:
                if "drivers_mds" not in l5_info[called_by]["outputs"][mds_label]:
                    l5_info[called_by]["outputs"][mds_label]["drivers_mds"] = []
                l5_info[called_by]["outputs"][mds_label]["drivers_mds"].append(fluxnet_label_map[pfp_driver])
        if len(l5_info[called_by]["outputs"][mds_label]["drivers"]) == 0:
            del l5_info[called_by]["outputs"][mds_label]
            continue
    return

def gfMergeSeries_createdict(cf, ds, info, label, called_by):
    """ Creates a dictionary in ds to hold information about the merging of gap filled
        and tower data."""
    merge_prereq_list = ["Fsd","Fsu","Fld","Flu","Ts","Sws"]
    # get the section of the control file containing the series
    section = pfp_utils.get_cfsection(cf, label, mode="quiet")
    # create the merge directory in the data structure
    if called_by not in info:
        info[called_by] = {}
    # check to see if this series is in the "merge first" list
    # series in the "merge first" list get merged first so they can be used with existing tower
    # data to re-calculate Fg, Fn and Fa
    merge_order = "standard"
    if label in merge_prereq_list:
        merge_order = "prerequisite"
    if merge_order not in list(info[called_by].keys()):
        info[called_by][merge_order] = {}
    # create the dictionary keys for this series
    info[called_by][merge_order][label] = {}
    # output series name
    info[called_by][merge_order][label]["output"] = label
    # merge source list
    section = pfp_utils.get_cfsection(cf, label, mode="quiet")
    src_list = pfp_utils.GetMergeSeriesKeys(cf, label, section=section)
    info[called_by][merge_order][label]["source"] = src_list
    # create an empty series in ds if the output series doesn't exist yet
    if label not in list(ds.root["Variables"].keys()):
        nrecs = int(float(ds.root["Attributes"]["nc_nrecs"]))
        var = pfp_utils.CreateEmptyVariable(label, nrecs)
        pfp_utils.CreateVariable(ds, var)
    return

def gfSOLO_createdict(cf, ds, l5_info, target_label, called_by, flag_code):
    """
    Purpose:
     Creates a dictionary in l5_info to hold information about the SOLO data
     used to gap fill the tower data.
    Usage:
    Side effects:
    Author: PRI
    Date: August 2014
    """
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    # make the L5 "description" attrubute for the target variable
    descr_level = "description_" + ds.root["Attributes"]["processing_level"]
    ds.root["Variables"][target_label]["Attr"][descr_level] = ""
    # create the solo settings directory
    if called_by not in list(l5_info.keys()):
        # create the GapFillUsingSOLO dictionary
        l5_info[called_by] = {"outputs": {}, "info": {}, "gui": {}}
        # only need to create the ["info"] dictionary on the first pass
        gfSOLO_createdict_info(cf, ds, l5_info, called_by)
        if ds.info["returncodes"]["value"] != 0:
            return
        # only need to create the ["gui"] dictionary on the first pass
        gfSOLO_createdict_gui(cf, ds, l5_info, called_by)
    # get the outputs section
    gfSOLO_createdict_outputs(cf, l5_info, target_label, called_by, flag_code)
    # add the summary plors section
    if "SummaryPlots" in cf:
        l5_info[called_by]["SummaryPlots"] = cf["SummaryPlots"]
    # create an empty series in ds if the SOLO output series doesn't exist yet
    outputs = list(cf["Fluxes"][target_label][called_by].keys())
    for output in outputs:
        l5_info["CheckL5Drivers"]["drivers"] += l5_info[called_by]["outputs"][output]["drivers"]
        if output not in list(ds.root["Variables"].keys()):
            # disable output to netCDF file for this variable
            l5_info["RemoveIntermediateSeries"]["not_output"].append(output)
            # create an empty variable
            variable = pfp_utils.CreateEmptyVariable(output, nrecs)
            # update the empty variable attributes
            target_variable = pfp_utils.GetVariable(ds, target_label)
            for vattr in ["long_name", "valid_range", "units", "standard_name"]:
                if vattr in target_variable["Attr"]:
                    variable["Attr"][vattr] = target_variable["Attr"][vattr]
            # create L5 specific vards.root["Variables"][label]iable attributes
            variable["Attr"]["drivers"] = l5_info[called_by]["outputs"][output]["drivers"]
            variable["Attr"]["target"] = target_label
            pfp_utils.append_to_attribute(variable["Attr"], {descr_level: "SOLO"})
            pfp_utils.CreateVariable(ds, variable)
    return

def gfSOLO_createdict_gui(cf, ds, l5_info, called_by):
    """
    Purpose:
     Get settings for the GapFillUsingSOLO routine from the [GUI]
     section of the L5 control file.
    Usage:
    Side effects:
    Author: PRI
    Date: August 2019
    """
    # local pointer to l5_info[called_by]
    l5s = l5_info[called_by]
    # local pointer to datetime
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    # populate the l5_info["gui"] dictionary with things that will be useful
    ts = int(float(ds.root["Attributes"]["time_step"]))
    # number of records per hour and per day
    l5s["gui"]["nperhr"] = int(float(60)/ts + 0.5)
    l5s["gui"]["nperday"] = int(float(24)*l5s["gui"]["nperhr"] + 0.5)
    # items from the [GUI] section of the control file
    sl = ["GUI", called_by]
    # window length
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "period_option", default="manual")
    if opt == "manual":
        l5s["gui"]["period_option"] = 1
    elif opt == "monthly":
        l5s["gui"]["period_option"] = 2
        opt = pfp_utils.get_keyvaluefromcf(cf, sl, "number_months", default=2)
        l5s["gui"]["number_months"] = int(opt)
    elif opt == "days":
        l5s["gui"]["period_option"] = 3
        opt = pfp_utils.get_keyvaluefromcf(cf, sl, "number_days", default=60)
        l5s["gui"]["number_days"] = int(opt)
    # overwrite option
    l5s["gui"]["overwrite"] = False
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "overwrite", default="no")
    if opt.lower() == "yes": l5s["gui"]["overwrite"] = True
    # show plots option
    if cf["Options"]["call_mode"].lower() == "interactive":
        l5s["gui"]["show_plots"] = True
    else:
        opt = pfp_utils.get_keyvaluefromcf(cf, sl, "show_plots", default="no")
        if opt.lower() == "no":
            l5s["gui"]["show_plots"] = False
    # show all plots option
    l5s["gui"]["show_all"] = False
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "show_all", default="no")
    if opt.lower() == "yes": l5s["gui"]["show_all"] = True
    # auto-complete option
    l5s["gui"]["auto_complete"] = True
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "auto_complete", default="yes")
    if opt.lower() == "no": l5s["gui"]["auto_complete"] = False
    # minimum percentage of good points required
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "min_percent", default=50)
    l5s["gui"]["min_percent"] = max(int(opt), 1)
    # nodes for SOFM/SOLO network
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "nodes", default="auto")
    l5s["gui"]["nodes"] = str(opt)
    # training iterations
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "training", default="500")
    l5s["gui"]["training"] = str(opt)
    # nda factor
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "nda_factor", default="5")
    l5s["gui"]["nda_factor"] = str(opt)
    # learning rate
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "learning_rate", default="0.001")
    l5s["gui"]["learning_rate"] = str(opt)
    # learning iterations
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "iterations", default="500")
    l5s["gui"]["iterations"] = str(opt)
    # start date of period to be gap filled
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "start_date", default="YYYY-MM-DD HH:mm")
    try:
        sd = dateutil.parser.parse(opt)
    except (ValueError, TypeError):
        sd = ldt["Data"][0].strftime("%Y-%m-%d %H:%M")
    l5s["gui"]["startdate"] = sd
    # end date of period to be gap filled
    opt = pfp_utils.get_keyvaluefromcf(cf, sl, "end_date", default="YYYY-MM-DD HH:mm")
    try:
        ed = dateutil.parser.parse(opt)
    except (ValueError, TypeError):
        ed = ldt["Data"][-1].strftime("%Y-%m-%d %H:%M")
    l5s["gui"]["enddate"] = ed
    return

def gfSOLO_createdict_info(cf, ds, l5_info, called_by):
    """
    Purpose:
     Create the "info" section of the l5_info[called_by] dictionary.
    Usage:
    Side effects:
    Author: PRI
    Date: Back in the day
          June 2019 - modified for new l5_info structure
    """
    l5s = l5_info[called_by]
    # reset the return message and code
    ds.info["returncodes"]["message"] = "OK"
    ds.info["returncodes"]["value"] = 0
    # time step
    time_step = int(float(ds.root["Attributes"]["time_step"]))
    # get the level of processing
    level = ds.root["Attributes"]["processing_level"]
    # local pointer to the datetime series
    ldt = ds.root["Variables"]["DateTime"]["Data"]
    # add an info section to the l5_info[called_by] dictionary
    #l5s["info"]["file_startdate"] = ldt[0].strftime("%Y-%m-%d %H:%M")
    #l5s["info"]["file_enddate"] = ldt[-1].strftime("%Y-%m-%d %H:%M")
    l5s["info"]["startdate"] = ldt[0].strftime("%Y-%m-%d %H:%M")
    l5s["info"]["enddate"] = ldt[-1].strftime("%Y-%m-%d %H:%M")
    l5s["info"]["called_by"] = called_by
    l5s["info"]["time_step"] = time_step
    l5s["info"]["site_name"] = ds.root["Attributes"]["site_name"]
    l5s["info"]["executable_suffix"] = pfp_utils.get_executable_suffix()
    # check to see if this is a batch or an interactive run
    call_mode = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "call_mode", default="interactive")
    l5s["info"]["call_mode"] = call_mode
    # truncate to last date in Imports?
    truncate = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "TruncateToImports", default="Yes")
    l5s["info"]["truncate_to_imports"] = truncate
    # number of records per day and maximum lags
    nperhr = int(float(60)/time_step + 0.5)
    l5s["info"]["nperday"] = int(float(24)*nperhr + 0.5)
    l5s["info"]["maxlags"] = int(float(12)*nperhr + 0.5)
    # get the plot path
    plot_path = pfp_utils.get_keyvaluefromcf(cf, ["Files"], "plot_path", default="./plots/")
    plot_path = os.path.join(plot_path, level, "")
    if not os.path.exists(plot_path):
        try:
            os.makedirs(plot_path)
        except OSError:
            msg = "Unable to create the plot path " + plot_path + "\n"
            msg = msg + "Press 'Quit' to edit the control file.\n"
            msg = msg + "Press 'Continue' to use the default path.\n"
            title = "Warning: " + level + " plot path"
            result = pfp_gui.MsgBox_ContinueOrQuit(msg, title=title)
            if result.clickedButton().text() == "Quit":
                # user wants to edit the control file
                msg = " Quitting " + level + " to edit control file"
                logger.warning(msg)
                ds.info["returncodes"]["message"] = msg
                ds.info["returncodes"]["value"] = 1
            else:
                plot_path = "./plots/"
                cf["Files"]["plot_path"] = "./plots/"
    l5s["info"]["plot_path"] = plot_path
    return

def gfSOLO_createdict_outputs(cf, l5_info, target, called_by, flag_code):
    """
    Purpose:
     Create the l5_info[called_by]["outputs"] dictionary.
    Usage:
    Side effects:
    Author: PRI
    Date: Back in the day
        August 2019 - rewritten for consistency and batch operation
    """
    level = cf["level"]
    iris = l5_info["RemoveIntermediateSeries"]
    so = l5_info[called_by]["outputs"]
    # loop over the outputs listed in the control file
    if level == "L5":
        section = "Fluxes"
        drivers = "Fn,Fg,SHD,SH,Ta,Ts"
        source = target
    elif level == "L6":
        section = "EcosystemRespiration"
        drivers = "Ta,Ts,Sws"
        #source = "Fco2"
        source = target
    else:
        msg = "Unrecognised control file level (must be L5 or L6)"
        logger.error(msg)
        return
    outputs = list(cf[section][target][called_by].keys())
    for output in outputs:
        # disable output to netCDF file for this variable
        iris["not_output"].append(output)
        # create the dictionary keys for this series
        so[output] = {}
        # get the target
        sl = [section, target, called_by, output]
        so[output]["target"] = pfp_utils.get_keyvaluefromcf(cf, sl, "target", default=target)
        so[output]["source"] = pfp_utils.get_keyvaluefromcf(cf, sl, "source", default=source)
        # add the flag_code
        so[output]["flag_code"] = flag_code
        # list of SOLO settings
        if "solo_settings" in cf[section][target][called_by][output]:
            src_string = cf[section][target][called_by][output]["solo_settings"]
            src_list = src_string.split(",")
            so[output]["solo_settings"] = {}
            so[output]["solo_settings"]["nodes_target"] = int(src_list[0])
            so[output]["solo_settings"]["training"] = int(src_list[1])
            so[output]["solo_settings"]["nda_factor"] = int(src_list[2])
            so[output]["solo_settings"]["learning_rate"] = float(src_list[3])
            so[output]["solo_settings"]["iterations"] = int(src_list[4])
        # list of drivers
        opt = pfp_utils.get_keyvaluefromcf(cf, sl, "drivers", default=drivers)
        so[output]["drivers"] = pfp_utils.string_to_list(opt)
        # fit statistics for plotting later on
        so[output]["results"] = {"startdate":[],"enddate":[],"No. points":[],"r":[],
                                 "Bias":[],"RMSE":[],"Frac Bias":[],"NMSE":[],
                                 "Avg (obs)":[],"Avg (SOLO)":[],
                                 "Var (obs)":[],"Var (SOLO)":[],"Var ratio":[],
                                 "m_ols":[],"b_ols":[]}
    return

# functions for GapFillFromClimatology
def GapFillFromClimatology(ds, l4_info, called_by):
    '''
    Gap fill missing data using data from the climatology spreadsheet produced by
    the climatology.py script.
    '''
    if called_by not in list(l4_info.keys()):
        return
    l4co = l4_info[called_by]["outputs"]
    # tell the user what we are going to do
    msg = " Reading climatology file and creating climatology series"
    logger.info(msg)
    # loop over the series to be gap filled using climatology
    cli_xlbooks = {}
    for output in list(l4co.keys()):
        # check to see if there are any gaps in "series"
        #index = numpy.where(abs(ds.root["Variables"][label]['Data']-float(c.missing_value))<c.eps)[0]
        #if len(index)==0: continue                      # no gaps found in "series"
        cli_filename = l4co[output]["file_name"]
        if not os.path.exists(cli_filename):
            logger.error(" GapFillFromClimatology: Climatology file %s doesn't exist", cli_filename)
            continue
        if cli_filename not in cli_xlbooks:
            cli_xlbooks[cli_filename] = xlrd.open_workbook(cli_filename)
        # local pointers to the series name and climatology method
        label = l4co[output]["target"]
        method = l4co[output]["method"]
        flag_code = l4co[output]["flag_code"]
        # do the gap filling
        cli_xlbook = cli_xlbooks[cli_filename]
        # choose the gap filling method
        if method == "interpolated daily":
            gfClimatology_interpolateddaily(ds, label, output, cli_xlbook, flag_code)
        else:
            logger.error(" GapFillFromClimatology: unrecognised method option for %s", label)
            continue

def gfClimatology_interpolateddaily(ds, series, output, xlbook, flag_code):
    """
    Gap fill using data interpolated over a 2D array where the days are
    the rows and the time of day is the columns.
    """
    # description string for this level of processing
    descr_level = "description_" + ds.root["Attributes"]["processing_level"]
    # gap fill from interpolated 30 minute data
    sheet_name = series + 'i(day)'
    if sheet_name not in xlbook.sheet_names():
        msg = " gfClimatology: sheet " + sheet_name + " not found, skipping ..."
        logger.warning(msg)
        return
    ldt = ds.root["Variables"]["DateTime"]["Data"]
    thissheet = xlbook.sheet_by_name(sheet_name)
    datemode = xlbook.datemode
    basedate = datetime.datetime(1899, 12, 30)
    nts = thissheet.ncols - 1
    ndays = thissheet.nrows - 2
    # read the time stamp values from the climatology worksheet
    tsteps = thissheet.row_values(1, start_colx=1, end_colx=nts+1)
    # read the data from the climatology workbook
    val1d = numpy.ma.zeros(ndays*nts, dtype=numpy.float64)
    # initialise an array for the datetime of the climatological values
    cdt = [None]*nts*ndays
    # loop over the rows (days) of data
    for xlRow in range(ndays):
        # get the Excel datetime value
        xldatenumber = int(thissheet.cell_value(xlRow+2, 0))
        # convert this to a Python Datetime
        xldatetime = basedate + datetime.timedelta(days=xldatenumber + 1462*datemode)
        # fill the climatology datetime array
        cdt[xlRow*nts:(xlRow+1)*nts] = [xldatetime+datetime.timedelta(hours=hh) for hh in tsteps]
        # fill the climatological value array
        val1d[xlRow*nts:(xlRow+1)*nts] = thissheet.row_values(xlRow+2, start_colx=1, end_colx=nts+1)
    # get the data to be filled with climatological values
    var = pfp_utils.GetVariable(ds, series)
    # get an index of missing values
    idx = numpy.where(numpy.ma.getmaskarray(var["Data"]) == True)[0]
    # there must be a better way to do this ...
    # simply using the index (idx) to set a slice of the data array to the gap filled values in val1d
    # does not seem to work (mask stays true on replaced values in data), the work around is to
    # step through the indices, find the time of the missing value in data, find the same time in the
    # gap filled values val1d and set the missing element of data to this element of val1d
    # actually ...
    # this may not be the fastest but it may be the most robust because it matches dates of missing data
    # to dates in the climatology file
    for ii in idx:
        try:
            jj = pfp_utils.find_nearest_value(cdt, ldt[ii])
            var["Data"][ii] = val1d[jj]
            var["Flag"][ii] = numpy.int32(flag_code)
        except ValueError:
            var["Data"][ii] = numpy.float64(c.missing_value)
            var["Flag"][ii] = numpy.int32(flag_code+1)
    # put the gap filled data back into the data structure
    var["Label"] = output
    var["Attr"][descr_level] = "climatology"
    pfp_utils.CreateVariable(ds, var)
    return

def gfClimatology_monthly(ds, series, output, xlbook):
    """ Gap fill using monthly climatology."""
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    Hdh = numpy.array([(d.hour + d.minute/float(60)) for d in ldt["Data"]])
    Month = numpy.array([d.month for d in ldt["Data"]])
    thissheet = xlbook.sheet_by_name(series)
    val1d = numpy.zeros_like(ds.root["Variables"][series]["Data"])
    values = numpy.zeros([48, 12])
    for month in range(1, 13):
        m = (month - 1)
        xlCol = m*5 + 2
        values[:, m] = thissheet.col_values(xlCol)[2:50]
    for i in range(len(ds.root["Variables"][series]["Data"])):
        h = int(2*Hdh[i])
        m = int(Month[i])
        val1d[i] = values[h, m-1]
    index = numpy.where(abs(ds.root["Variables"][output]["Data"] - c.missing_value) < c.eps)[0]
    ds.root["Variables"][output]["Data"][index] = val1d[index]
    ds.root["Variables"][output]["Flag"][index] = numpy.int32(460)

# functions for GapFillUsingInterpolation
def GapFillUsingInterpolation(ds, info):
    """
    Purpose:
     Gap fill variables in the data structure using interpolation.
     All variables in the [Variables], [Drivers] and [Fluxes] section
     are processed.
    Usage:
     pfp_gf.GapFillUsingInterpolation(cf,ds)
     where cf is a control file object
           ds is a data structure
    Author: PRI
    Date: September 2016
    """
    # get the maximum gap length to be filled by interpolation
    max_length_hours = info["GapFillUsingInterpolation"]["MaxGapInterpolate"]
    # bug out if interpolation disabled in control file
    if max_length_hours == 0:
        msg = " Gap fill by interpolation disabled in control file"
        logger.info(msg)
        return
    # get list of variables from control file
    targets = info["GapFillUsingInterpolation"]["targets"]
    # get the interpolation type
    int_type = info["GapFillUsingInterpolation"]["InterpolateType"]
    # tell the user what we are doing
    msg = " Using " + int_type +" interpolation (max. gap = " + str(max_length_hours) +" hours)"
    logger.info(msg)
    # do the business
    pfp_ts.InterpolateOverMissing(ds, targets, max_length_hours=max_length_hours, int_type=int_type)
    return

# miscellaneous L4 routines
def gf_getdiurnalstats(DecHour,Data,ts):
    nInts = 24*int((60/ts)+0.5)
    Num = numpy.ma.zeros(nInts,dtype=int)
    Hr = numpy.ma.zeros(nInts,dtype=float)
    for i in range(nInts):
        Hr[i] = float(i)*ts/60.
    Av = numpy.ma.masked_all(nInts)
    Sd = numpy.ma.masked_all(nInts)
    Mx = numpy.ma.masked_all(nInts)
    Mn = numpy.ma.masked_all(nInts)
    if numpy.size(Data)!=0:
        for i in range(nInts):
            li = numpy.ma.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-float(c.missing_value))>c.eps))
            Num[i] = numpy.size(li)
            if Num[i]!=0:
                Av[i] = numpy.ma.mean(Data[li])
                Sd[i] = numpy.ma.std(Data[li])
                Mx[i] = numpy.ma.max(Data[li])
                Mn[i] = numpy.ma.min(Data[li])
    return Num, Hr, Av, Sd, Mx, Mn

def ImportSeries(ds, info):
    cfg = info["cfg"]
    # check to see if there is an Imports section
    if "Imports" not in list(cfg.keys()):
        return
    info["ImportSeries"] = {}
    # number of records
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    # get the start and end datetime
    ldt = ds.root["Variables"]["DateTime"]["Data"]
    start_date = ldt[0]
    end_date = ldt[-1]
    # loop over the series in the Imports section
    for label in list(cfg["Imports"].keys()):
        import_filename = pfp_utils.get_keyvaluefromcf(cfg, ["Imports", label],
                                                       "file_name", default="")
        if import_filename == "":
            msg = " ImportSeries: import filename not found in control file, skipping ..."
            logger.warning(msg)
            continue
        var_name = pfp_utils.get_keyvaluefromcf(cfg, ["Imports", label],
                                                "var_name", default="")
        if var_name == "":
            msg = " ImportSeries: variable name not found in control file, skipping ..."
            logger.warning(msg)
            continue
        ds_import = pfp_io.NetCDFRead(import_filename)
        if ds_import.info["returncodes"]["value"] != 0:
            return
        if var_name not in list(ds_import.root["Variables"].keys()):
            msg = " Requested variable not found in imported data"
            logger.warning(msg)
            continue
        msg = "  Importing variable " + label
        logger.info(msg)
        ts_import = int(float(ds_import.root["Attributes"]["time_step"]))
        ldt_import = ds_import.root["Variables"]["DateTime"]["Data"]
        si = pfp_utils.GetDateIndex(ldt_import, start_date, ts=ts_import,
                                    default=0, match="exact")
        ei = pfp_utils.GetDateIndex(ldt_import, end_date, ts=ts_import,
                                    default=len(ldt_import)-1, match="exact")
        var_import = pfp_utils.GetVariable(ds_import, var_name, start=si, end=ei)
        var_import["Attr"]["time_coverage_start"] = ldt_import[0].strftime("%Y-%m-%d %H:%M")
        var_import["Attr"]["time_coverage_end"] = ldt_import[-1].strftime("%Y-%m-%d %H:%M")
        ldt_import = ldt_import[si:ei+1]
        indainb, indbina = pfp_utils.FindMatchingIndices(ldt_import, ldt)
        var = pfp_utils.CreateEmptyVariable(label, nrecs, attr=var_import["Attr"])
        var["Data"][indbina] = var_import["Data"][indainb]
        var["Flag"][indbina] = var_import["Flag"][indainb]
        pfp_utils.CreateVariable(ds, var)
        info["ImportSeries"][label] = {"start": ldt[indbina[0]],
                                       "end": ldt[indbina[-1]]}
    return
