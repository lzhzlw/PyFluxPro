# standard modules
import logging
# 3rd party
import numpy
# PFP modules
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def Ws_from_Ux_Uy(ds, Ws_out, Ux_in, Uy_in):
    """
    Purpose:
     Function to calculate wind speed from the horizontal components.
    Usage:
     pfp_func_transforms.Ws_from_Ux_Uy(Ux_in, Uy_in)
    Author: PRI
    Date: August 2023
    """
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    Ux = pfp_utils.GetVariable(ds, Ux_in)
    Uy = pfp_utils.GetVariable(ds, Uy_in)
    attr = {"long_name": "Wind speed", "units": "m/s",
            "statistic_type": "average"}
    Ws = pfp_utils.CreateEmptyVariable(Ws_out, nrecs, attr=attr)
    Ws["Data"] = numpy.ma.sqrt(Ux["Data"]*Ux["Data"] + Uy["Data"]*Uy["Data"])
    Ws["Flag"] = numpy.ma.where(numpy.ma.getmaskarray(Ws["Data"]), ones, zeros)
    pfp_utils.CreateVariable(ds, Ws)
    return 1
def Wd_from_Ux_Uy(ds, Wd_out, Ux_in, Uy_in):
    """
    Purpose:
     Function to calculate wind direction from the horizontal components.
    Usage:
     pfp_func_transforms.Wd_from_Ux_Uy(Ux_in, Uy_in)
    Author: PRI
    Date: August 2023
    """
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    Ux = pfp_utils.GetVariable(ds, Ux_in)
    Uy = pfp_utils.GetVariable(ds, Uy_in)
    attr = {"long_name": "Wind direction", "units": "degrees",
            "statistic_type": "average"}
    Wd = pfp_utils.CreateEmptyVariable(Wd_out, nrecs, attr=attr)
    if ((Ux["Attr"]["instrument"] in ["WindMaster Pro"]) and
        (Uy["Attr"]["instrument"] == Ux["Attr"]["instrument"])):
        Wd_sonic = numpy.ma.mod(360 - numpy.degrees(numpy.ma.arctan2(Uy["Data"], Ux["Data"])), 360)
        Wd["Attr"]["instrument"] = Ux["Attr"]["instrument"]
    elif ((Ux["Attr"]["instrument"] in ["CSAT", "CSAT3A", "CSAT3B"]) and
          (Uy["Attr"]["instrument"] == Ux["Attr"]["instrument"])):
        Wd_sonic = numpy.ma.mod(180 - numpy.degrees(numpy.ma.arctan2(Uy["Data"], Ux["Data"])), 360)
        Wd["Attr"]["instrument"] = Ux["Attr"]["instrument"]
    else:
        msg = " Unrecognised sonic anemometer type (" + Ux["instrument"]
        msg += Uy["instrument"] + ")"
        logger.warning(msg)
    # we want the direction FROM which the wind blows
    Wd["Data"] = numpy.ma.mod(Wd_sonic + 180, 360)
    Wd["Flag"] = numpy.ma.where(numpy.ma.getmaskarray(Wd["Data"]), ones, zeros)
    pfp_utils.CreateVariable(ds, Wd)
    return 1
