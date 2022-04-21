"""
Prepare boundary condition files for ANUGA simulation
Grabs relevant USGS/NOAA data based on information in settings.py
"""
import pandas as pd
import os
from .settings import *
import noaa_coops as nc
from nwis import Nwis
nwis_data = Nwis()

# DISCHARGE-------------------------------------------------------------------
# Retrieve data
WLO = nwis_data.get_data(site=ID_Calumet, start_date=CDT_dates[0], 
                         end_date=CDT_dates[1], data_type='iv').to_dataframe()
ATC = nwis_data.get_data(site=ID_MorganCity, start_date=CDT_dates[0], 
                          end_date=CDT_dates[1], data_type='iv').to_dataframe()
FRA = nwis_data.get_data(site=ID_Franklin, start_date=CDT_dates[0], 
                         end_date=CDT_dates[1], data_type='iv').to_dataframe()

# Convert to metric
WLO['Q'] = WLO['00060']*0.028316847 # CFS to m3/s
ATC['Q'] = ATC['00060']*0.028316847 # CFS to m3/s
try:
    # Franklin only has real Q data up to 2019
    FRA['Q'] = -1*FRA['00060']*0.028316847 # CFS to m3/s
except:
    # If unavailable, compute Q from WL rating curve
    FRA['WL'] = (FRA['00065']-1.21)*0.3048 # ft to m navd88
    FRA['Q'] = -1*(751.58*FRA['WL'] - 47.64) # regression from historical Q data

# Filter time
WLO = WLO[WLO.index >= UTC_starttime][['Q']]
ATC = ATC[ATC.index >= UTC_starttime][['Q']]
FRA = FRA[FRA.index >= UTC_starttime][['Q']]

# Save to disk
WLO.to_csv(name_WLO, index_label='datetime')
ATC.to_csv(name_ATC, index_label='datetime')
FRA.to_csv(name_FRA, index_label='datetime')

# TIDES-------------------------------------------------------------------------
NOAA_Tides = nc.Station(ID_NOAA_Tides) # Retrieve data

NOAA_Tides = NOAA_Tides.get_data(
    begin_date=UTC_starttime[0:16].replace('-',''),
    end_date=CDT_dates[1].replace('-',''),
    product="water_level",
    datum=datum,
    units="metric",
    time_zone="gmt")

NOAA_Tides = NOAA_Tides[['water_level']]
NOAA_Tides.to_csv(name_Tides, index_label='datetime') # Save

# WIND--------------------------------------------------------------------------
NOAA_Wind = nc.Station(ID_NOAA_Wind) # Retrieve data

NOAA_Wind = NOAA_Wind.get_data(
    begin_date=UTC_starttime[0:16].replace('-',''),
    end_date=CDT_dates[1].replace('-',''),
    product="wind",
    units="metric",
    time_zone="gmt")

NOAA_Wind = NOAA_Wind[['spd','dir']].rename(
    columns={'spd':'Speed','dir':'Direction'})
NOAA_Wind.to_csv(name_Wind, index_label='datetime') # Save

