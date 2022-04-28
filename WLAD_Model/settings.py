"""
Settings for running the SPRING DELTA-X CAMPAIGN in ANUGA
"""
# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
modelname = 'WLAD2_2103_7day'
Q_steady = False # Use constant or varying discharge
T_steady = False # Use constant or varying tidal water level
wind_on = False # Is wind on or off?
hot_start = False # Is this model run a hot or cold start?
save_new_IC = True # Save new IC from this run to use later?
timestep = 900.  # Time in seconds between timesteps
finaltime = 7*24*3600+timestep # End time in seconds

UTC_starttime = '2021-03-20 00:00:00' # Actual UTC start time for model
CDT_dates = ("2021-03-19", "2021-04-03") # Date range to grab BC data in CDT/CST

# BC GAUGE IDs------------------------------------------------------------------
ID_Calumet = '07381590' # USGS
ID_MorganCity = '07381600' # USGS
ID_Franklin = '07381670' # USGS
ID_LAWMA = '8764227' # NOAA
ID_Eugene = '8764314' # NOAA

# DISCHARGE---------------------------------------------------------------------
# Wax Lake Outlet (Calumet inlet)
Q_WLO_loc = [[653808.0,3298573.0],[654788.0,3298573.0]]
name_WLO = 'BoundaryConditions/Q_Calumet_%s-%s.csv' % (
           UTC_starttime[0:10].replace('-',''),CDT_dates[1].replace('-',''))
# Atchafalaya (Morgan City inlet)
Q_ATC_loc = [[674112.0,3290260.0],[672944.0,3290260.0]]
name_ATC = 'BoundaryConditions/Q_MorganCity_%s-%s.csv' % (
           UTC_starttime[0:10].replace('-',''),CDT_dates[1].replace('-',''))
# Franklin (GIWW west outlet)
Q_FRA_loc = [[646980.0,3285310.0],[646980.0,3284810.0]]
name_FRA = 'BoundaryConditions/Q_Franklin_%s-%s.csv' % (
           UTC_starttime[0:10].replace('-',''),CDT_dates[1].replace('-',''))
# Avoca (GIWW east outlet) (LACKS DISCHARGE GAUGE)
Q_AVO_loc = [[672735.0,3269176.0],[672235.0,3269176.0]]
# Avg of historical USGS ADCP measurements by month as follows:
# Q_AVO_data = [-29,-213,-459,-533,-655,-666,-499,-420,9,121,117,11]
Q_AVO = -511.84 # SPRING (avg of all Mar-Apr)

# TIDES-------------------------------------------------------------------------
steady_T_elev = 0.26 # WL elevation for steady tide runs, normal ~0.26
# Uncomment for LAWMA:
# tide_source = 'lawma'
# tide_timeoffset = "100m"
# tide_vertoffset = -0.05
# Uncomment for EUGENE:
tide_source = 'eugene'
tide_timeoffset = "45m"
tide_vertoffset = -0.3

ID_NOAA_Tides = int(ID_LAWMA if tide_source=='lawma' else ID_Eugene)
datum = 'NAVD' if tide_source=='lawma' else 'MLLW'
name_Tides = 'BoundaryConditions/Tides_%s_%s_%s-%s.csv' % (tide_source, datum,
             UTC_starttime[0:10].replace('-',''),CDT_dates[1].replace('-',''))

# WIND--------------------------------------------------------------------------
wind_source = 'eugene' # 'lawma' or 'eugene'
ID_NOAA_Wind = int(ID_LAWMA if wind_source=='lawma' else ID_Eugene)
name_Wind = 'BoundaryConditions/Wind_%s_%s-%s.csv' % (wind_source,
             UTC_starttime[0:10].replace('-',''),CDT_dates[1].replace('-',''))

# MESH SETTINGS-----------------------------------------------------------------
base_res = 625 # 25m res, background resolution outside interior polygons

# INITIAL CONDITIONS------------------------------------------------------------
# Input file for topographic data
topography_source = r'BathymetryPolygons/WLAD_topo.asc'
# Hot start files from a previous run
IC_depth_source = r'BoundaryConditions/IC_depth.csv'
IC_xmom_source = r'BoundaryConditions/IC_xmom.csv'
IC_ymom_source = r'BoundaryConditions/IC_ymom.csv'
IC_source_run = r'220209_1551_WLAD2_2103' # Indicate source of IC files for reproducibility
# Name of new IC files if saving for future restart
IC_depth_dest = r'BoundaryConditions/IC_depth_restart.csv'
IC_xmom_dest = r'BoundaryConditions/IC_xmom_restart.csv'
IC_ymom_dest = r'BoundaryConditions/IC_ymom_restart.csv'

# FRICTION----------------------------------------------------------------------
# Define friction parameters for classes 1-6, i.e.
# [ocean, channels, small channels, subtidal, intertidal, supratidal]
friction_loc = r'BoundaryConditions/FrictionMap.tif'
n_array = [0.015, 0.028, 0.005,     0,     0,     0]
m_array = [    0,     0,     0,   120,   200,   200]
h_array = [    0,     0,     0,   0.5,   1.0,    10]
D_array = [    0,     0,     0,  0.01,  0.01, 0.015]
