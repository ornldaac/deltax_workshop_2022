"""
Background functions used to setup ANUGA simulation
"""
# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------
from __future__ import division
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from osgeo import gdal
import os
import glob
import anuga
from settings import *

# ------------------------------------------------------------------------------
# Boundary condition functions
# ------------------------------------------------------------------------------
def GenerateTideGauge(filename=name_Tides,
                      time_offset=tide_timeoffset,
                      vert_offset=tide_vertoffset):
    """
    Function to generate a tidal BC from gauge data.
    Specify inputs in settings.py. Returns a function of time.
    """
    # Measured water levels
    tides = pd.read_csv(filename, header=0, names=['datetime', 'WL'])
    # Some correction factors to apply due to using inland gauge as BC out in bay
    time_offset = pd.Timedelta(time_offset) # Correct tide arrival time
    # Convert date string to datetime, create new column epoch for int numeric time
    tides['datetime'] = pd.to_datetime(tides['datetime'])
    tides['epoch'] = (tides['datetime'] - tides['datetime'][0] - time_offset) // pd.Timedelta("1s")
    fBC_tides = interp1d(tides['epoch'], tides['WL'] + vert_offset, kind='linear')
    return fBC_tides

def GenerateTideCosine(amplitude=0.25, period=7.2722e-5, phase=0, offset=0.26):
    """
    Function to generate an artificial tidal BC from a sinusoid. Default is
    cosine with same tidal range as LAWMA, freq of N2.
    Returns a function of time.
    """
    # Default is pseudo tide with same tidal range as LAWMA, freq of N2. Time in sec
    def fBC_tides(t):
        return amplitude*np.cos(t*period + phase) + offset
    return fBC_tides

def GenerateWind(filename=name_Wind):
    """
    Function to generate a wind forcing from gauge data.
    Specify inputs in settings.py. Returns two functions of time.
    In output, 0 degrees is from the West, 90 degrees is from the North.
    """
    # NOAA Amerada Pass, [m/s], direction FROM (0 is from N, 90 is from E, etc)
    met = pd.read_csv(filename, header = 0,
                      names=['datetime', 'speed', 'dir'])
    # Convert to anuga's weird wind direction convention
    wind_dir = -1*(np.array(met['dir']) + 90)
    met['datetime'] = pd.to_datetime(met['datetime'])
    met['epoch'] = (met['datetime'] - met['datetime'][0]) // pd.Timedelta("1s")
    # Linearly interpolate wind speeds at each time:
    fBC_windspeed = interp1d(met['epoch'], met['speed'], kind='linear')
    # Use previous 6-min direction value over window:
    fBC_winddir = interp1d(met['epoch'], wind_dir, kind='previous')
    return fBC_windspeed, fBC_winddir

def GenerateHydrograph(filename, steady=Q_steady):
    """
    Function to generate a hydrograph from USGS gauge data.
    Specify inputs in settings.py. Returns a function of time.
    When calling, filename either name_WLO, name_ATC, name_FRA
    """
    Q = pd.read_csv(filename, header=0, 
                    names = ['datetime','Q'])
    if steady:
        def fQ(t):
            return np.mean(Q['Q'])
    else:
        Q['datetime'] = pd.to_datetime(Q['datetime'])
        Q['epoch'] = (Q['datetime'] - Q['datetime'][0]) // pd.Timedelta("1s")
        fQ = interp1d(Q['epoch'], Q['Q'], kind='linear')
    return fQ

# ------------------------------------------------------------------------------
# Domain settings functions
# ------------------------------------------------------------------------------
def Raster2Mesh(meshX, meshY, raster_filename=friction_loc):
    """
    Function to grab raster value at each mesh cell centroid.
    Used for friction assignment. Specify raster in settings.py.
    Returns a numpy.ndarray
    """
    src = gdal.Open(raster_filename)
    data_array = np.array(src.GetRasterBand(1).ReadAsArray())
    # Get geographic info
    transform = src.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]

    # Loop and grab all values at mesh coords
    meshVal = np.zeros(len(meshX), dtype=int)
    for ii in list(range(0, len(meshX))):
        col = int((meshX[ii] - xOrigin) / pixelWidth)
        row = int((yOrigin - meshY[ii] ) / pixelHeight)
        try:
            meshVal[ii] = data_array[row][col]
        except IndexError:
            meshVal[ii] = 0
    # Return values
    return meshVal

def AssignFricValue(FricVal,
                    n_array=n_array,
                    m_array=m_array,
                    h_array=h_array,
                    D_array=D_array):
    """
    Using Friction ID (output of Raster2Mesh), assign friction parameters.
    All Mannings (n) and Baptist (m, hv, D) values returned.
    Specify coefficients for each class in settings.py
    """
    FricVal = FricVal.astype(int)-1 # Subtract 1 assuming map ID's start at 1
    n = np.zeros_like(FricVal, dtype=float)
    m = np.zeros_like(FricVal, dtype=float)
    hv = np.zeros_like(FricVal, dtype=float)
    D = np.zeros_like(FricVal, dtype=float)
    
    for ii in list(range(len(FricVal))):
        n[ii] = n_array[FricVal[ii]]
        m[ii] = m_array[FricVal[ii]]
        hv[ii] = h_array[FricVal[ii]]
        D[ii] = D_array[FricVal[ii]]
    return n, m, hv, D

def GenerateDomainGeometry():
    """
    Function to generate inputs for the domain instance.
    No inputs required, all these details are hard-coded.
    """
    # ---------------Create boundaries-------------------
    bounding_polygon = anuga.read_polygon('BathymetryPolygons/WLAD_Boundary.csv')
    boundary_orig = [min([item[0] for item in bounding_polygon]), 
                     min([item[1] for item in bounding_polygon])]

    boundary_tags={'bay': [9], 
                   'sides': [0,1,2,3,4,5,6,7,8,10,11,12,13,14,15,
                             16,17,18,19,20,21,22,23,24,25,26,27]}

    # ---------------Load in polygons-------------------
    polygon_files = glob.glob(r'BathymetryPolygons/*Reg*.csv')
    inside_regions = []
    for poly in polygon_files:
        polyres = int(poly.split('_Res')[-1].replace('.csv',''))
        inside_regions.append([anuga.read_polygon(poly), polyres])

    # ---------------Define Geo Reference-------------------
    geo_reference = anuga.Geo_reference(zone=15,
                                        datum='wgs84',
                                        projection='UTM',
                                        false_easting=500000,
                                        false_northing=0)

    return bounding_polygon, boundary_tags, inside_regions, geo_reference

# ------------------------------------------------------------------------------
# Post-processing functions
# ------------------------------------------------------------------------------
def SaveOutputs(swwvals, mydir):
    """
    Function to extract and save individual output files in binary format.
    Inputs are the 'mydir' defined in the run script and reloaded .sww
    values, i.e. the output of anuga.utilities.plot_utils.get_centroids()
    """
    # Query values: time, x, y, stage, elev, height, xmom, ymom, xvel, yvel, friction, vel, etc
    np.save(mydir+'/data/time.npy', swwvals.time.data, allow_pickle=False)
    np.save(mydir+'/data/topo.npy', swwvals.elev.data, allow_pickle=False)
    np.save(mydir+'/data/x.npy', swwvals.x.data, allow_pickle=False)
    np.save(mydir+'/data/y.npy', swwvals.y.data, allow_pickle=False)
    np.save(mydir+'/data/xmom.npy', swwvals.xmom, allow_pickle=False)
    np.save(mydir+'/data/ymom.npy', swwvals.ymom, allow_pickle=False)
    np.save(mydir+'/data/xvel.npy', swwvals.xvel, allow_pickle=False)
    np.save(mydir+'/data/yvel.npy', swwvals.yvel, allow_pickle=False)
    np.save(mydir+'/data/depth.npy', swwvals.height, allow_pickle=False)
    np.save(mydir+'/data/stage.npy', swwvals.stage, allow_pickle=False)
    # Extract friction values
    FricVal = Raster2Mesh(swwvals.x, swwvals.y, friction_loc)
    n, m, hv, D = AssignFricValue(FricVal)
    a2 = (D*m/(2*9.81)) # Second lumped coefficient
    chezy = np.zeros_like(swwvals.height)
    for ii in list(range(len(swwvals.time))): 
        idn = n>0
        idm = m>0
        chezy[ii,idn] = (swwvals.height[ii,idn]**(1/6))/n[idn]
        chezy[ii,idm] = (65.**-2 + a2[idm]*1.68*np.minimum(swwvals.height[ii,idm], hv[idm]))**-0.5 + \
                        7.8289*np.log(np.maximum(swwvals.height[ii,idm], hv[idm]) / hv[idm])
    chezy[chezy==0] = 65.
    np.save(mydir+'/data/chezy.npy', chezy, allow_pickle=False)
    return

def SaveInitialConditions(swwvals, mydir, time_id=-1,
                          name_depth=IC_depth_source,
                          name_xmom=IC_xmom_source,
                          name_ymom=IC_ymom_source):
    """
    Function to extract and save a specific timestep to use for future IC
    Inputs are reloaded .sww file and IC filenames, defaults set
    in settings.py. Default time_id is last timestep.
    """
    import shutil
    # Replace old IC files with new ones
    np.savetxt(name_depth, swwvals.height[time_id,:], delimiter=',')
    np.savetxt(name_xmom, swwvals.xmom[time_id,:], delimiter=',')
    np.savetxt(name_ymom, swwvals.ymom[time_id,:], delimiter=',')
    # Then copy them into folder to keep a record of IC files
    newdir = os.path.join(mydir,'New_IC_Files')
    os.makedirs(newdir)
    shutil.copy2(name_depth,os.path.join(newdir,name_depth.split(os.sep)[-1]))
    shutil.copy2(name_xmom,os.path.join(newdir,name_xmom.split(os.sep)[-1]))
    shutil.copy2(name_ymom,os.path.join(newdir,name_ymom.split(os.sep)[-1]))
    return
