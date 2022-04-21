"""
Code to run the ANUGA simulation in parallel on TACC. Settings loaded
from settings.py. Assumes both prepare_*.py scripts have been run.
Submit to slurm using parallel.job
"""
# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------
from __future__ import division, print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cmocean
import os
import datetime
import anuga
from anuga.utilities import animate
from anuga import Inlet_operator, distribute, myid, numprocs, finalize, barrier
from anuga.operators.baptist_operator import Baptist_operator
from .settings import *
from .tools import *

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

# ------------------------------------------------------------------------------
# Setup boundary conditions
# ------------------------------------------------------------------------------
# TIDES-------------------------------------------------------------------------
if not T_steady:
    fBC_tides = GenerateTideGauge()
    # fBC_tides = GenerateTideCosine()

# WIND--------------------------------------------------------------------------
if wind_on:
    fBC_windspeed, fBC_winddir = GenerateWind()

# DISCHARGE---------------------------------------------------------------------
fQ_WLO = GenerateHydrograph(name_WLO) # Calumet inlet
fQ_ATC = GenerateHydrograph(name_ATC) # Morgan City inlet
fQ_FRA = GenerateHydrograph(name_FRA) # Franklin outlet
# Avoca outlet discharge constant and set by Q_AVO

# INITIAL CONDITIONS------------------------------------------------------------
if hot_start:
    init_depth = np.loadtxt(IC_depth_source, delimiter=',')
    init_xmom = np.loadtxt(IC_depth_source, delimiter=',')
    init_ymom = np.loadtxt(IC_depth_source, delimiter=',')

# ------------------------------------------------------------------------------
# Do the domain creation on processor 0
# ------------------------------------------------------------------------------
if myid == 0:
    # ------------------------------------------------------------------------------
    # Make working directory here so we don't repeat the process on each cpu
    # ------------------------------------------------------------------------------
    mydir = os.path.join('/Outputs/', 
                         datetime.datetime.now().strftime('%y%m%d_%H%M_' + modelname))
    os.makedirs(mydir)
    os.makedirs(mydir + '/figs')
    os.makedirs(mydir + '/data')

    # ------------------------------------------------------------------------------
    # Create domain and mesh
    # ------------------------------------------------------------------------------
    bounding_polygon, boundary_tags, inside_regions, geo_reference = GenerateDomainGeometry()

    domain = anuga.create_domain_from_regions(bounding_polygon, boundary_tags,
                                              maximum_triangle_area=base_res,
                                              interior_regions=inside_regions,
                                              poly_geo_reference=geo_reference,
                                              mesh_geo_reference=geo_reference,
                                              mesh_filename = 'WLAD2.msh')
    domain.set_name('WLAD2')
    domain.set_datadir(mydir)
    domain.set_flow_algorithm('DE1')
    domain.set_low_froude(1)  # Use low-froude DE1 to reduce flux damping
    domain.set_minimum_allowed_height(0.005)  # Only store heights > 0.5 cm

    # Plot mesh
    fig = plt.figure(figsize=(10, 12), dpi=200, facecolor='w', edgecolor='k')
    dplotter = animate.Domain_plotter(domain)
    plt.triplot(dplotter.triang, linewidth=0.1);
    plt.axis('scaled')
    plt.savefig(mydir + '/mesh.png')
    plt.close()

    print(domain.statistics())

    # ------------------------------------------------------------------------------
    # Set Inital Conditions
    # ------------------------------------------------------------------------------
    # ---------------Load pre-established elevation----------------------
    topo = np.loadtxt(topography_source.replace('.asc','.csv'), delimiter=",")

    stage = topo.copy()  # Initialize stage as = topography
    if hot_start:
        print('Source of IC files is %s' % IC_source_run)
        stage = topo + init_depth  # Add depth to topo
        # stage[topo <= fBC_tides(0)] = fBC_tides(0)
        # Initialize momentums if we have them:
        domain.set_quantity('xmomentum', init_xmom, location='centroids')
        domain.set_quantity('ymomentum', init_ymom, location='centroids')
    elif not T_steady:
        stage[topo <= fBC_tides(0)] = fBC_tides(0) # Used for cold start
    else:
        stage[topo <= steady_T_elev] = steady_T_elev
    domain.set_quantity('elevation', topo, location='centroids') 
    domain.set_quantity('stage', stage, location='centroids')  # Initialize depth

    # Plot the centroid elevation of the mesh cells
    msl = 0.3
    fig = plt.figure(figsize=(10, 10), dpi=400, facecolor='w', edgecolor='k')
    plt.tripcolor(dplotter.triang, facecolors=dplotter.elev,
                  vmax=msl+5, vmin=msl-5, cmap='cmo.topo')
    plt.colorbar();
    plt.title("Elevation");
    plt.axis('scaled')
    plt.savefig(mydir + '/elev.png')
    plt.close()

else:
    domain = None

# ------------------------------------------------------------------------------
# Produce parallel domain
# ------------------------------------------------------------------------------
domain = distribute(domain)
domain.set_store_vertices_uniquely(False)

#---------------------------------------------------------------------------
# Assign Friction
#---------------------------------------------------------------------------
x = domain.quantities['x'].centroid_values
y = domain.quantities['y'].centroid_values
FricVal = Raster2Mesh(x, y) # Interpolate from raster onto grid
n, m, hv, D = AssignFricValue(FricVal) # Assign to each cell

# Assign Mannings
domain.set_quantity('friction', n, location = 'centroids')
# Assign Baptist
bapt_op = Baptist_operator(domain,
                           veg_diameter=D,
                           veg_density=m,
                           veg_height=hv,
                           use_diffusivity=False)

# ------------------------------------------------------------------------------
# Apply forcings
# ------------------------------------------------------------------------------
if T_steady:
    Bout = anuga.Dirichlet_boundary([steady_T_elev, 0.0, 0.0]) # Constant (No Tides)
else:
    Bout = anuga.Time_boundary(domain, function = lambda t: [fBC_tides(t), 0.0, 0.0])
Br = anuga.Reflective_boundary(domain)  # Solid walls

domain.set_boundary({'bay': Bout, 'sides': Br})

if wind_on:
    W = anuga.Wind_stress(fBC_windspeed(0), fBC_winddir(0))
    domain.forcing_terms.append(W) # Apply wind

# Setup inlets
inlet_WLO = Inlet_operator(domain, Q_WLO_loc, Q = fQ_WLO(0))
inlet_ATC = Inlet_operator(domain, Q_ATC_loc, Q = fQ_ATC(0))
outlet_FRA = Inlet_operator(domain, Q_FRA_loc, Q = fQ_FRA(0))
outlet_AVO = Inlet_operator(domain, Q_AVO_loc, Q = Q_AVO)

# ------------------------------------------------------------------------------
# Evolve system through time
# ------------------------------------------------------------------------------
barrier()
# Evolve the domain
for n, t in enumerate(domain.evolve(yieldstep=timestep, finaltime=finaltime)):
    if myid == 0:
        domain.print_timestepping_statistics()
    
    if not Q_steady:
        # Update discharge
        try:
            inlet_WLO.Q = fQ_WLO(t)
            print('Updated WLO Discharge on Node %s' % myid)
        except:
            pass
        try:
            inlet_ATC.Q = fQ_ATC(t)
            print('Updated ATC Discharge on Node %s' % myid)
        except:
            pass
        try:
            inlet_FRA.Q = fQ_FRA(t)
            print('Updated FRA Discharge on Node %s' % myid)
        except:
            pass

    if wind_on:
        # Update wind field
        W.speed = fBC_windspeed(t+timestep/2) # Use value centered in each yieldstep
        W.phi = fBC_winddir(t+timestep/2)
barrier()

# ---------------------------------------------------------------
# Finalize
# ---------------------------------------------------------------
domain.sww_merge(delete_old=True)
finalize()

# ------------------------------------------------------------------------------
# Save a couple extra outputs on processor 0
# ------------------------------------------------------------------------------
if myid == 0:
    swwvals = anuga.utilities.plot_utils.get_centroids(os.path.join(mydir,'WLAD2.sww'),
                                                       timeSlices = 'all')
    # Query values: time, x, y, stage, elev, height, xmom, ymom, xvel, yvel, friction, vel, etc
    SaveOutputs(swwvals, mydir) # Save individual output files
    if save_new_IC:
        SaveInitialConditions(swwvals, mydir, -1,
                              IC_depth_dest,IC_xmom_dest,IC_ymom_dest) # Save new IC
