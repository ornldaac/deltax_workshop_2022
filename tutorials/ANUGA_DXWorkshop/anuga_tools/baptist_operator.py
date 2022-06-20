"""
Vegetation operators
"""
from __future__ import division, absolute_import, print_function
import numpy as num
from numpy import sqrt, minimum, maximum, log

# from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator

# from math import sqrt, log
# from anuga.config import epsilon, g

# import anuga.utilities.log as log
# from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            # netcdf_float

# import os
# from scipy.interpolate import NearestNDInterpolator

#===============================================================================
# Vegetation operator applying drag to the flow
#===============================================================================
class Baptist_operator(Operator):
    """
    Baptist operator that applies a drag on the flow due to the presence of veg.
    Methodology is based on Baptist et al, 2007 for emergent and submerged veg.
    Modified Chezy coefficient is computed as:
    Cv = (Cb^-2 + (Cd*m*D/(2g))*min(h,hv))^-0.5 + (g^0.5/k)*ln(max(h,hv)/hv)
    
    (Cv is vegetated Chezy, Cb is bed Chezy, Cd is drag coefficient, m is stem
    density [m^-2], D is stem diameter [m], g is gravity, h is flow depth [m], 
    hv is stem height [m], k is von karman constant)
    
    However, we precompute some of these terms together in the init for speed.
    The result looks like:
    Cv = (a1 + a2*CD*min(h,hv))^-0.5 + a3*ln(max(h,hv)/hv)
    
    Operator uses explicit form to update velocities, according to:
    d/dt(uh) = -g*uh*sqrt(uh^2 + vh^2)/(Cv^2*h^2)
    """

    def __init__(self, 
                 domain, 
                 veg_diameter=None,
                 veg_density=None,
                 veg_height=None,
                 bed_friction=65.0,
                 use_diffusivity=False,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
        """

        Initialize vegetation characteristics

        **Inputs** :

            domain : `object`
                ANUGA domain instance

            veg_diameter : `float or np.ndarray`
                Vegetation stem diameters, given in meters. Input can be either
                one constant value to be applied everywhere, or an array giving
                one value per cell centroid. If None is given, we check to see
                if values were previously defined in domain.quantities()

            veg_density : `float or np.ndarray`
                Vegetation stem density, given in number per unit area [#/m^2].
                Input can be either one constant value to be applied everywhere,
                or an array giving one value per cell centroid. If None is
                given, we check to see if values were previously defined in
                domain.quantities()

            veg_height : `float or np.ndarray`
                Vegetation stem height, given in meters. Input can be either
                one constant value to be applied everywhere, or an array giving
                one value per cell centroid. If None is given, we check to see
                if values were previously defined in domain.quantities()

            bed_friction : `float or np.ndarray`, optional
                Bed friction Chezy coefficient. Default is 65. Input can be either
                one constant value to be applied everywhere, or an array giving
                one value per cell centroid.

            use_diffusivity : `bool`
                Indicates whether to turn on vegetated kinematic viscosity, or
                to just use friction. Default is true.

        **Outputs** :

            all_walk_data : `dict`
                Dictionary of all x and y locations and travel times, with
                details same as input previous_walk_data

        """
        Operator.__init__(self, domain, description, label, logging, verbose)

        #-----------------------------------------------------
        # Pull domain information
        self.depth = self.stage_c - self.elev_c
        self.g = self.domain.g

        #-----------------------------------------------------
        # Populate some variables
        if veg_diameter is None:
            try:
                self.veg_diameter = self.domain.get_quantity('veg_diameter')
            except:
                raise ValueError('Vegetation diameter not defined')
        else:
            if isinstance(veg_diameter, float):
                veg_diameter = num.ones_like(self.depth, dtype=float)*veg_diameter
            self.veg_diameter = veg_diameter
            Quantity(domain, name='veg_diameter', register=True)
            self.domain.quantities['veg_diameter'].\
                set_values(self.veg_diameter, location = 'centroids')

        if veg_density is None:
            try:
                self.veg_density = self.domain.get_quantity('veg_density')
            except:
                raise ValueError('Vegetation spacing not defined')
        else:
            if isinstance(veg_density, float):
                veg_density = num.ones_like(self.depth, dtype=float)*veg_density
            self.veg_density = veg_density
            Quantity(domain, name='veg_density', register=True)
            self.domain.quantities['veg_density'].\
                set_values(self.veg_density, location = 'centroids')

        if veg_height is None:
            try:
                self.veg_height = self.domain.get_quantity('veg_height')
            except:
                raise ValueError('Vegetation height not defined')
        else:
            if isinstance(veg_height, float):
                veg_height = num.ones_like(self.depth, dtype=float)*veg_height
            self.veg_height = veg_height
            Quantity(domain, name='veg_height', register=True)
            self.domain.quantities['veg_height'].\
                set_values(self.veg_height, location = 'centroids')

        self.use_diffusivity = use_diffusivity
        if self.use_diffusivity:
            self.diffusivity = num.zeros((self.num_cells,))
            self.mix_length = num.zeros((self.num_cells,))

            try:
                diff = self.domain.get_quantity('diffusivity')
            except:
                Quantity(domain, name='diffusivity', register=True)

        self.domain.set_use_kinematic_viscosity(self.use_diffusivity)

        #----------------------------------------------------
        # Start populating/computing coefficients
        self.ind = (self.veg_density > 0) # Vegetated indices

        self.Cd = 1.68 # Assumed for now, replace later with equation
        if isinstance(bed_friction, float): # Check if Cb is const or array
            self.bed_friction = num.ones_like(self.depth, dtype=float)*bed_friction
        self.a1 = self.bed_friction**-2 # First lumped coefficient

        self.a2 = (self.veg_diameter * self.veg_density / (2*self.g)) # Second lumped coefficient

        self.a3 = 7.8289 # Third lumped coefficient



    def __call__(self):
        """
        Apply vegetation drag according to veg_diameter and veg_density quantities
        """
        # Get the timestep for explicit update
        self.dt = self.get_timestep()
        
        self.depth = self.stage_c - self.elev_c
        self.ind = (self.veg_density > 0) & (self.depth > 0.01) # Update indices

        self.update_quantities()
        
        
        
    def update_quantities(self):
        """
        Calculate the drag that vegetation imparts on the flow
        and update momentum quantities
        """
    
        if sum(self.ind)>0:
            # Cut down some variables to just vegetated areas
            self.depth_w = self.depth[self.ind]
            hv_w = self.veg_height[self.ind]
            Cv_w = num.zeros_like(self.depth_w) # Vegetated Chezy
            a1_w = self.a1[self.ind]
            a2_w = self.a2[self.ind]
            xmom_w = self.xmom_c[self.ind]
            ymom_w = self.ymom_c[self.ind]
                    
            if self.use_diffusivity:
                self.mix_length[:] = 0
                self.diffusivity[:] = 0
                self.calculate_diffusivity()

            # calculate discharge in the cell
            qcell_w = sqrt(xmom_w*xmom_w + ymom_w*ymom_w)

            # Calculate Chezy
            Cv_w = (a1_w + a2_w * self.Cd * minimum(self.depth_w, hv_w))**-0.5 + self.a3 * log(maximum(self.depth_w, hv_w) / hv_w)

            # Compute friction slope
            Sf_x = self.g * xmom_w * qcell_w / (Cv_w**2 * self.depth_w**2 + 1e-6)
            Sf_y = self.g * ymom_w * qcell_w / (Cv_w**2 * self.depth_w**2 + 1e-6)

            self.xmom_c[self.ind] = xmom_w - Sf_x * self.dt
            self.ymom_c[self.ind] = ymom_w - Sf_y * self.dt



    def calculate_drag_coefficient(self):
        """
        Calculate the drag coefficient Cd as a function of ad using
        the curve fitted to Figure 6 in Nepf (1999)
        """
        
        self.Cd = (56.11 * self.ad**2
                   - 15.28 * self.ad
                   + 1.3
                   - 0.0005465 / self.ad)
        self.Cd[self.ad < 0.006] = 1.2
        
        self.Cd_veg = 0.5 * self.Cd * self.alpha
        

        
    def calculate_diffusivity(self):
        """
        Calculates a value for the quantity diffusivity for use with
        kinematic viscosity operator.
        
        Based on the formulation of Nepf (1999)
        
        For all cells, the transition from mix_length = depth to a linear
        approximation of the mixing length happens at veg_density = depth,
        which is ad = diameter**2 / depth**2
        """
        
        ad_deltaS = self.veg_d_w**2 / self.depth_w**2
        
        mix_length_slope = (self.veg_d_w - self.depth_w) / (0.01 - ad_deltaS)
        
        self.mix_length = (mix_length_slope * self.ad_w +
                          (self.depth_w - mix_length_slope * ad_deltaS))
        
        ind = self.veg_s_w > self.depth_w
        self.mix_length[ind] = self.depth_w[ind]
        
        ind = self.ad_w > 0.01        
        self.mix_length[ind] = self.veg_d_w[ind]

        # turbulent kinetic energy
        Cb = 0.001
        k = ((1 - self.ad_w) * Cb + (self.Cd[self.ind] * self.ad_w)**0.66) * self.velocity**2
        
        
        # total diffusivity
        self.diffusivity[self.ind] = (num.sqrt(k)**0.5 * self.mix_length +
                                      self.ad_w * self.velocity * self.veg_d_w)
                    
                    
        self.domain.quantities['diffusivity'].\
                set_values(self.diffusivity, location = 'centroids')
            


    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Veg_operator, time '
        message += str(self.get_time())
        return message