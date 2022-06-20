"""
Prepare topography files for ANUGA simulation
Interpolates .ASCII raster onto unstructured mesh 
"""
# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------
import numpy as np
import anuga
from .settings import topography_source, base_res
from .tools import GenerateDomainGeometry

# FILE CONVERSIONS--------------------------------------------------------------
# Note: If starting from GeoTIFF, need to convert to ASCII. Optional command below:
# gdal_translate -of AAIGrid  input_bathy.tif input_bathy.asc
# ASCII file needs attributes ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value
# Convert from ASCII to DEM
anuga.asc2dem(topography_source, use_cache=False, verbose=True)
# Convert from DEM to PTS
anuga.dem2pts(topography_source.replace('.asc','.dem'),
              use_cache=False, verbose=True)

# CREATE MESH-------------------------------------------------------------------
bounding_polygon, boundary_tags, inside_regions, geo_reference = GenerateDomainGeometry()

domain = anuga.create_domain_from_regions(bounding_polygon, boundary_tags,
                                          maximum_triangle_area=base_res,
                                          interior_regions=inside_regions,
                                          poly_geo_reference=geo_reference,
                                          mesh_geo_reference=geo_reference,
                                          mesh_filename = 'WLAD2.msh')

print domain.statistics()
# Save list of cell areas
np.savetxt('areas.csv', domain.areas, delimiter=",")

# INTERPOLATE-------------------------------------------------------------------
domain.set_quantity('elevation', 
                    filename=topography_source.replace('.asc','.pts'),
                    use_cache=False,
                    verbose=True,
                    alpha=0.1)

topo = domain.quantities['elevation'].centroid_values
# Save to disk
np.savetxt(topography_source.replace('.asc','.csv'), topo, delimiter=",")
