# 2024-02-16, Caroline Haslebacher
# open each polygon that we analyse,
# calculate number of sub-polygons (not trivial with rings and holes!)
# add a unique fixed ID
# save again

#%% import

import os
from pathlib import Path
import matplotlib.pyplot as plt
import json
from geojson import dump
# from osgeo import gdal
import numpy as np
from datetime import datetime
import argparse

from MAIN_loop1_to_loop4 import get_polygeopaths
from science_with_manual_segs_algo import load_tiff_arr

#%% define path to titania

# current = os.getcwd()
# titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

#%% get polygon paths

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Get basepath.')

    parser.add_argument('--basepath', required=True,
                    metavar="where to find /geotiffs/for_analysis and /polygons/for_analysis",
                    help='as posix')
    parser.add_argument('--savefoldername', required=True,
                    metavar="where to save changed geojson from /polygons/for_analysis",
                    help='string, e.g. Equirectangular')

    args = parser.parse_args()
    ''' for debugging:
    class args_init:
        basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/azimuth_analysis/Viviane_fullRegmaps/data'
        savefoldername = 'Equirectangular'
    args = args_init()
    current = os.getcwd()
    titaniach = Path(current.split('Caroline')[0]) / 'Caroline'
    '''

    geopoly_basepath = Path(args.basepath) # titaniach / 'lineament_detection/galileo_manual_segmentation/'
    # equirectangular
    geo_orig = geopoly_basepath / 'geotiffs/for_analysis'
    poly_orig = geopoly_basepath / 'polygons/for_analysis'
    geopaths, polypaths = get_polygeopaths(poly_orig, geo_orig)

    savegeojson = geopoly_basepath / 'polygons/for_analysis' / args.savefoldername
    os.makedirs(savegeojson, exist_ok=True)

    ###### assign fixed id (date and miliseconds?)

    for polypath in polypaths:
        print(polypath)
        # open geojson
        with open(polypath) as jsf:
            geoms = json.load(jsf)
        # loop through features
        for i, feature in enumerate(geoms['features']):
            # assign a unique, fixed ID (UFID)
            UFID = '{}_{}'.format(polypath.stem, i)# NOTE: old format took only first entry of the polypath (format(polypath.stem.split('_')[0], i)). Now we take all to allow for several polygon files for one geotiff
            feature['properties']['UFID'] = UFID

            holes_count = 0 # count number of holes!
            # print(feature['geometry'].keys())
            try:
                for fts_idx in range(len(feature['geometry']['coordinates'])):
                    # print(fts_idx)
                    holes = False # set holes index to False for each part
                    
                    p = feature['geometry']['coordinates'][fts_idx]
                    # print(p)
                    # print(len(p))
                    # print(UFID)
                    # print('\n')
                    # for numpair in p:
                    #     print(numpair)

                    if len(p) == 1:
                        p = p[0] # then, we need to reduce the dimension
                        # now, len(p) are the number of points
                        # print('length after reduction: {}'.format(len(p)))
                    # Handling holes (from conversion of shapefiles to geojson)
                    # if len(p) still 2 or 3, then, this means it has 1 or 2 holes.
                    # we can better check the shape of parr

                    # convert to array
                    parr = np.asarray(p, dtype='object') # against deprectation Warning
                    # HOLES: in the case of holes, parr.shape is for example (3,) instead of e.g. (5,2).
                    # also, numpy throws a warning: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray.
                    if len(parr.shape) == 1:
                        # change index
                        holes = True
                        # consider the holes count
                        holes_count += parr.shape[0] - 1 # e.g. 3 parts, means 2 holes (1 is the outer part)
            except TypeError:
                print('error. I try to skip this')
                # for example, for empty features
                continue

            # number of subpolygons, minus number of holes
            num_subpolygons = len(feature['geometry']['coordinates']) - 1 + holes_count
            # new entry in the geojson data:
            feature['properties']['num_subpolygons'] = num_subpolygons

        # save and close polygon again
        with open(savegeojson.joinpath(polypath.name), 'w') as f:
            dump(geoms, f)



# %%
