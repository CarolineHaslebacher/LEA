# very quick and dirty script for tiling up an input geotiff into smaller tifs based on a cut_size.
# adapted from LineaMapper_to_img_cutsize.py
# Additionally to the version in Z:\Groups\PIG\Caroline\lineament_detection\pytorch_maskrcnn\projects\code, 
# I added the ability to tile up feature layers
# I needed this script to tile up 17ESREGMAP02, the huge map, because it did not fit into the memory for an analysis

#%% import
from osgeo import gdal, gdalnumeric, ogr, osr # needed to be put first, otherwise it led to DLL import errors
import os
from pathlib import Path

# import albumentations as A
from datetime import datetime

#%matplotlib inline
import matplotlib.pyplot as plt
from PIL import ImageColor, ImageDraw, ImageFont
import torchvision.transforms.functional as F

import numpy as np
import time

import json
import numpy as np
import skimage.io
import math
import scipy

from PIL import Image
import cv2

import pickle
import argparse
import warnings


from geojson import MultiPolygon, Feature, FeatureCollection, dump, Point

from geojson.geometry import Point

from MAIN_loop1_to_loop4 import get_polygeopaths
from science_with_manual_segs_algo import load_tiff_arr, coordm_to_lon

#%%

# current = os.getcwd()
# titaniach = Path(current.split('Caroline')[0]) / 'Caroline'



#%% functions

def get_geotransform(dataset):
    # from here: https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
    # ulx, uly is the upper left corner, lrx, lry is the lower right corner
    ulx, xres, xskew, uly, yskew, yres  = dataset.GetGeoTransform()
    lrx = ulx + (dataset.RasterXSize * xres)
    lry = uly + (dataset.RasterYSize * yres)

    return ulx, xres, xskew, uly, yskew, yres, lrx, lry

def geotiff_to_arr(geopath):
    '''
        Reads a full path to a TIFF file.
        Opens it with GDAL as a dataset.
        Normalizes the array
        Sets fill values to -0.01
    '''
    geotiff_path = geopath
    print(geotiff_path)
    # open geotiff
    dataset = gdal.Open(geotiff_path.as_posix(), gdal.GA_ReadOnly)  #, gdal.GA_ReadOnly not to change the file!

    # transform to numpy array
    arr = dataset.ReadAsArray()
    arr = np.array(arr, dtype='float32')
    # there are fill values (=-3.4028227e+38) we want to get rid of
    # now it can happen that 0 is 'the wrong minimum', and the contrast gets small
    # therefore, we can't simply use arr[arr<0] = 0
    # but we want to figure out the minimum and the maximum and rescale that to 0-1, in the best case.
    arrmin = np.min(arr[arr>0])
    # first, get rid of fill values
    arr[arr<0] = arrmin - 0.01 # the minus 0.01 makes sure that we are smaller than the minimum
    # lets' rescale now
    arr = (arr - arr.min())/(arr-arr.min()).max()

    if len(arr.shape) > 2:
        print('WARNING: I only take the first band. Fingers crossed it is the right one. Otherwise please implement functionality of RGB channels.')
        arr = arr[0]
        
    return arr, dataset


def moving_window_tiling_ref(arr, row_step_n, col_step_n, tile_size, sort_zero=False, pad=False, xoffset=0, yoffset=0):
    '''
        pad: If True, the output tiles in the last row and column, where usually not a full tile shift fits in, is padded with zeros (in order to stack them).
             If False, the output tile size is conserved and more overlap is generated for the last row and column.
    example for tests:
    row_step_n = 100
    col_step_n = 100
    tile_size = 110
    sort_zero=False
    pad=True
    '''
    img_length = arr.shape[1]
    img_height = arr.shape[0]
    # calculate row and col partitioning
    # img_height - 1 not to run into 'out of bounds' index error --> but I think this leads later to other problems!
    rowsp = (img_height   - tile_size)/row_step_n # e.g. = 2.97, then we shift two full times plus one time only 0.97 of row_step
    colsp = (img_length  - tile_size)/col_step_n

    # initialize list of arrays
    tiles = []
    positions = [] # to store the x- and y-position of the moved file
    # end-of-row indicator
    end_of_row = False
    end_of_col = False

    for row_idx_n in range(math.ceil(rowsp)+1):
        # print(row_idx_n)
        # round up and test in every loop if row_idx is more than 1 away of 'rowsp' value
        if row_idx_n < rowsp:
            # print('row index: {}'.format(row_idx))
            row_idx = row_idx_n
        else:
            # print('end of row')
            end_of_row = True # set indicator to true
            if pad:
                row_idx = row_idx_n
            else:
                row_idx = rowsp # int((rowsp - row_idx_n + 1)*100)
                # print('row_idx for last row: {}'.format(row_idx))

        for col_idx_n in range(math.ceil(colsp)+1): # goes from 1,2,3 if colsp is 2.79 for example
            # print(col_idx_n)
            if col_idx_n < colsp:
                # print(col_idx_n)
                col_idx = col_idx_n
            else:
                # print('end of column')
                end_of_col = True # set indicator to true
                if pad:
                    col_idx = col_idx_n
                else:
                    # print('col_index: {}'.format(col_idx))
                    col_idx = colsp # int((colsp - col_idx + 1)*100)
                    # print('col_idx for last column: {}'.format(col_idx))

            positions.append((round(row_idx*row_step_n) + xoffset, round(col_idx*col_step_n) + yoffset))
            # print(positions)

            if (end_of_col or end_of_row) and pad:
                # for every end of row and if we would like to pad, treat specially with np.pad
                # print('padding')
                # I am sure there is a slightly better way of doing this, for example by using a variable for the 'tile_size'  
                if end_of_row and not end_of_col:
                    temp_tile = arr[round(row_idx*row_step_n):, round(col_idx*col_step_n):round(col_idx*col_step_n + tile_size)]
                    pad_rows = tile_size - temp_tile.shape[0] # the number of rows
                    pad_cols = 0
                elif end_of_col and not end_of_row:
                    temp_tile = arr[round(row_idx*row_step_n):round(row_idx*row_step_n + tile_size), round(col_idx*col_step_n):]
                    pad_rows = 0
                    pad_cols = tile_size - temp_tile.shape[1] # the number of cols
                elif end_of_row and end_of_col:
                    temp_tile = arr[round(row_idx*row_step_n):, round(col_idx*col_step_n):]
                    pad_rows = tile_size - temp_tile.shape[0] # the number of rows
                    pad_cols = tile_size - temp_tile.shape[1] # the number of cols
                # pad now
                # print(pad_rows, pad_cols)
                tile = np.pad(temp_tile, ((0, pad_rows), (0, pad_cols)), 'constant', constant_values=(0, 0) ) # Number of values padded to the edges of each axis. ((before_1, after_1), ... (before_N, after_N))
            else: # normal indexing
                tile = arr[round(row_idx*row_step_n):round(row_idx*row_step_n + tile_size),
                        round(col_idx*col_step_n):round(col_idx*col_step_n + tile_size) ]
            

            if sort_zero == True:
                # check if tile is totally empty. We do not want zero arrays! just do not append them and give a warning
                if np.equal(tile, np.zeros((tile_size, tile_size))).all() == False: # then, array is not equal to np.zeros() and we do append
                    tiles.append(tile)
            else:
                # just append
                # print(tile.shape)
                tiles.append(tile)

            # set col indicator to False again
            end_of_col = False
        # set row indicator to False again
        end_of_row = False

    # return an array with the tiles in the first axis
    if len(tiles) == 0:
        raise Warning
    else:
        return np.stack(tiles), positions
    
#%%


# # simply provide the full path 
# # to the GEOTIFF
# basepath_geotiff = titaniach / 'lineament_detection/RegionalMaps_CH_NT_EJL_LP/mapping/data/geotiff'
# filename = '17ESREGMAP02_GaliloSSI_Equi-cog_training_region'
# input_geotiff = basepath_geotiff.joinpath(filename + '.tif')
# # FEATURE layer
# basepath_featurelayer= titaniach / 'lineament_detection/RegionalMaps_CH_NT_EJL_LP/mapping/data/geojson'
# # same filename!
# input_featurelayer = basepath_featurelayer.joinpath(filename + '.geojson')

# # provide the path to the folder where to save (as a string)
# savepath_geotiff = basepath_geotiff / 'training_region_tiled_1000px'
# savepath_featurelayer_geojson = basepath_featurelayer / 'training_region_tiled_1000px/geojson'
# savepath_featurelayer_shp = basepath_featurelayer / 'training_region_tiled_1000px/shape_files'
# os.makedirs(savepath_geotiff, exist_ok=True)
# os.makedirs(savepath_featurelayer_geojson, exist_ok=True)
# os.makedirs(savepath_featurelayer_shp, exist_ok=True)
# # provide the cut_size
# cut_size = 500

#%%


if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Get basepath.')

    parser.add_argument('--basepath', required=True,
                    metavar="where to find /geotiffs/for_analysis and /polygons/for_analysis",
                    help='as posix')
    # parser.add_argument('--jsonpath', required=True,
    #                 metavar="where to find the crs_file to read in",
    #                 help='as posix')
    parser.add_argument('--crs_add', required=False,
                        default=False,
                        type=bool,
                    metavar="State if you want a coordinate system to be added to a feature layer. NOTE: only if this argument is omitted will it be false!",
                    help='True or False (0)')   
    parser.add_argument('--crs_file', required=False,
                        default='Equirectangular_EUROPA.prj',
                    metavar="the name of the crs file. Folder is jsonpath. No reprojection is happening.",
                    help='jsonpath.joinpath(equifile) has to be valid. crs_add has to be True.')   
    parser.add_argument('--savefoldername', required=False,
                        default='Equirectangular',
                    metavar="the name of sub-folder to save",
                    help='e.g. Equirectangular')  
    # instead of bool: action='store_true': This sets the argument to True if the flag is present on the command line and False if it is omitted. This is suitable for flags that enable a feature, where the default state is False.
    parser.add_argument('--equal_names', action='store_true', 
                    help='use this flag if geotiffs and polygons have the exact same name. if not, the first string subset of the polygon name (before an underline) should be unique. e.g. Rhadamanthys_poly1_....goejson would be mapped to a geotiff that also has Rhadamanthys in it.') 
    parser.add_argument('--tesseldeg', required=False, 
                        default=2.0,
                        type=float,
                    help='degrees used for tessellation (the lonbox parameter in this script).') 


    args = parser.parse_args()

    basepath = Path(args.basepath) # titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    if args.crs_add == True:
        print('CRS ADD IS TRUE')
        print(args.crs_add)
        # otherwise, no crs_file is specified.
        equifile = args.crs_file
    savefolder = args.savefoldername
    '''
    current = os.getcwd()
    titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

    basepath =  titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/analyse_LM1_1'
    jsonpath = basepath / 'polygons/for_analysis/Equirectangular'
    equifile = 'Equirectangular_EUROPA.prj'
    savefolder = 'UFID'
    '''
    # get tessellation degrees
    lon_box = args.tesseldeg # e.g. 2 degrees longitude splits
    # equirectangular
    geo_orig = basepath / 'geotiffs/for_analysis'
    poly_orig = basepath / 'polygons/for_analysis'
    geopaths, polypaths = get_polygeopaths(poly_orig, geo_orig)
    # each geopath has one polypath (the polypath are all different, while the geopaths can be identical (they have the same base image))

    # where to find equifile
    jsonpath = poly_orig / savefolder
    #
    '''
        We loop through all polypaths in poly_orig
        For each, we identify the correct geotiff to read in. We load the geotiff array and define the longitude extent. The longitude box of 2° could be made an input parameter in a later version.
        The overlap is set to 0, because for analysis, we do not want any overlap.
        With gdal_translate, we cut the geotiff array into boxes.
        With ogr2ogr -clip, we 'cut', or rather clip, the feature layers into the same boxes and save them with the correct coordinate system, if specified as an input. (for LineaMapper output, this is already specified correctly.)
    '''
    # fixed bug with:
    # polypath =  Path('z:/Groups/PIG/Caroline/publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/analyse_LM1_1/polygons/for_analysis/C0449974429R_5_202407060352.geojson')
    for pidx, polypath in enumerate(polypaths): 
        # note: pidx is no longer used
        # define the filename:
        filename = polypath.stem
        # define a savepath:
        savepath_geotiff = geo_orig / savefolder / 'tessellation' / polypath.stem
        clippedpolytiff_savepath = geo_orig / 'clipped'
        os.makedirs(savepath_geotiff, exist_ok=True)
        os.makedirs(clippedpolytiff_savepath, exist_ok=True)
        savepath_featurelayer_geojson = poly_orig / savefolder / 'tessellation' / polypath.stem / 'geojson'
        savepath_valid_geojson = poly_orig / savefolder / 'valid'
        savepath_featurelayer_shp = poly_orig / savefolder / 'tessellation' / polypath.stem / 'shp'
        os.makedirs(savepath_featurelayer_geojson, exist_ok=True)
        os.makedirs(savepath_featurelayer_shp, exist_ok=True)
        os.makedirs(savepath_valid_geojson, exist_ok=True)

        if args.equal_names:
            # we simply deduct the geopath from the polypath
            geotiff = geopaths[pidx]
            print(polypath.stem)
            print(geotiff)
        else:
            # then, we match the first subset of the polygon to the geotiff path
            for gidx, geopath in enumerate(geopaths):
                # this looks for the original full geotiff
                if polypath.stem.split('_')[0] in str(geopath):
                    geotiff = geopath
            print(polypath.stem, geotiff.stem)


        # load the array
        arr, dataset = geotiff_to_arr(geotiff)

        # first: get area weight for subimg and store it
        # How much of the area is covered by the current polypath? --> clip raster by polygon, read in later in Merge...py
        # clip current polypath to identified geotiff
        polyfile = ogr.Open(polypath.as_posix())
        layer = polyfile.GetLayer()
        min_x, max_x, min_y, max_y = layer.GetExtent() # e.g. (230067.951154, 483692.433656, -212937.187313, -23388.202028)
        # Clip the raster
        dt = gdal.Translate(clippedpolytiff_savepath.joinpath(filename + '.tif') .as_posix(), geotiff.as_posix(), projWin=[min_x, max_y, max_x, min_y])
        # Close the layer dataset
        polyfile = None
        dt.FlushCache()
        dt = None

        # calculation of cut_size: (not straightforward, because equirectangular projection, which is increasingly distorted towards the poles, will be converted to orthographic)
        # BUT: for longitude, it must always be 2°, so I could just always take 2° boxes!
        px_width = arr.shape[1] # e.g. (568, 253), then 568 is latitude, and 253 is longitude! This is where it is mixed.
        # lon_box = 2 # 2 degrees longitude splits --> is now a parameter
        ulx, xres, xskew, uly, yskew, yres, lrx, lry = get_geotransform(dataset)
        lon_extent = abs( coordm_to_lon(ulx) - coordm_to_lon(lrx))
        cut_size = int(px_width * lon_box/lon_extent) # this is the cut size in longitudinal direction, but this is also ok for latitudinal, since this is an equirectangular grid (!) and we keep lon_box constant at 2°.
        # you could make it more complex, if you would calculate depending on tesselation grid (Moseley et al. 2020)
        # I quickly checked in Moseley et al. (2020) and will use the following tessellation windows for given latitudes:
        # until 60deg: 2deg
        # above 60deg: 4deg
        # above 70deg: 8deg
        # above 80 deg: 15 deg



        overlap = 0 # no overlap wanted!
        # new: split array into subimages into 2° longitude boxes (TODO: and stitch output together in the end)
        if np.max(arr.shape) > cut_size:
            print('I have to cut the input image into sub-images.')
            # with the row_step, we can ensure a small overlap, which is important for connecting the features
            # I adapt the moving_window_tiling_ref in a way that allows for the last row and column to not repeat too much. I can simply cut them to ensure an overlap
            row_step_subi = cut_size - overlap
            col_step_subi = cut_size - overlap
            subimgs, sub_positions = moving_window_tiling_ref(arr, row_step_subi, col_step_subi, cut_size, sort_zero=False, pad=True)

            # we produce a valid geopackage.
            # NOTE: this step is only needed once for each filename, therefore, I moved this outside the subimages loop
            # EDIT on 2024-02-21: I added -makevalid to make polygons valid. 
            # EDIT on 2024-02-22: I realised that I was not taking the manipulated (generated with add_fixed_ID_and_numSubPolygons.py) geojson files in 'Equirectangular' with inVectorPath = polypath.as_posix(). I changed this now.
            inVectorPath = jsonpath.joinpath(filename + '.geojson').as_posix() # input vector
            validVectorPath = savepath_valid_geojson.joinpath(filename + '.gpkg').as_posix()
            # first step: make this valid! and transform into geopackage
            valgeojson = 'ogr2ogr -nlt MULTIPOLYGON --config PG_USE_COPY YES -overwrite ' + validVectorPath + ' ' + inVectorPath + ' -makevalid'
            os.system(valgeojson)
            print('Geojson is now valid')

            for sidx, (subimg, subpos) in enumerate(zip(subimgs, sub_positions)):
                # loop through subimages
                
                # if sidx == 2:
                #     return # for development/debugging
                appimg = subimg
                apppos = subpos
                sidx = sidx

                # save as geotiff:
                # Actually! The subimgs I don't need. I only need the subpositions and the input image
                # https://gis.stackexchange.com/questions/65840/subsetting-geotiff-with-python 
                inDS = geotiff.as_posix() # input raster
                outDS = savepath_geotiff.joinpath(filename + '_{}.tif'.format(str(sidx))) # output raster
                print('outDS: {}'.format(outDS))
                print('\n\n')
                xoff = apppos[1]
                yoff = apppos[0]
                xsize = cut_size # arr.shape[1]
                ysize = cut_size # arr.shape[0]
                # xoff yoff xsize ysize
                translate = 'gdal_translate -srcwin %s %s %s %s -co NUM_THREADS=8 --config GDAL_CACHEMAX 512 %s %s' %(xoff, yoff, xsize, ysize, inDS, outDS.as_posix())
                os.system(translate)

                # read the new geotiff in, take coordinates for cutting the FEATURE layer
                dataset_cut = gdal.Open(outDS.as_posix(), gdal.GA_ReadOnly)
                arr_cut = dataset_cut.ReadAsArray()
                # check here if this is an empty array
                if arr_cut.max() < 0.01:
                    # print('deleting')
                    del dataset_cut, arr_cut
                    # then, delete the geotiff and continue with the next tile, without generating a geojson/shapefile
                    outDS.unlink()
                    continue
                # get geotransform
                ulx, xres, xskew, uly, yskew, yres, lrx, lry = get_geotransform(dataset_cut)
                print(ulx)
                # check if it is empty! There are some border empty tiles!
                
                print('I am processing sub-image Nr. {} out of {}.'.format(sidx, len(subimgs)))

                # clip FEATURE layer to raster extent: 
                # 'ogr2ogr ' + outVectorPath + ' ' + inVectorPath + ' -clipsrc ' + extent + ' ' + ' -skipfailures '
                # https://gis.stackexchange.com/questions/389083/clipping-shapefile-to-raster-extent-using-gdal-with-python 

                outVectorPath = savepath_featurelayer_geojson.joinpath(filename + '_{}.geojson'.format(str(sidx))).as_posix() # output 
                extent = '{} {} {} {}'.format(ulx, lry, lrx, uly)
                # also, we want to define an output coordinate system, with -a_srs (without reprojecting, unlike -t_srs)
                if args.crs_add:
                    print('Yes, args.crs_add is true.')
                    clip = 'ogr2ogr -nlt MULTIPOLYGON --config PG_USE_COPY YES -overwrite ' + outVectorPath + ' ' + validVectorPath + ' -clipsrc ' + extent + ' ' + ' -makevalid -skipfailures -a_srs ' + jsonpath.joinpath(equifile).as_posix()
                else: # no crs is added
                    clip = 'ogr2ogr -nlt MULTIPOLYGON --config PG_USE_COPY YES -overwrite ' + outVectorPath + ' ' + validVectorPath + ' -clipsrc ' + extent + ' ' + ' -makevalid -skipfailures'
                os.system(clip)
                print('Geojson is now clipped')
                # to shapefile
                shp_path = savepath_featurelayer_shp.joinpath(filename + '_{}.shp'.format(str(sidx))).as_posix() # output 
                if args.crs_add:
                    to_shapefile = 'ogr2ogr -nlt MULTIPOLYGON --config PG_USE_COPY YES -overwrite -skipfailures -makevalid -a_srs {} {} {}'.format(jsonpath.joinpath(equifile).as_posix(), shp_path, outVectorPath) # here, outvectorpath is the new input
                else: # no coordinate system is added extra.
                    to_shapefile = 'ogr2ogr -nlt MULTIPOLYGON --config PG_USE_COPY YES -overwrite -skipfailures -makevalid {} {}'.format(shp_path, outVectorPath) # here, outvectorpath is the new input

                os.system(to_shapefile)




#     # else, we have no job to do
# C0449961800R_0_202407051145.geojson C0449961800R_0_202407051145_clean.geojson
# ogr2ogr -nlt MULTIPOLYGON C0449961800R_0_202407051145_clean.geojson C0449961800R_0_202407051145.geojson -dialect sqlite -sql "select ST_MakeValid(geometry) as geometry, * from C0449961800R_0_202407051145"

#update, makevalid
# ogr2ogr -f GeoJSON  -append -update  -nln valid_now  -makevalid C0449961802R_0_202407051158.geojson C0449961802R_0_202407051158.geojson
# ogr2ogr -nlt MULTIPOLYGON C0449961802R_0_202407051158_valid.geojson C0449961802R_0_202407051158.geojson -makevalid

# ogr2ogr -f GeoJSON -dialect sqlite -sql "select ST_MakeValid(geometry),STATEFP10,ZCTA5CE10,GEOID10,CLASSFP10,MTFCC10,FUNCSTAT10,ALAND10,AWATER10,INTPTLAT10,INTPTLON10,PARTFLG10 from C0449961800R_0_202407051145.geojson C0449961800R_0_202407051145_clean.geojson
