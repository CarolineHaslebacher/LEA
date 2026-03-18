# Caroline Haslebacher
# 2024-02-20
# this script runs all scripts necessary for the extraction algorithm
# rename to MAIN_loop1_to_loop4.py for functionality
# original name: MAIN_loop1_to_loop4_GalileoDTMs.py

# how to run:
# cd into/the/current/working/directory
# python MAIN_loop1_to_loop4.py

#%% import

import os
import subprocess
import glob
from pathlib import Path
# from osgeo import gdal, gdalnumeric, ogr, osr
# from datetime import datetime
# import time


#%% define path to titania

current = Path.cwd()
titaniach = current.parents[0]

#%%

# we define here the function get_polygeopaths, because it might be specific to a certain dataset
def get_polygeopaths(polypath_orig, geopath_orig):
    # take code snippets from 'preparing_Galileo_tiles_from_qgis.py'
    # get geojson and geotiff paths

    # os.makedirs(basepath / 'data/polygons/for_analysis', exist_ok=True)
    polypaths = sorted((polypath_orig).glob('*.geojson'))

    # construct geotiff from polypath
    geopaths = []
    for polypath in polypaths:
        # For this to work, the geotiff and the polygon need to have the exact same name, except for the line ending
        geopath = geopath_orig.joinpath('{}.tif'.format(polypath.stem))

        geopaths.append(geopath)

    return geopaths, polypaths

#%%
if __name__ == '__main__':

    '''
    basepath: where to find /data/geotiffs/for_analysis and /data/polygons/for_analysis
    crs_file: the file 'Equirectangular_EUROPA.prj' must be copied to jsonpath
    jsonpath: for example 'polygons/for_analysis/Equirectangular'
    '''
    basepath = titaniach # WindowsPath('c:/Users/chaslebacher/OneDrive - SWRI/working/Galileo_DTM')
    # workaround on Laptop for openMP library issue
    # os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    
    # loop1: add attributes unique, fixed ID (UFID) and number of subpolygons to the geojson files
    command = 'python add_fixed_ID_and_numSubPolygons.py --basepath="{}" --savefoldername={}'.format((basepath / 'data').as_posix(), 'Equirectangular')
    print(command)
    # os.system(command)
    subprocess.run(command)

    # loop2: 
    # 2.1 cut polygon and geotiffs into smaller pieces
    command = 'python Cut_geotiff_to_smaller_geotiffs.py --basepath="{}"  --savefoldername={} --crs_add=True --crs_file={} --equal_names --tesseldeg=0.3'.format((basepath / 'data').as_posix(), 'Equirectangular', 'Equirectangular.prj')
    print(command)
    # os.system(command)
    subprocess.run(command)
    # 2.2 project to orthographic projection 
    command = 'python gdal_reproject_orthographic.py --basepath="{}" --savefoldername={} --crs_file={}'.format((basepath / 'data').as_posix(), 'Equirectangular', 'Equirectangular.prj')
    print(command)
    # os.system(command)
    subprocess.run(command)

    # loop3: extract azimuth, perhaps cut into smaller segments for better azimuth analysis, extract width, length, area, etc.
    command = 'python science_with_manual_segs_algo.py --basepath="{}" --crsfoldername={} --crs_file={} --output_foldername={} --SSI_id=true --illumination_info'.format((basepath / 'data').as_posix(), 'Equirectangular', 'Equirectangular.prj', 'output/extraction_loop3_tessel03_AJik')
    print(command)
    # os.system(command)
    subprocess.run(command)

    # loop4: merging 
    command = 'python Merge_polygon_calculations.py  --basepath="{}" --input_foldername={} --output_foldername={} --crsfoldername={} --crs_file={} --illumination_info'.format((basepath / 'data').as_posix(), 'output/extraction_loop3_tessel03_AJik', 'output/extraction_loop4_tessel03_AJik', 'Equirectangular', 'Equirectangular.prj')
    print(command)
    # os.system(command)
    subprocess.run(command)

    # and perhaps prepare input for vMF
    command = 'python Generate_input_for_vMF_SUNAZI.py  --basepath="{}" --input_foldername={} --output_foldername={}'.format((basepath / 'data').as_posix(), 'output/extraction_loop4_tessel03_AJik', 'output/extraction_loop4_tessel03_AJik/for_vMF')
    print(command)
    # os.system(command)
    subprocess.run(command)

# %%
