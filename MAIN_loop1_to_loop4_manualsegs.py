# Caroline Haslebacher
# 2024-02-20
# this script runs all scripts necessary for the extraction algorithm
# rename to MAIN_loop1_to_loop4.py for functionality
# original name: MAIN_loop1_to_loop4_manualsegs.py

# how to run:
# cd into/the/current/working/directory
# python MAIN_loop1_to_loop4.py

#%% import

import os
import glob
from pathlib import Path
# from datetime import datetime
# import time


#%% define path to titania

current = os.getcwd()
titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

#%%

# we define here the function get_polygeopaths, because it might be specific to a certain dataset
def get_polygeopaths(polypath_orig, geopath_orig):
    # take code snippets from 'preparing_Galileo_tiles_from_qgis.py'
    # get geojson and geotiff paths

    # os.makedirs(basepath / 'data/polygons/for_analysis', exist_ok=True)
    polypaths = sorted((polypath_orig).glob('*.geojson'))

    # # crs dict
    # crsd = {'Equi-cog': 'Europa_Mosaics_Equirectangular', 'NPola-cog': 'Europa_Mosaics_NPolar', 'SPola-cog': 'Europa_Mosaics_SPolar'}
    # geopaths = []

    # for polypath in polypaths:
    #     # For this to work, the geotiff and the polygon need to have the exact same name, except for the line ending
    #     geopath = geopath_orig.joinpath('{}.tif'.format(polypath.stem))

    #     geopaths.append(geopath)

    # construct geotiff from polypath
    geopaths = []
    # crs dict
    crsd = {'Equi-cog': 'Europa_Mosaics_Equirectangular', 'NPola-cog': 'Europa_Mosaics_NPolar', 'SPola-cog': 'Europa_Mosaics_SPolar'}

    for polypath in polypaths:
        # separate filename by '_'
        # e.g. filename is G7ESTYRMAC02_GalileoSSI_Equi-cog_fixed.geojson
        namesep = polypath.stem.split('_')
        # first entry in list is the ID
        id = namesep[0]
        # third entry is the coordinate system (equirectangular or polar),
        # with which we can get the folder name
        # (the folder we expect to find our geotiff)
        # crs_folder = crsd[namesep[2]]

        # name of the geotiff file (gets appended to the geopath)
        # (this is probably quite an ugly way of doing this)
        geoname = id + '_' + namesep[1] + '_' + namesep[2] + '.tif'
        # full geotiff path
        geopath = (geopath_orig).joinpath(geoname)

        geopaths.append(geopath)

    return geopaths, polypaths

#%%
if __name__ == '__main__':

    '''
    basepath: where to find /data/geotiffs/for_analysis and /data/polygons/for_analysis
    crs_file: the file 'Equirectangular_EUROPA.prj' must be copied to jsonpath
    jsonpath: for example 'polygons/for_analysis/Equirectangular'
    '''
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/' # manual segs analysis
    # basepath = titaniach / 'lineament_detection/RegionalMaps_CH_NT_EJL_LP/mapping/'

    # loop1: add attributes unique, fixed ID (UFID) and number of subpolygons to the geojson files
    command = 'python add_fixed_ID_and_numSubPolygons.py --basepath={} --savefoldername={}'.format((basepath / 'data').as_posix(), 'Equirectangular')
    print(command)
    os.system(command)

    # loop2: 
    # 2.1 cut polygon and geotiffs into smaller pieces
    command = 'python Cut_geotiff_to_smaller_geotiffs.py --basepath={}  --savefoldername={} --crs_add=True --crs_file={}'.format((basepath / 'data').as_posix(), 'Equirectangular', 'Equirectangular_EUROPA.prj')
    print(command)
    os.system(command)
    # 2.2 project to orthographic projection 
    command = 'python gdal_reproject_orthographic.py --basepath={} --savefoldername={} --crs_file={}'.format((basepath / 'data').as_posix(), 'Equirectangular', 'Equirectangular_EUROPA.prj')
    print(command)
    os.system(command)

    # loop3: extract azimuth, perhaps cut into smaller segments for better azimuth analysis, extract width, length, area, etc.
    command = 'python science_with_manual_segs_algo.py --basepath={} --crsfoldername={} --crs_file={} --output_foldername={} --illumination_info'.format((basepath / 'data').as_posix(), 'Equirectangular', 'Equirectangular_EUROPA.prj', 'output/extraction_loop3')
    print(command)
    os.system(command)

    # loop4: merging 
    command = 'python Merge_polygon_calculations.py  --basepath={} --input_foldername={} --output_foldername={} --crsfoldername={} --crs_file={} --illumination_info'.format((basepath / 'data').as_posix(), 'output/extraction_loop3', 'output/extraction_loop4', 'Equirectangular', 'Equirectangular_EUROPA.prj')
    print(command)
    os.system(command)

    # and perhaps prepare input for vMF
    command = 'python Generate_input_for_vMF.py  --basepath={} --input_foldername={} --output_foldername={}'.format((basepath / 'data').as_posix(), 'output/extraction_loop4', 'output/extraction_loop4/for_vMF')
    print(command)
    os.system(command)

# %%
