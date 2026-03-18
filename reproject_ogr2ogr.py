
# 2024-07-12
# Caroline Haslebacher
# this script is basically a general reprojection script.
# uses ogr2ogr. mind that I found no way to add the coordinate system to the output geojson!
# quickly reproject all orthographic, valid geojson to equirectangular
# this allows them to be processed with 'Merge_polygon_calculations.py'. We need equirectangular projection there




#%%

from pathlib import Path
import os
import json
from osgeo import gdal
import numpy as np
import argparse

#%% 

if __name__ == '__main__':

    '''
    This script projects geojson from a source projection to a target projection
    either a path or a .prj file need be specified.
    in case of a path, each file needs to have the same name as the geojson

    input:
    inputpath (all geojson found within this directory are processed)
    outputpath
    src_prj_path or file
    dst_prj_path or file
    '''

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Get arguments.')

    parser.add_argument('--inputpath', required=True,
                    metavar="where to find geojson",
                    help='as posix')
    parser.add_argument('--outputpath', required=True,
                    metavar="the name of full folder to save",
                    help='as posix')  
    parser.add_argument('--src_prj', required=True,
                    metavar="source projection",
                    help='the name of the file or path to find source projection. As Pathlib Path.')  
    parser.add_argument('--dst_prj', required=True,
                    metavar="destination projection",
                    help='the name of the file or path to find destination projection. As Pathlib Path.')  

    args = parser.parse_args()

    inputpath = Path(args.inputpath)
    outputpath = Path(args.outputpath)
    src_prj = Path(args.src_prj)
    dst_prj = Path(args.dst_prj)
    '''
    test/debug with:
    current = os.getcwd()
    titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

    inputpath = titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/do_analysis/polygons/for_analysis/UFID/valid'
    outputpath = titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/do_analysis/polygons/for_analysis/Equirectangular/valid'
    src_prj = titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/do_analysis/orig_prj_files'
    dst_prj = titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/do_analysis/polygons/for_analysis/Equirectangular/Equirectangular_EUROPA.prj'
    '''
    #%%

    # make outputpath
    os.makedirs(outputpath, exist_ok=True)
    # read in geojson
    geojfiles = sorted(inputpath.glob('*.geojson'))

    # if source projection is a path, not a file, find corresponding prj files (with geojfiles stem)
    if src_prj.is_file(): 
        # then, source projection is a single file
        # make a list by repeating
        source_projs = [src_prj for x in range(len(geojfiles))]
    elif src_prj.is_dir():
        # then, it is a directory, and we get all projection files
        source_projs = []
        for geojf in geojfiles:
            # we construct projection file!
            source_projs.append(src_prj.joinpath(geojf.stem + '.prj'))


    # same for destination projection
    if dst_prj.is_file(): 
        # then, source projection is a single file
        # make a list by repeating
        dest_projs = [dst_prj for x in range(len(geojfiles))]
    elif dst_prj.is_dir():
        # then, it is a directory, and we get all projection files
        dest_projs = []
        for geojf in geojfiles:
            # we construct projection file!
            dest_projs.append(dst_prj.joinpath(geojf.stem + '.prj'))

    # the lists are all sorted, so we loop through with zip
    for geojf, src_crs, dst_crs in zip(geojfiles, source_projs, dest_projs):
        # construct output path 
        outgeo = outputpath.joinpath(geojf.name)

        # write command
        command = 'ogr2ogr -nlt MULTIPOLYGON -overwrite -skipfailures -makevalid -s_srs {} -t_srs {} {} {}'.format(src_crs.as_posix(), dst_crs.as_posix(), outgeo.as_posix(), geojf.as_posix())
        # print command
        print(command)
        # exectue command
        os.system(command)   
