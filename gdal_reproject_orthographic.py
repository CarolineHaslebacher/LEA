# python script to reproject .tif files and .geojson files to mercator and sinusoidal projection
# Caroline Haslebacher
# 2023-12-05

#%%

from pathlib import Path
import os
import json
from osgeo import gdal
import numpy as np
import argparse

from science_with_manual_segs_algo import coordm_to_lon, coordm_to_lat

#%% define basepath of data
# current = os.getcwd()
# titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

# this got substituted! basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'




if __name__ == '__main__':

    '''
    This script projects tessellated geotiff and geojson to orthographic projection 
    '''

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Get basepath.')

    parser.add_argument('--basepath', required=True,
                    metavar="where to find /data/geotiffs/for_analysis and /data/polygons/for_analysis",
                    help='as posix')
    parser.add_argument('--crs_file', required=False,
                        default='Equirectangular_EUROPA.prj',
                    metavar="the name of the crs file. Folder is jsonpath.",
                    help='jsonpath.joinpath(equifile) has to be valid.')   
    parser.add_argument('--savefoldername', required=False,
                        default='Equirectangular',
                    metavar="the name of sub-folder to save",
                    help='e.g. Equirectangular')  
    args = parser.parse_args()

    basepath = Path(args.basepath) # titaniach / 'lineament_detection/galileo_manual_segmentation/data'

    equifile = args.crs_file
    savefolder = args.savefoldername
    # e.g.
    # basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    # jsonpath = basepath / 'polygons/for_analysis/Equirectangular'
    # equifile = 'Equirectangular_EUROPA.prj'

    # equirectangular
    geo_orig = basepath / 'geotiffs/for_analysis'
    poly_orig = basepath / 'polygons/for_analysis'

    # where to find equifile
    jsonpath = poly_orig / savefolder

    #%% gdal commands for .tif files

    # asign path of Geotiff files
    tifpath = geo_orig

    # get paths where tessellated data live
    tessell_geotiff = geo_orig / savefolder / 'tessellation'
    tessell_featurelayer_geojson = poly_orig / savefolder / 'tessellation'

    # get the geotiffs in the subdirs (polypath.stem for each observation)
    geotiffs = sorted(tessell_geotiff.glob('*/*.tif'))
    # equipath = basepath / 'polygons/for_analysis/Equirectangular'
    orthopath = basepath / 'geotiffs/for_analysis/Orthographic'

    # filenames of projection files
    # equifile = 'Equirectangular_EUROPA.prj' # this must be located in the equipath directory


    #%%

    # define where to save projection files
    ortho_prj_path = basepath / 'orthographic_projection'
    os.makedirs(ortho_prj_path, exist_ok=True)

    for geosin in geotiffs:
        # print(geosin)
        dataset = gdal.Open(geosin.as_posix(), gdal.GA_ReadOnly)
        ulx, xres, xskew, uly, yskew, yres  = dataset.GetGeoTransform()
        # print(dataset.GetProjection())
        # print(dataset.GetGeoTransform())
        # get central coordinates of dataset
        # central_merid = '{:.3f}'.format(360-coordm_to_lon(ulx + (0.5*dataset.RasterXSize * xres))) # x is longitude
        centre_lat = '{:.4f}'.format(coordm_to_lat(uly + (0.5*dataset.RasterYSize * yres))) #  y is latitude
        centre_lon = '{:.4f}'.format(360-coordm_to_lon(ulx + (0.5*dataset.RasterXSize * xres))) # x is longitude
        # print(centre_lat)
        # print(centre_lon)
        # define projection file (copied from S:\Groups\PIG\Caroline\lineament_detection\galileo_manual_segmentation\data\polygons\for_analysis\Sinusoidal\Sinusoidal_EUROPA_0.prj)
        orthographic = f'PROJCS["Orthographic_EUROPA_template_2",GEOGCS["New Geographic Coordinate System",DATUM["<Custom>",SPHEROID["<Custom>",1560800.0,0.0]],PRIMEM["Reference_Meridian",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Orthographic"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Longitude_Of_Center",{centre_lon}],PARAMETER["Latitude_Of_Center",{centre_lat}],UNIT["Meter",1.0]]'
        # write to file
        print(orthographic)
        with open(ortho_prj_path.joinpath("{}.prj".format(geosin.stem)), "w") as prf_file:
            prf_file.write(orthographic)
        # close dataset
        del dataset

    #%% EDIT on 2023-02-13: I deleted "sinuspath.joinpath(sinusfile).as_posix()", and had to add -s_srs {} 
    #  sinuspath.joinpath("{}.prj".format(geotif.stem))

    for geotif in geotiffs:
        # reproject to mercator and orthographic
        # example command: gdalwarp -t_srs ./Mercator/Mercator_EUROPA_0.prj input.tif output.tif

        # # Mercator
        # command = 'gdalwarp -s_srs {} -t_srs {} {} {}'.format(jsonpath.joinpath(equifile).as_posix(), mercatorpath.joinpath(mercatorfile).as_posix(), geotif, mercatorpath.joinpath(geotif.stem + '.tif').as_posix())
        # print(command)
        # os.system(command)   

        # orthographic
        orthoind_path = (orthopath / '_'.join(geotif.stem.split('_')[:-1])) # retrieve only the original stem; from '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed_0', we want '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed'
        os.makedirs(orthoind_path, exist_ok=True)
        command = 'gdalwarp -s_srs {} -t_srs {} {} {}'.format(jsonpath.joinpath(equifile).as_posix(), ortho_prj_path.joinpath("{}.prj".format(geotif.stem)).as_posix(), geotif, orthoind_path.joinpath(geotif.stem + '.tif').as_posix())
        print(command)
        os.system(command)   


    #%%

    # find all geojson files
    tessel_geojsonpath = jsonpath / 'tessellation'
    geojfiles = sorted(tessel_geojsonpath.glob('*/geojson/*.geojson'))

    orthopath = poly_orig / 'Orthographic'
    # equipath = basepath / 'polygons/for_analysis/Equirectangular'

    # reproject and save as geojson, and also convert to shapefiles
    for geojfile in geojfiles:
        # reproject to mercator and sinusoidal
        # example command: ogr2ogr -s_srs EPSG:4326 -t_srs EPSG:3857 output.gpkg input.gpkg

        # make orthographic paths:
        orthoind_path_geojson = (orthopath / 'geojson' / '_'.join(geojfile.stem.split('_')[:-1])) # retrieve only the original stem for defining the directory; from '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed_0', we want '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed'
        os.makedirs(orthoind_path_geojson, exist_ok=True)
        orthoind_path_shp = (orthopath / 'shapefiles' / '_'.join(geojfile.stem.split('_')[:-1])) # retrieve only the original stem; from '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed_0', we want '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed'
        os.makedirs(orthoind_path_shp, exist_ok=True)

        # orthographic
        command = 'ogr2ogr -nlt MULTIPOLYGON -overwrite -skipfailures -makevalid -s_srs {} -t_srs {} {} {}'.format(jsonpath.joinpath(equifile).as_posix(), ortho_prj_path.joinpath("{}.prj".format(geojfile.stem)).as_posix(), (orthoind_path_geojson).joinpath(geojfile.stem + '.geojson').as_posix(), geojfile)
        print(command)
        os.system(command)   

        # orthographic, shapefile
        command = 'ogr2ogr -nlt MULTIPOLYGON -overwrite -skipfailures -makevalid -s_srs {} -t_srs {} {} {}'.format(jsonpath.joinpath(equifile).as_posix(), ortho_prj_path.joinpath("{}.prj".format(geojfile.stem)).as_posix(), (orthoind_path_shp).joinpath(geojfile.stem + '.shp').as_posix(), geojfile)
        print(command)
        os.system(command)   


#%%




