# python script to reproject .tif files and .geojson files to orthographic projection
# Caroline Haslebacher
# 2023-12-05

#%%

from pathlib import Path
import os
import json
from osgeo import gdal
import numpy as np
import argparse

#%% define basepath of data
# current = os.getcwd()
# titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

# basepath = titaniach / 'lineament_detection/RegionalMaps_CH_NT_EJL_LP/extraction_from_LineaMapper/data'



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
    parser.add_argument('--savefoldername', required=False,
                        default='Equirectangular',
                    metavar="the name of sub-folder to save",
                    help='e.g. Equirectangular')  
    args = parser.parse_args()

    basepath = Path(args.basepath) # e.g. basepath = titaniach / 'lineament_detection/RegionalMaps_CH_NT_EJL_LP/extraction_from_LineaMapper/data'
    savefolder = args.savefoldername
    # e.g.
    # jsonpath = basepath / 'polygons/for_analysis/Equirectangular'

    # equirectangular
    geo_orig = basepath / 'geotiffs/for_analysis'
    poly_orig = basepath / 'polygons/for_analysis'

    #%% gdal commands for .tif files

    # asign path of Geotiff files
    # tifpath = geo_orig

    # get paths where tessellated data live
    tessell_geotiff = geo_orig / savefolder / 'tessellation'
    tessell_featurelayer_geojson = poly_orig / savefolder / 'tessellation'

    # get the geotiffs in the subdirs (polypath.stem for each observation)
    geotiffs = sorted(tessell_geotiff.glob('*/*.tif'))
    # equipath = basepath / 'polygons/for_analysis/Equirectangular'
    orthopath = basepath / 'geotiffs/for_analysis/Orthographic'

    #%% gdal commands for .tif files

    # e.g.:
    # # Geotiff files
    # tifpath = basepath / 'geotiffs/for_analysis'
    # # get paths where tessellated data live$
    # geo_orig = basepath / 'geotiffs/for_analysis'
    # poly_orig = basepath / 'polygons/for_analysis'
    # tessell_geotiff = geo_orig / 'cuts/tessellation'
    # tessell_featurelayer_geojson = poly_orig / 'cuts/tessellation'
    # # get the geotiffs in the subdirs (polypath.stem for each observation)
    # geotiffs = sorted(tessell_geotiff.glob('*/*.tif'))
    # orthopath = basepath / 'geotiffs/for_analysis/Orthographic'

    #%%
    # find all geojson files
    # tessel_geojsonpath = basepath / 'polygons/for_analysis/cuts/tessellation' --> outdated
    # tessell_featurelayer_geojson = poly_orig / savefolder / 'tessellation'

    #%% EDIT on 2023-02-13: projection must be changed so that the center longitude of the observation is equal to the central meridian

    # define where to save projection files
    ortho_prj_path = basepath / 'orthographic_projection'
    os.makedirs(ortho_prj_path, exist_ok=True)

    ORIG_prj_path = basepath / 'orig_prj_files'
    os.makedirs(ORIG_prj_path, exist_ok=True) 
    # get list of original geojson files   
    geojfiles_orig = sorted(poly_orig.glob('*.geojson'))
    # get coordinate system from geojson, save to .prj file (otherwise this lead to errors of 'missing [')
    for geojfile in geojfiles_orig:
        print(geojfile)
        with open(geojfile) as jsf_orig:
            geoms_orig = json.load(jsf_orig)
        source_crs = geoms_orig['crs']['properties']['name']

        centre_lon = source_crs.split('"central_meridian",')[1].split(']')[0] # center longitude of projection
        centre_lat = source_crs.split('"latitude_of_origin",')[1].split(']')[0] # central latitude, e.g. -11.4713, split from the proj string
        # e.g 'PROJCS["Orthographic_EUROPA_template_2",GEOGCS["New Geographic Coordinate System",DATUM["Custom",SPHEROID["<Custom>",1560800,0]],PRIMEM["Reference_Meridian",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]]],PROJECTION["Orthographic"],PARAMETER["latitude_of_origin",-11.4713],PARAMETER["central_meridian",133.1926],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
        # print('LAT: {}'.format(centre_lat))
        # print('LON: {}\n'.format(centre_lon))
        # define projection file (copied from S:\Groups\PIG\Caroline\lineament_detection\galileo_manual_segmentation\data\polygons\for_analysis\Sinusoidal\Sinusoidal_EUROPA_0.prj)
        source_crs_newstring = f'PROJCS["Orthographic_EUROPA_template_2",GEOGCS["New Geographic Coordinate System",DATUM["<Custom>",SPHEROID["<Custom>",1560800.0,0.0]],PRIMEM["Reference_Meridian",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Orthographic"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Longitude_Of_Center",{centre_lon}],PARAMETER["Latitude_Of_Center",{centre_lat}],UNIT["Meter",1.0]]'
        # write to file
        # write to file
        with open(ORIG_prj_path.joinpath("{}.prj".format(geojfile.stem)), "w") as prf_file:
            prf_file.write(source_crs_newstring)        

    #%%

    # from orthographic!
    def ortho_coordm_to_lon_lat(x, y, lambda_0, phi_1): 
        '''
        input: 
        - longitude (x) in meters, from Orthographic projection
        - latitude (y) in meters, from Orthographic projection
        - lambda_0: center longitude of projection, in radians
        - phi_1: center latitude of projection, in radians
        output:
        - longitude and latitude in degrees
        fixed:
        - radius
        Projection information from https://neacsu.net/docs/geodesy/snyder/5-azimuthal/sect_20/ , equations (20-14, 20-15, 20-18, 20-19)
        '''
        radius = 1560800 # metres, assuming sphere
        rho = np.sqrt( x**2 + y**2 ) # = radius * sin(c)
        c = np.arcsin( rho/radius) # angular distance from center of projection

        phi_lat = np.arcsin( np.cos(c)*np.sin(phi_1) + ( y*np.sin(c)*np.cos(phi_1)/rho ) ) # latitude
        lambda_lon = lambda_0 + np.arctan( x*np.sin(c) / ( rho*np.cos(phi_1)*np.cos(c) - y*np.sin(phi_1)*np.sin(c) ) ) # longitude

        return lambda_lon/np.pi*180, phi_lat/np.pi*180 # lon, lat

    # project the geotiff (orthographic to orthographic)
    for geosin in geotiffs:
        # print(geosin)
        dataset = gdal.Open(geosin.as_posix(), gdal.GA_ReadOnly)
        ulx, xres, xskew, uly, yskew, yres  = dataset.GetGeoTransform()
        # print(ulx)
        # get central coordinates of dataset
        lambda_0 = np.pi/180 *float(dataset.GetProjection().split('"central_meridian",')[1].split(']')[0]) # center longitude of projection
        phi_1 = np.pi/180 *float(dataset.GetProjection().split('"latitude_of_origin",')[1].split(']')[0]) # central latitude, e.g. -11.4713, split from the proj string
        # e.g 'PROJCS["Orthographic_EUROPA_template_2",GEOGCS["New Geographic Coordinate System",DATUM["Custom",SPHEROID["<Custom>",1560800,0]],PRIMEM["Reference_Meridian",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]]],PROJECTION["Orthographic"],PARAMETER["latitude_of_origin",-11.4713],PARAMETER["central_meridian",133.1926],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
        lono, lato = ortho_coordm_to_lon_lat(ulx + (0.5*dataset.RasterXSize * xres), uly + (0.5*dataset.RasterYSize * yres), lambda_0, phi_1 )
        centre_lat = '{:.4f}'.format(lato) # y is latitude
        centre_lon = '{:.4f}'.format(lono) # x is longitude
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

        # reproject to orthographic
        # example command: gdalwarp -t_srs ./Mercator/Mercator_EUROPA_0.prj input.tif output.tif
        # orthographic to orthographic
        # source_crs from file
        source_crs = ORIG_prj_path.joinpath("{}.prj".format(geosin.stem))
        # target savepath:
        orthoind_path = (orthopath / '_'.join(geosin.stem.split('_')[:-1])) # retrieve only the original stem; from '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed_0', we want '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed' for the original projection
        os.makedirs(orthoind_path, exist_ok=True)
        command = 'gdalwarp -t_srs {} {} {}'.format(ortho_prj_path.joinpath("{}.prj".format(geosin.stem)).as_posix(), geosin, orthoind_path.joinpath(geosin.stem + '.tif').as_posix()) # source_crs,  ortho_prj_path.joinpath("{}.prj".format(geosin.stem)).as_posix()
        # -s_srs {}  (lead to errors)
        print(command)
        print('\n')
        os.system(command)   


    #%%

    # find all geojson files
    # tessel_geojsonpath = basepath / 'polygons/for_analysis/cuts/tessellation' --> outdated
    # tessell_featurelayer_geojson = poly_orig / savefolder / 'tessellation'
    geojfiles = sorted(tessell_featurelayer_geojson.glob('*/geojson/*.geojson'))

    # re-deifne orthopath
    orthopath = basepath / 'polygons/for_analysis/Orthographic'

    #%% EDIT on 2023-02-13: I deleted "sinuspath.joinpath(sinusfile)"
    # reproject and save as geojson, and also convert to shapefiles

    for geojfile in geojfiles:
        # reproject to orthographic
        # example command: ogr2ogr -s_srs EPSG:4326 -t_srs EPSG:3857 output.gpkg input.gpkg

        # source_crs from file
        source_crs = ORIG_prj_path.joinpath("{}.prj".format('_'.join(geojfile.stem.split('_')[:-1]))) # e.g. C0449961800R_0_202407051145_12.geojson
        # print(source_crs)
        # make orthographic paths:
        orthoind_path_geojson = (orthopath / 'geojson' / '_'.join(geojfile.stem.split('_')[:-1])) # retrieve only the original stem for defining the directory; from '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed_0', we want '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed'
        os.makedirs(orthoind_path_geojson, exist_ok=True)
        orthoind_path_shp = (orthopath / 'shapefiles' / '_'.join(geojfile.stem.split('_')[:-1])) # retrieve only the original stem; from '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed_0', we want '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed'
        os.makedirs(orthoind_path_shp, exist_ok=True)

        # orthographic, add coordinate system here! somehow, this does not do what I want... I delete -a_srs again
        command = 'ogr2ogr -nlt MULTIPOLYGON -overwrite -skipfailures -makevalid -s_srs {} -t_srs {} {} {}'.format(source_crs, ortho_prj_path.joinpath("{}.prj".format(geojfile.stem)).as_posix(), (orthoind_path_geojson).joinpath(geojfile.stem + '.geojson').as_posix(), geojfile)
        print(command)
        os.system(command)   

        # orthographic, shapefile
        command = 'ogr2ogr -nlt MULTIPOLYGON -overwrite -skipfailures -makevalid -s_srs {} -t_srs {} {} {}'.format(source_crs, ortho_prj_path.joinpath("{}.prj".format(geojfile.stem)).as_posix(), (orthoind_path_shp).joinpath(geojfile.stem + '.shp').as_posix(), geojfile)
        print(command)
        os.system(command)   


    #%%
