# 2024-02-27
# Caroline Haslebacher
# this script generates input for the R-Studio code written by Lutz Dümbgen for regression of directional data

#%% import
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import numpy as np
from osgeo import gdal
from matplotlib.colors import to_hex
import time
import argparse

from MAIN_loop1_to_loop4 import get_polygeopaths #, get_geotransform, load_tiff_arr, gaussian_kernel
# from lineana_funcpool import get_colors_and_markers, get_resolution, load_df, read_list

#%%
# current = os.getcwd()
# titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

#%%

def load_df_polystem(polystem, input_path):
    df = pd.read_csv((input_path / 'individual_polystems').joinpath('{}_parents_lineaments.csv'.format(polystem)))
    return df

def load_df_ALL(input_path):
    path = input_path.joinpath('ALL_parents_lineaments.csv')
    df = pd.read_csv(path) #, skip_blank_lines=True, na_values=np.nan, converters={"profiles": read_list(), "stdprofs": read_list(), "bprofiles": read_list(), "bstdprofiles": read_list()})
    return df

def load_df_polystem_CHILDREN(polystem, input_path):
    df = pd.read_csv((input_path / 'individual_polystems_children').joinpath('{}_children_lineaments.csv'.format(polystem)))
    return df

def load_df_ALL_CHILDREN(input_path):
    path = input_path.joinpath('ALL_children_lineaments.csv')
    df = pd.read_csv(path) #, skip_blank_lines=True, na_values=np.nan, converters={"profiles": read_list(), "stdprofs": read_list(), "bprofiles": read_list(), "bstdprofiles": read_list()})
    return df

#%%

def project_north(theta, phi):
    '''
    lon (theta), lat (phi)
    this function projects (0,0,1) onto the tangential plane defined by lon (theta), lat (phi)
     we use

    lat = np.pi/180 * (-1*sel_df['lat_d'])
    X1 = np.cos(lat)*np.cos(lon)
    X2 = np.cos(lat)*np.sin(lon)
    X3 = -np.sin(lat)

    '''

    N1 = -np.sin(phi)*np.cos(theta)
    N2 = -np.sin(phi)*np.sin(theta)
    N3 = np.cos(phi)

    return np.array([N1, N2, N3])

def construct_east(theta, phi):
    '''
    lon (theta), lat (phi)
    this function constructs the east direction by using the normal vector of the tangent and the north projected vector.
    Calculating the scalar product of North x normal lets us express the east vector with.

    '''

    E1 = -np.sin(theta)
    E2 = np.cos(theta)
    E3 = 0 * theta # we multiply with theta in case it is an array as input

    return np.array([E1, E2, E3])

def alpha_to_vec(alpha, lon, lat):
    '''
    input: angle alpha defined in clockwise direction from North
    output: a normalized vector of the direction that corresponds to angle alpha in the tangential plane 

    alpha = 33
    lon = np.pi/180 *40
    lat = np.pi/180 *(-20)
    '''

    vec = np.cos(alpha)*project_north(lon, lat) + np.sin(alpha)*construct_east(lon, lat) # shape (3, N)
    # shape (3,)
    # normalize by dividing by its length (Betrag) (absolute length --> absl)
    absl = np.sqrt(np.sum((vec**2), axis=0)) # shape (N,)

    return vec/absl # output shape (3, N)

#%% I select a suitable test region
# I take the 'children' data, since I believe for azimuth analysis, this is more accurate.

if __name__ == '__main__': # this allows me to import functions defined above!

    '''
    This script converts the output azimuth, location, and time into vectors on the sphere
    and into coordinates on a 2D map.
    '''
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Get basepath.')

    parser.add_argument('--basepath', required=True,
                    metavar="where to find /data/geotiffs/for_analysis and /data/polygons/for_analysis",
                    help='as posix')
    parser.add_argument('--input_foldername', required=False,
                        default='output/extraction_loop4_csv',
                    metavar="the name of sub-folder to read",
                    help='e.g. output/extraction_loop4_csv') 
    parser.add_argument('--output_foldername', required=False,
                        default='output/extraction_loop4_csv/for_vMF',
                    metavar="the name of sub-folder to save",
                    help='e.g. output/extraction_loop4_csv/for_vMF')  


    args = parser.parse_args()

    basepath = Path(args.basepath) # titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    input_folder = args.input_foldername
    output_folder = args.output_foldername
    # e.g.
    # basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    '''
    current = os.getcwd()
    titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

    basepath = titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/do_analysis'
    # azimuth:
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/azimuth_analysis/LineaMapper_data/data/'
    input_folder = 'output/extraction_loop4'
    output_folder = 'output/extraction_loop4/for_vMF'

    ## manualsegs:
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/TEST_extraction_workflow/data'
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    input_folder = 'output/extraction_loop4'
    output_folder = 'output/extraction_loop4/for_vMF'  

    # testset 17ESREG 0 < azimuth < 30 deg
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham/data/tests/auxiliary_files'
    input_folder = 'auxiliary_files'
    output_folder = 'auxiliary_files/for_vMF'  

    '''
    # equirectangular
    geo_orig = basepath / 'geotiffs/for_analysis'
    poly_orig = basepath / 'polygons/for_analysis'
    #####

    geopaths, polypaths = get_polygeopaths(poly_orig, geo_orig)

    savepath_csv = basepath.parents[0] / output_folder # e.g. titaniach / 'lineament_detection/galileo_manual_segmentation/output/extraction_loop4_csv/for_vMF'
    os.makedirs(savepath_csv, exist_ok=True) 

    input_path = basepath.parents[0] / input_folder

    #%%
    # first, load the full dataframe (dff) to determine the maximum FPK value, so that we can adjust the bins, for a fair comparison of 'ages'

    test = False

    if test:

        # TEST SET = 17ESREGMAP02 TRAINING REGION
        trainingregionpath = {'path': polypaths[8], 'id_int_selection': [2]}
        # elevenESREG:
        elevenESREGpath = {'path': polypaths[1], 'id_int_selection': [2,4]}
        selectedpaths = [trainingregionpath, elevenESREGpath]

        for seldict in selectedpaths:
            selpath = seldict['path']
            print(selpath)
            # load dataframe for conamara
            df = load_df_polystem_CHILDREN(selpath.stem.split('_')[0], input_path)
            # filter for double ridges
            sel_df_list = []
            for id_int in seldict['id_int_selection']:
                sel_df_list.append(df[df['id_int'] == id_int])
            # concat dataframe
            sel_df = pd.concat(sel_df_list)

            # discard FPK = -1
            sel_df = sel_df[sel_df['FPK'] >= 0]

            # construct Y1, Y2, X1, X2, Time dataframe
            # Y1 = cos(2 alpha)
            Y1 = np.cos(2* np.pi/180* sel_df['azimuth'])
            # Y2 = sin(2 alpha)
            Y2 = np.sin(2* np.pi/180* sel_df['azimuth'])
            # X1 = longitude in radians
            X1 = np.pi/180 * sel_df['lon_d']
            # X2 = latitude in radians
            X2 = np.pi/180 * sel_df['lat_d']
            # 'Time'
            T = sel_df['FPK']

            df_traininreg_vMF = pd.DataFrame([Y1, Y2, X1, X2, T]).T
            df_traininreg_vMF.columns = ['Y1', 'Y2', 'X1', 'X2', 'T']

            # save
            df_traininreg_vMF.to_csv(savepath_csv.joinpath('{}_{}_test.csv'.format(selpath.stem.split('_')[0], seldict['id_int_selection'])))




    # %% full dataframe
    # run on 2024-03-04


    # make sure the paths are correct!
    # print(polypaths)

    dff_parents = load_df_ALL(input_path)
    dff_children = load_df_ALL_CHILDREN(input_path)
    # test: dff_parents['group_index'] = dff_parents.groupby(['emission']).ngroup()

    # these two lines were only needed for the testset 17ESREG 0 < azimuth < 30 deg
    # dff_children_list = [x[1] for x in dff_children.iterrows() if x[1]['UFID'] in list(dff_parents['UFID'])]
    # dff_children = pd.DataFrame(dff_children_list)

#%%
    # add length property to children:
    dff_children['parent_length'] = [dff_parents[dff_parents['UFID'] == x[1]['UFID']]['length_km'].values[0] for x in dff_children.iterrows()] # [dff_parents[dff_parents[x['UFID']]] for x in dff_children]
    # e.g. dff_parents[dff_parents['UFID'] == '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed_10']['length_km'].values[0] returns 50.52
    # we save data for:
    # - ridge complexes and bands, curved and straight
    # - double ridges and undiff lineae, curved and straight

    dfs = []

    # custom design for azimuth thesis draft
    for dff, dffname, youngold in [(dff_parents, 'parents', False), (dff_children, 'children', True)]:
        # create group index
        dff['cluster'] = dff.groupby(['emission']).ngroup()

        bands_df = dff[dff['id_int']==1]
        doublridge_df = dff[dff['id_int']==2]
        ridgecomplex_df = dff[dff['id_int']==3]
        undifflinea_df = dff[dff['id_int']==4]
        
        # append tuples (df, dfname, FPKFilter)$
        # NOTE: FPK filter is false for these.
        dfs.append( (bands_df, f'Band_{dffname}', False))
        dfs.append( (ridgecomplex_df, f'RC_{dffname}', False))
        dfs.append( (doublridge_df, f'DR_{dffname}', False))
        dfs.append( (undifflinea_df, f'UL_{dffname}', False))

    # technically, we do now everything without an FPK filter. But it's here if needed
    for df_idx, (sel_df, df_name, FPKfilter) in enumerate(dfs):
        # # discard FPK = -1 --> this filters for length
        # if FPKfilter:
        #     # only filter if FPKfilter is True
        #     sel_df = sel_df[sel_df['FPK'] >= 0]

        # no 2D datasets needed, therefore some lines are commented out
        # transformation needed:
        # so far, the coordinate system was defined so that North is at 0° (on the y-axis)
        # and the positive direction is clockwise, so that east is at 90° and south at 180°
        # however, we now want to transform coordinates so that they are in line with a unit circle coord system (0° is east, 90° is North (anti-clockwise))
        azimuth_transf = 90-sel_df['azimuth']
        # furthermore, we want the longitude to be east-positive instead of west-positive:
        lon_transformed = 360 - sel_df['lon_d']
        # # construct Y1, Y2, X1, X2, Time dataframe
        # # Y1 = cos(2 alpha)
        # Y1 = np.cos(2* np.pi/180* azimuth_transf)
        # # Y2 = sin(2 alpha)
        # Y2 = np.sin(2* np.pi/180* azimuth_transf)
        # # X1 = longitude in radians
        # X1 = np.pi/180 * lon_transformed
        # # X2 = latitude in radians
        # X2 = np.pi/180 * sel_df['lat_d']
        # # 'Time'
        # T = sel_df['FPK']
        # the cluster
        cluster = sel_df['cluster']

        # let's add the sun azimuth for each observation. can be used to calculate the sun azimuth for one gridpoint
        sunazi = 90 - sel_df['sunazimuth'] # transformed to mathematical coord system

        # df_traininreg_vMF = pd.DataFrame([Y1, Y2, X1, X2, T, sunazi, cluster]).T
        # df_traininreg_vMF.columns = ['Y1', 'Y2', 'X1', 'X2', 'T', 'sun', 'cluster']

        # # save
        # df_traininreg_vMF.to_csv(savepath_csv.joinpath('{}_for_vMF.csv'.format(df_name)))




        ################# Spherical coordinates
        # construct (X1, X2, X3), (Y1, Y2, Y3), (V1, V2, V3)
        # longitude in radians, to east-positive!
        # latitude from -90 to 90 --> change the sign by multiplying by 
        # X1 = cos(lat)*cos(lon)
        # X2 = cos(lat)*sin(lon)
        # X3 = -sin(lat)

        lon = np.pi/180 * (360 - sel_df['lon_d']) # transform longitude to east-positive
        lat = np.pi/180 * (sel_df['lat_d'])
        X1 = np.cos(lat)*np.cos(lon)
        X2 = np.cos(lat)*np.sin(lon)
        X3 = np.sin(lat)
        # 'Time'
        # T = sel_df['FPK']

        alpha = np.pi/180* sel_df['azimuth'] # is expected in North-East (x,y) system. north is 0deg, east is 90 deg. Therefore, no transform.
        # about the shapes:
        # project_north(lon, lat).shape = (3, N)
        # alpha.shape (N,)
        # solution to allow broadcasting --> add dimension:
        alpha = np.array(alpha)[None, :] # transformation to array so that vector can get added
        V1, V2, V3 = alpha_to_vec(alpha, lon, lat) # np.cos(alpha)*project_north(lon, lat) + np.sin(alpha)*construct_east(lon, lat)
        # doppelter Winkel: * 2
        Y1, Y2, Y3 = alpha_to_vec( 2*alpha , lon, lat)
        # Skalarprodukt zum Check:
        scalarvec = V1*X1+V2*X2+V3*X3
        # round to three decimal points
        scalarvec = np.round(scalarvec, 3)

        # Column names are (sadly weird):
        # Visibility  0 deg
        # visibili_1  45 deg
        # vis_90deg  90 deg
        # vis_135deg  135 deg
        if 'visibility' in sel_df.columns:
            # SUNAZI Special:
            J00 = sel_df['visibility']
            J45 = sel_df['visibili_1']
            J90 = sel_df['vis_90deg']
            J135 = sel_df['vis_135deg']
        else:
            continue
        # calculate angle between azimuth 
        A00 = np.pi/180*abs(0-sel_df['azimuth'])
        A45 = np.pi/180*abs(45-sel_df['azimuth'])
        A90 = np.pi/180*abs(90-sel_df['azimuth'])
        A135 = np.pi/180*abs(135-sel_df['azimuth'])

        id_int = sel_df['id_int']

        # # check with for loop:
        # alpha = np.pi/180* sel_df['azimuth']
        # for (aX1, aX2, aX3, alph, alon, alat) in (zip(X1, X2, X3, alpha, lon, lat)):
        #     # print(alpha_to_vec(alph, alon, alat))
        #     V1, V2, V3 = alpha_to_vec(alph, alon, alat)
        #     # check scalar product:
        #     scalar_product = aX1*V1 + aX2*V2 + aX3*V3
        #     # print(scalar_product)
        #     if scalar_product > 0.01:
        #         print('fatal')

        # ATTENTION: if I mix pandas dataframe series, which are indexed, and since they are filtered, the indices do not match up with the numpy arrays, it is a mess
        # therefore, re-index or convert to numpy array
        df_traininreg_vMF = pd.DataFrame([np.array(X1), np.array(X2), np.array(X3), Y1, Y2, Y3, V1, V2, V3, alpha[0], np.array(lon), np.array(lat), scalarvec, sunazi, cluster, id_int, J00, J45, J90, J135, A00, A45, A90, A135]).T
        df_traininreg_vMF.columns = ['X1', 'X2', 'X3', 'Y1', 'Y2', 'Y3', 'V1', 'V2', 'V3', 'alpha', 'lon', 'lat', 'scalarproduct_XY', 'sun', 'cluster', 'id_int', 'J00', 'J45', 'J90', 'J135', 'A00', 'A45', 'A90', 'A135']

        # save
        df_traininreg_vMF.to_csv(savepath_csv.joinpath('{}_for_vMF_spherical.csv'.format(df_name)))




    # %%

    # once dff all without filter
    # start fresh
    dff_all = load_df_ALL(input_path)
    dff_children = load_df_ALL_CHILDREN(input_path)

    for (dff, df_name) in [(dff_all, "dff_all"), (dff_children, "dff_children")]:
        # create group index
        dff['cluster'] = dff.groupby(['emission']).ngroup()
        # the cluster
        cluster = dff['cluster']
        sel_df = dff
        ################# Spherical coordinates
        # construct (X1, X2, X3), (Y1, Y2, Y3), (V1, V2, V3)
        # longitude in radians, to east-positive!
        # latitude from -90 to 90 --> change the sign by multiplying by 
        # X1 = cos(lat)*cos(lon)
        # X2 = cos(lat)*sin(lon)
        # X3 = -sin(lat)
        lon = np.pi/180 * (360 - sel_df['lon_d']) # transform longitude to east-positive
        lat = np.pi/180 * (sel_df['lat_d'])
        X1 = np.cos(lat)*np.cos(lon)
        X2 = np.cos(lat)*np.sin(lon)
        X3 = np.sin(lat)
        # 'Time'
        # T = sel_df['FPK']

        # alpha = np.pi/180* (90 - sel_df['azimuth']) # transform again
        # # wait, why 90-azi here? I thought the transformation is based on north-east (x,y) system!!
        alpha = np.pi/180* sel_df['azimuth'] # transform again
        # about the shapes:
        # project_north(lon, lat).shape = (3, N)
        # alpha.shape (N,)
        # solution to allow broadcasting --> add dimension:
        alpha = np.array(alpha)[None, :] # transformation to array so that vector can get added
        V1, V2, V3 = alpha_to_vec(alpha, lon, lat) # np.cos(alpha)*project_north(lon, lat) + np.sin(alpha)*construct_east(lon, lat)
        # doppelter Winkel: * 2
        Y1, Y2, Y3 = alpha_to_vec( 2*alpha , lon, lat)
        # Skalarprodukt zum Check:
        scalarvec = V1*X1+V2*X2+V3*X3
        # round to three decimal points
        scalarvec = np.round(scalarvec, 3)

        # let's add the sun azimuth for each observation. can be used to calculate the sun azimuth for one gridpoint
        sunazi = 90 - sel_df['sunazimuth'] # transformed to mathematical coord system

        if 'visibility' in sel_df.columns:
            # SUNAZI Special:
            J00 = sel_df['visibility']
            J45 = sel_df['visibili_1']
            J90 = sel_df['vis_90deg']
            J135 = sel_df['vis_135deg']
        else:
            continue
        # calculate angle between azimuth 
        A00 = np.pi/180*abs(0-sel_df['azimuth'])
        A45 = np.pi/180*abs(45-sel_df['azimuth'])
        A90 = np.pi/180*abs(90-sel_df['azimuth'])
        A135 = np.pi/180*abs(135-sel_df['azimuth'])

        id_int = sel_df['id_int']

        # ATTENTION: if I mix pandas dataframe series, which are indexed, and since they are filtered, the indices do not match up with the numpy arrays, it is a mess
        # therefore, re-index or convert to numpy array
        df_traininreg_vMF = pd.DataFrame([np.array(X1), np.array(X2), np.array(X3), Y1, Y2, Y3, V1, V2, V3, alpha[0], np.array(lon), np.array(lat), scalarvec, sunazi, cluster, id_int, J00, J45, J90, J135, A00, A45, A90, A135]).T
        df_traininreg_vMF.columns = ['X1', 'X2', 'X3', 'Y1', 'Y2', 'Y3', 'V1', 'V2', 'V3', 'alpha', 'lon', 'lat', 'scalarproduct_XY', 'sun', 'cluster', 'id_int', 'J00', 'J45', 'J90', 'J135', 'A00', 'A45', 'A90', 'A135']


        # save
        df_traininreg_vMF.to_csv(savepath_csv.joinpath('{}_nofilter_for_vMF_spherical.csv'.format(df_name)))



# %%
