# 2024-02-22
# CAroline Haslebacher
# This script finally merges the polygons back and calculates properties of the parent lineament,
# which is identifiable with the Unique, fixed ID

# It also produces a second table with all child azimuths, lengths, widths, centre_longitude/latitude left untouched, 
# but with an added crossings per length (=fragmentations per km FPK), which allows us to sort into age bins,
# because one lineament should have the same age (assumption)

#%% import modules

import pandas as pd
from pathlib import Path
import json
import os
import numpy as np
from geojson import dump
from sklearn.utils.extmath import weighted_mode
import argparse
import matplotlib.pyplot as plt
from PIL import Image

from MAIN_loop1_to_loop4 import get_polygeopaths
from Cut_geotiff_to_smaller_geotiffs import geotiff_to_arr

#%%
# current = os.getcwd()
# titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

#%%

def check_df_col(df1, colname):
    '''
    This function simply checks if only one value is in df1[colname].
    For this, we use the .unique() function. If the length is 1, there is only one value.

    Other design possibility: If the mean is exactly equal to each field, it is True
    '''
    colval = df1[colname]
    # if there are many decimal places, sometimes the unique function doesn't work (for example for the resolution_mpx)
    if colval.dtype == 'float64':
        colval = round(colval, 4)
    return len(colval.unique()) == 1

def check_uniqueness(myarr):
    '''
    This function simply checks if the given array has only 1 unique value and returns this unique value
    Else, it raises a Warning
    '''
    if len(myarr.unique()) == 1:
        return myarr.unique()[0]
    else:
        raise Warning('Not unique values in array {}.'.format(myarr))

def weighted_mean(myarr, weights=None):
    '''
    input: numpy array with float or int values. Perhpas weights, same size of myarr or 1D
    output: a weighted mean of the values, rounded to the number of decimal points of the input array
    formula: mean = sum( wi * xi) / sum( wi) 
    or from numpy: avg = sum(a * weights) / sum(weights)
    '''
    # for the number of decimal points, we simply take the first entry, convert it to a string, split it by '.' and count the length of the second entry (e.g. 45.2843 --> ['45', '2843'] --> len('2843)=4)
    num_decimal_points = len(str(myarr[0]).split('.')[-1])
    return round(np.average(myarr, weights=weights), num_decimal_points)


def weighted_mean_error(myarr, weights):
    '''
    input: numpy array with float or int values (entries are called 'std' in description, because for example, myarr are the standard deviations) and corresponding weights, same shape than myarr.
    output: the error (Gaussian error calculation) of the weighted mean, rounded to the number of decimal points of the input array
    formula: error_of_mean = sqrt( sum( (wi*sdti)**2 ) ) / sum( wi) 
    '''
    # for the number of decimal points, we simply take the first entry, convert it to a string, split it by '.' and count the length of the second entry (e.g. 45.2843 --> ['45', '2843'] --> len('2843)=4)
    num_decimal_points = len(str(myarr[0]).split('.')[-1])
    # 
    error_of_weighted_mean = 1/weights.sum() * np.sqrt( np.sum(myarr**2 * weights**2) )
    return round(error_of_weighted_mean, num_decimal_points)

def cast_FPK_to_minus1_based_on_length(FPK, length):
    '''
    This simple function returns -1, if FPK is zero AND the length is smaller than 15 km. We assume that for these lineaments, we cannot say anything about their age.
    For all other cases, it simply returns the input FPK
    '''
    if (FPK == 0) and (length < 15): # 15 km
            return -1
    return FPK

def Merge(dict1, dict2):
    res = {**dict1, **dict2}
    return res

#%%

if __name__ == '__main__': # this allows me to import functions defined above!

    '''
    This script merges the scientifically analysed tessels back together (the feature layers)
    and saves them in the specified location
    This also generates a full merged csv of all individual images!

    test/Debug with:
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    equifile = 'Equirectangular'
    savefolder = 'Equirectangular_EUROPA.prj'
    output_foldername = 'output/testing_scwmanualsegs_algo'
    '''
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Get basepath.')

    parser.add_argument('--basepath', required=True,
                    metavar="where to find /data/geotiffs/for_analysis and /data/polygons/for_analysis",
                    help='as posix')
    parser.add_argument('--crsfoldername', required=False,
                        default='Equirectangular',
                    metavar="the name of sub-folder to save",
                    help='e.g. Equirectangular')  
    parser.add_argument('--crs_file', required=False,
                        default='Equirectangular_EUROPA.prj',
                    metavar="the name of the crs file. Folder is jsonpath.",
                    help='jsonpath.joinpath(equifile) has to be valid.')
    parser.add_argument('--input_foldername', required=False,
                        default='output/extraction_loop3',
                    metavar="the name of sub-folder to read",
                    help='e.g. output/extraction_loop3') 
    parser.add_argument('--output_foldername', required=False,
                        default='output/extraction_loop4_csv',
                    metavar="the name of sub-folder to read",
                    help='e.g. output/extraction_loop4_csv')  
    # instead of bool: action='store_true': This sets the argument to True if the flag is present on the command line and False if it is omitted. This is suitable for flags that enable a feature, where the default state is False.
    parser.add_argument('--illumination_info', action='store_true', 
                    help='use this flag if you want to attach illumination information from a Galileo SSI image of Europa.') 


    args = parser.parse_args()

    basepath = Path(args.basepath) # titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    crsfolder = args.crsfoldername
    crs_file = args.crs_file
    input_fol = args.input_foldername
    output_fol = args.output_foldername
    illum_info = args.illumination_info
    # e.g.
    # basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    '''
    to debug, use for example:
    current = os.getcwd()
    titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data_phoc_sunAzi'
    crsfolder = 'Equirectangular' 
    crs_file = 'Equirectangular_EUROPA.prj'
    input_fol = 'output_phocSunAzi/extraction_loop3'
    output_fol = 'output_phocSunAzi/extraction_loop4'
    illum_info = True    
    
    basepath = titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/do_analysis'
    crsfolder = 'UFID'
    
    basepath = titaniach / 'lineament_detection/RegionalMaps_CH_NT_EJL_LP/mapping/data'
    crsfolder = 'Equirectangular' 
    input_fol = 'output/extraction_loop3'
    output_fol = 'output/extraction_loop4'

    # LineaMapper real debug
    basepath = titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/analyse_LM1_1'
    crsfolder = 'Equirectangular/valid'
    input_fol = 'output/LM1_1_extraction_loop3'
    output_fol = 'output/LM1_1_extraction_loop4'
    SSI_id = 'true'
    polystem_0 = 'C0449974429R_0_202407060009'
    polystem_1 = 'C0449974429R_1_202407060119'
    polystem_2 = 'C0449974429R_2_202407060158'
    polystem_4 = 'C0449974429R_4_202407060234'
    polystem_5 = 'C0449974429R_5_202407060352'
    polystems = [polystem_0, polystem_1, polystem_2, polystem_4, polystem_5]
    # to test if the result does not change for others
    polystem = 'C0449961813R_0_202407051201'
    '''
    #%%
    # equirectangular
    geo_orig = basepath / 'geotiffs/for_analysis'
    poly_orig = basepath / 'polygons/for_analysis'
    #####

    # define here where to save csv files in the end
    savepath_csv = basepath.parents[0] / output_fol # e.g. output/extraction_loop4_csv
    os.makedirs(savepath_csv, exist_ok=True) 
    os.makedirs(savepath_csv / 'individual_polystems', exist_ok=True)
    os.makedirs(savepath_csv / 'individual_polystems_children', exist_ok=True)
    os.makedirs(savepath_csv / 'lineament_density', exist_ok=True)


    import_basepath = basepath.parents[0] / input_fol

    polystems = [f.name for f in os.scandir(import_basepath) if f.is_dir()] # ['11ESCOLORS01-01', '11ESREGMAP01EXCERPT1', ...

    # read projection file:
    crspath = poly_orig / crsfolder
    with open(crspath.joinpath(crs_file)) as prjfile:
        proj = prjfile.read()
        print(proj)

    #%% get original (Equirectangular) polygon paths
    # take the 'valid' folder, if it exists

    geopaths, polypaths = get_polygeopaths(poly_orig / crsfolder, geo_orig)

    # this is for saving the original geojson files with the appended calculated fields
    savegeojson = basepath / 'polygons/analysed' / crsfolder
    os.makedirs(savegeojson, exist_ok=True)
    # we also save the children merged geojson
    savegeojson_CHILDREN = basepath / f'polygons/analysed/{crsfolder}_CHILDREN'
    os.makedirs(savegeojson_CHILDREN, exist_ok=True)

    #%% 

    ''''
    Table 1: Merged parent lineaments only

    # FOR PARENT LINEAMENT:
    # weighted average width. weighted by child length.
    # sum length
    # sum area
    # weighted mode (!) of azimuth. weighted by child length.
    # error propagation of RMS/width_px from azimuth weighted mean
    # mean of longitude
    # mean of latitude
    # weighted average of meanIF. weighted by child area.
    # error propagation of stdIF from meanIF weighted mean

    # finally, FPK is calculated, but set to -1 if the FPK is zero AND the length is smaller than 15 km

    # calculate Fragmentations per Kilometer (FPK). This is important for Table 2.

    Description of loop system logic:
    A) we loop through all folders in //titania.unibe.ch/Space/Groups/PIG/Caroline/lineament_detection/galileo_manual_segmentation/output/extraction_loop3/
    B) We loop through all geojson files in the current directory, for example //titania.unibe.ch/Space/Groups/PIG/Caroline/lineament_detection/galileo_manual_segmentation/output/extraction_loop3/11ESREGMAP01EXCERPT1/
    , read each geojson file and append it to one dataframe per subdirectory. In there are all children which were analyzed in science_with_manual_segs_algo.py.
    C) Now, we merge the children. For this, we group the dataframe by Unique Fixed ID (UFID), which is the same for each child, 
    and calculate the properties of the parent lineament out of the children properties.
    We append to the list df_parents_list, which contains all parent lineaments in all subdirectories together.
    We also append each child with an added UFID, num_subpolygons, and FPK, to a special df_childs_list, because taking the weighted average of the azimuth is not appropriate for an azimuth analysis.

    We add the parent properties back to the original untessellated equirectangular geojson, which facilitates validation in ArcGIS. But with the UFID, we can also cross-check directly with the csv table

    We also save the merged children geojson. This helps to see which parts belong to one lineament.

    In addition, we sum the lineament total and fractional areas from the csv files

    '''

    # initialize lists
    df_parents_list = []
    df_childs_list = []

    ########## A)
    for polystem in polystems:
        print(polystem)
        dfs_list = []
        df_tosave_individuals = [] # parents
        df_tosave_individual_children = []
        geoms_children = {}
        # somehow, the 'crs' wasn't copied correctly, therefore, we enter this manually.
        # FOr the children, this is all correct, because they are read in from 'back_to_prj' folder!
        geoms_children['crs'] = {'type': 'name', 'properties': {'name': proj}} 
        # was once hardcoded to: 'PROJCS["Equirectangular_EUROPA",GEOGCS["GCS_EUROPA",DATUM["D_EUROPA",SPHEROID["EUROPA_localRadius",1560800.0,0.0]],PRIMEM["Reference_Meridian",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Equidistant_Cylindrical"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",180.0],PARAMETER["Standard_Parallel_1",0.0],UNIT["Meter",1.0]]'
        geoms_children['features'] = [] # children geojson dict to append/update
        # for each folder, we first get a list of all saved geojson files. NOTE: here, we take the back-projected geojson files. (problem: all sub-features have a different orthographic coordinate system!) The attributes are taken from orthographic, while the coordinates are equirectangular for merging
        calcs_paths = sorted((import_basepath / polystem / 'back_to_prj/json_files').glob('*.geojson'))
        # we load each geojson file, convert it to a table 
        ########## B)
        for calc_path in calcs_paths:
            with open(calc_path) as jsf:
                geoms = json.load(jsf)
                # geoms_children['features'].append(geoms['features'].copy())
                geoms_children['features'] = geoms['features'].copy() + geoms_children['features'].copy()
                # print('length: {}'.format(len(geoms['features'])))
            # make 'UFID' and 'num_subpolygons' its own column, so that we can use it to group by.
            geoms_df = pd.DataFrame(geoms['features'])
            geoms_df['UFID'] = [x['properties']['UFID'] for _, x in geoms_df.iterrows()]
            geoms_df['num_subpolygons'] = [x['properties']['num_subpolygons'] for _, x in geoms_df.iterrows()]
            # tested with: geoms_df['properties'][0]['UFID'] as an example for the first row
            # now, delete 'UFID' from properties
            for _, x in geoms_df.iterrows():
                del x['properties']['UFID']
                del x['properties']['num_subpolygons']
            # convert to DataFrame (UFID is the identifier )
            dfs_list.append(geoms_df)
        # children geoms: append also other information. Keys of geoms are dict_keys(['type', 'name', 'features', 'crs'])
        # I reproject to Equirectangular in science_with_manual_segs_algo.py!
        geoms_children['type'] = geoms['type']

        # read in csv files with areas
        # 1 polystem must lead to one summed area csv
        # the good news is that dataframes with equal structure of columns and rows can easily be summed with df1+df2
        # however, for the fractions, we need the average
        # therefore, we simply divide the sum 
        csv_paths = sorted((import_basepath / polystem / 'csv_area_files').glob('*_linea_area.csv'))
        # loop through all csv files (of the children)
        area_semantic_dict = {'1': {'fraction': 0}, 
                        '2': {'fraction': 0}, 
                        '3': {'fraction': 0}, 
                        '4': {'fraction': 0}, 
                        'total_linea': {'fraction': 0},
                        'total_area': {'fraction': 0},
                       }
        
        # retrieve the number of the subimage (from LineaMapper_to_img.py)
        # for example 5 from 'C0449974429R_5_202407060352'
        subimg_num = polystem.split('_')[1]
        # check now if there are other subimages, or if this image was not split in LineaMapper_to_img.py
        # we introduce a helper boolean to indicate whether or not we should substitute the total area (by taking it from the saved images in geo_orig / 'clipped', with Cut_geotiffs....py)
        subst_area = False
        # we test if we find multiple polystems with the same polystem SSI id
        ssi_id = polystem.split('_')[0] # e.g. C0449974429R from 'C0449974429R_5_202407060352'
        if len([x for x in polystems if ssi_id in x]) > 1:
            # if the length of the list with all polystems that contain the ssi_id of the current polystem is greater than 1, there are subimages.
            # therefore, we set the boolean to True and substitute the area.
            subst_area = True

        if subst_area == True:
            # if this number is greater than 0, we need to adjust the total_area
            # OR: we just do this for all...!? for most of them, it should not make a difference
            # read in the dataset and extract the area:
            subimg_path = (geo_orig / 'clipped').joinpath(polystem + '.tif')
            subimg, dataset = geotiff_to_arr(subimg_path)
            ulx, xres, xskew, uly, yskew, yres  = dataset.GetGeoTransform()
            total_area_subimg_px = len(subimg[subimg>0])
            total_area_subimg_km2 = abs(total_area_subimg_px*xres*yres*1e-6)


        for cs_idx, csvf in enumerate(csv_paths):
            # read csv
            # we need to specify index_col=0 to make sure 'fraction' and 'km2' are read as indexes, otherwise they get added up, this leads to 'fractionfractionfractionfractionfractionfracti...' and 'km2km2km2km2km2km...'
            dfarea = pd.read_csv(csvf, index_col=0)
            # print(dfarea)
            # check up:
            # checktot = round(dfarea['total_linea']['km2'] / dfarea['total_area']['km2'], 5) == round(dfarea['total_linea']['fraction'], 5)
            # if checktot == False:
            #     print(csvf)
            #     print(dfarea)
            #     print(round(dfarea['total_linea']['km2'] / dfarea['total_area']['km2'], 5))
            #     print(round(dfarea['total_linea']['fraction'], 5))
            #     raise ValueError
            # first dataframe is intialisation, after, we sum
            if cs_idx == 0:
                df_area_polystem = dfarea.copy()
                continue # else, we count the first twice
            # sum
            df_area_polystem += dfarea
        
        # calculate fractional average
        # Can it happen that one dataframe is entirely zero? yes, and it should be counted. Empty tiles are sorted very early in the process, no problem with that.
        # df_area_polystem = df_area_polystem.apply(lambda x: x/len(csv_paths) if x.name == 'fraction' else x, axis = 1)

        # update total area with clipped geotiff
        if subst_area == True:
            print('old total area was: {}'.format(df_area_polystem.loc['km2', 'total_area']))
            df_area_polystem.loc['km2', 'total_area'] = total_area_subimg_km2
            print('new area is: {}'.format(df_area_polystem.loc['km2', 'total_area']))
        # else, it is already correct

        # calculate fractions from totals. Mind that else, we would need to apply a weighted average for summing the fractions!!!
        # therefore, drop the original fraction:
        for fracidx in range(1,5):
            area_semantic_dict[str(fracidx)]['fraction'] = df_area_polystem[str(fracidx)]['km2']/ df_area_polystem['total_area']['km2']
        # total fraction
        area_semantic_dict['total_linea']['fraction'] =  df_area_polystem['total_linea']['km2'] / df_area_polystem['total_area']['km2']
        area_semantic_dict['total_area']['fraction'] = df_area_polystem['total_area']['fraction'] / len(csv_paths) # should be 1
        # now, we can drop the fraction
        df_area_polystem.drop('fraction', inplace=True)
        # area semantic dict to dataframe to replace 'fraction' row
        # to do this, we need to have the exact same column names (might not be the case for string and int, for example)
        df_fraction = pd.DataFrame(area_semantic_dict)
        # concat and overwrite:
        df_area_polystem =pd.concat([df_fraction, df_area_polystem])
        # now, we can save the csv
        df_area_polystem.to_csv((savepath_csv / 'lineament_density').joinpath('{}_parents_area.csv'.format(polystem)))
        # and directly make a pie-chart
        # below IS WRONG DUE TO NON-FULLY SEMANTIC CATEGORICAL FRACTIONs
        # alllineamentfracs = df_area_polystem[['1', '2', '3', '4']].loc['fraction']*area_semantic_dict['total_linea']['fraction'] # multiply by total linea fraction. 
        # # add 'non-lineament' fraction
        # non_lineamentfrac = pd.Series(data=1-area_semantic_dict['total_linea']['fraction'], index=['5'])
        # sizes = pd.concat([alllineamentfracs, non_lineamentfrac])
        # check that 1-area_semantic_dict['total_linea']['fraction'] is positive. else, take 0 instead.
        # also, we set 1 if area_semantic_dict['total_linea']['fraction'] is larger than 1 (this can happen if the polygon exceeds the image borders)
        sizes = [min(1, area_semantic_dict['total_linea']['fraction']), max(1-area_semantic_dict['total_linea']['fraction'], 0)]
        # print(polystem)
        # print(sizes)
        # print(df_area_polystem)
        fig_pie, ax = plt.subplots(figsize=(7, 7))
        # tx = ax.pie(sizes, labels=['band', 'double ridge', 'ridge complex', 'undiff. linea', 'non-lineament'],
        #     colors=['#7F4A9D', '#ED9A22', '#8ED311', '#00FFC5', 'maroon'], autopct='%1.1f%%', radius=1.5,
        #     pctdistance=0.8, labeldistance=1.1)
        # NOTE: because the total linea area fraction is NOT the sum of the categorical fractions due to non-semantic overlaps, we only plot fully-semantic lineament vs non-lineament fraction
        tx = ax.pie(sizes, labels=['lineament', 'other'],
            colors=['maroon', 'grey'], autopct='%1.1f%%', radius=1.5,
            pctdistance=0.5, labeldistance=1.1, textprops={'fontsize': 20})
        # save figure
        fig_pie.savefig((savepath_csv / 'lineament_density').joinpath('{}_parents_lindens_pie.png'.format(polystem)), dpi=300, bbox_inches='tight')
        # make barchart with categorical fractions:
        fig_bar, ax = plt.subplots(figsize=(7, 11))
        colors=['#7F4A9D', '#ED9A22', '#8ED311', '#00FFC5']
        labels=['band', 'double\nridge', 'ridge\ncomplex', 'undiff.\nlinea']
        alllineamentfracs = df_area_polystem[['1', '2', '3', '4']].loc['fraction']
        for color, fracidx, lab in zip(colors, np.arange(0,4), labels):
            barval = ax.bar(1+fracidx, alllineamentfracs.iloc[fracidx], width=0.5, color=color, label=lab)
            ax.bar_label(barval, fontsize=15, fmt='%.2f', padding=0.5, color='black') 
        ax.set_ylim(0, 1)
        ax.set_xticks([x+1 for x in np.arange(0,4)])
        ax.set_xticklabels(labels, fontsize=15)
        # save
        fig_bar.savefig((savepath_csv / 'lineament_density').joinpath('{}_parents_lindens_bar.png'.format(polystem)), dpi=300, bbox_inches='tight')      
        # # combined: --> decided not to do this, since quality of drawn canvas is bad for the fonts
        # fig_bar.canvas.draw()
        # pil_bar = Image.frombytes('RGB', fig_bar.canvas.get_width_height(),fig_bar.canvas.tostring_rgb())
        plt.close('all')
        ##################

        # join the dataframes
        df = pd.concat(dfs_list)
        # TODO: load original geojson in equirectangular projection (for visualization)
        # After parent calculation, append analysis information by matching the UFID
        # match polypath:
        load_polypath = 'null'
        for json_orig_path in polypaths:
            if polystem == json_orig_path.stem: # ATTENTION, THEY NEED TO BE EQUAL!
                load_polypath = json_orig_path 
        if load_polypath == 'null':
            # then, we did not find any json_path!
            print(polystem)
            raise Warning('no json_orig_path found.')
        # now we can load the matched 'load_polypath'
        with open(load_polypath) as jsf_orig:
            geoms_orig = json.load(jsf_orig)

        # TODO: this is the mistake for orthographic!!!! go find original crs!
        geoms_orig['crs'] = {'type': 'name', 'properties': {'name': 'PROJCS["Equirectangular_EUROPA",GEOGCS["GCS_EUROPA",DATUM["D_EUROPA",SPHEROID["EUROPA_localRadius",1560800.0,0.0]],PRIMEM["Reference_Meridian",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Equidistant_Cylindrical"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",180.0],PARAMETER["Standard_Parallel_1",0.0],UNIT["Meter",1.0]]'}} 
        # we append the original name to the merged children geojson
        # but there does not have to be a name
        try:
            geoms_children['name'] = geoms_orig['name']
        except KeyError:
            geoms_children['name'] = load_polypath.stem # construct
        ########## C)
        # now, we groupy by the Unique Fixed id and loop through the 'sub-dataframes'
        for gridx, dfgr in enumerate(df.groupby('UFID')):
            df_props_list = []
            for index, row in dfgr[1].iterrows():
                # print(row['properties'], row['UFID'])
                df_props_list.append(pd.DataFrame(row['properties'], index=[index]))
            # concat df_props
            df_props = pd.concat(df_props_list)
            # instantiate a dictionary for the parent. This makes it easy to transfer to a dictionary later
            parent_dict = {}
            # first, add constant fields to the parent dict: (also serves as validity check, therefore, we do this first!)
            if illum_info:
                constantcols = ['id_int', 'comments', 'unit', 'incidence', 'phase', 'emission', 'res_mpx']
                # 2025-12-02: 'sunazimuth' was removed from this list due to introduction of phocube extraction. added further below
            else:
                constantcols = ['id_int', 'comments', 'unit', 'res_mpx']
            for colname in constantcols:
                if check_df_col(df_props, colname) == False: # function output should normally be True (except for resolution)
                    # if colname == 'res_mpx':
                    #    print('resolution not equal') # THE RESOLUTION IS NOT EQUAL DUE TO GDAL FUNCTIONS, I ASSUME.
                    if colname != 'res_mpx':
                        print(df_props[colname])
                        # I had one case where the unit for one child was assigned 'NaN'. is this because of the predicted empty category?
                        # therefore, try to drop NAN
                        df_props.dropna(axis=0, inplace=True) # --> this does nothing really. I would have to drop it in parent_dict! Can't think right now how to safely do this. So I just drop NaN later
                        # check again
                        if check_df_col(df_props, colname) == False: 
                            print(df_props[colname])
                            raise Warning('A value that should be constant was not constant throughout the childs!')
                        else:
                            print('all good now')
                # append to parent lineament. Because they are all equal, we can take now the unique value
                if colname == 'id_int':
                    parent_dict[colname] = float(df_props[colname].unique()[0]) # otherwise geojson dump is unhappy
                else:
                    parent_dict[colname] = df_props[colname].unique()[0]
            # calculate properties of parent lineament:
            # first, add UFID and num_subpolygons
            UFID = dfgr[0]
            length = df_props['length_km'].sum()
            if length == 0.0:
                print('WARNING: FOUND EMPTY LENGTH LINEAMENT.')
                # TODO:  I will not save this one, as most likely this is an empty case. how to do this safely? see concern above.
            fragmentations = check_uniqueness(df['num_subpolygons'][df['UFID'] == dfgr[0]]).copy()
            FPK = round(cast_FPK_to_minus1_based_on_length(fragmentations/length, length), 3)
            parent_dict['UFID'] = UFID # this is the first entry of the tuple (since we grouped by UFID, no need to check uniqueness)
            parent_dict['fragmentations'] = float(fragmentations) # we mask the dataframe by asking for the subdataframe of only the Unique Fixed ID and then taking the number of subpolygons from this
            # for each parent, calculate FPK
            parent_dict['FPK'] = FPK
            # then, add calculated fields
            parent_dict['length_km'] = round(length, 2)
            parent_dict['area_km2'] = round(df_props['area_km2'].sum(), 2)
            parent_dict['azimuth'] = round(weighted_mode(np.array(df_props['azimuth']), np.array(df_props['length_km']))[0][0], 2) # output example: (array([150.08]), array([33.9566])), we take the first entry of the first array
            parent_dict['azimuth_std'] = round(np.nanstd(np.array(df_props['azimuth'])), 2) # no weights here
            parent_dict['RMS_pwpx'] = weighted_mean_error(np.array(df_props['RMS_pwpx']), weights=np.array(df_props['length_km']))
            # print('mean: {}, mode: {}'.format(weighted_mean(np.array(df_props['azimuth']), np.array(df_props['length_km'])), weighted_mode(np.array(df_props['azimuth']), np.array(df_props['length_km']))))
            parent_dict['RMS'] = weighted_mean_error(np.array(df_props['RMS']), weights=np.array(df_props['length_km']))
            parent_dict['width_km'] = weighted_mean(np.array(df_props['width_km']), weights=np.array(df_props['length_km']))
            parent_dict['lon_d'] = weighted_mean(np.array(df_props['lon_d'])) # we use the function without weights
            parent_dict['lat_d'] = weighted_mean(np.array(df_props['lat_d']))
            parent_dict['meanIF'] = weighted_mean(np.array(df_props['meanIF']), weights=np.array(df_props['area_km2']))
            parent_dict['stdIF'] = weighted_mean_error(np.array(df_props['stdIF']), weights=np.array(df_props['area_km2']))
            if illum_info:
                # length-weighted sunazimuth
                # # debugging:
                # if weighted_mean(np.array(df_props['sunazimuth']), weights=np.array(df_props['length_km'])) > 290:
                #     print(df_props['sunazimuth'])
                #     print(df_props['length_km'])
                parent_dict['sunazimuth'] = weighted_mean(np.array(df_props['sunazimuth']), weights=np.array(df_props['length_km']))

            # now we construct the dataframe from the dict
            df_parent = pd.DataFrame(parent_dict, index=[gridx])
            # print(df_parent)
            # drop rows that are zero throughout:
            df_parent = df_parent.loc[~(df_parent == 0).all(axis=1)]
            df_props = df_props.loc[~(df_props == 0).all(axis=1)]

            # append to final parents list
            df_parents_list.append(df_parent)
            # and to the list which we save individually (yes, we could make this shorter, but it's fine)
            df_tosave_individuals.append(df_parent)

            # add UFID, fragmentations, and FPK to childs dataframes
            df_props.insert(0, "UFID", [UFID for x in range(len(df_props))], True)
            df_props.insert(1, "fragmentations", [fragmentations for x in range(len(df_props))], True)
            df_props.insert(2, "FPK", [FPK for x in range(len(df_props))], True)
            # append to childs list
            df_childs_list.append(df_props)
            # and to save individually
            df_tosave_individual_children.append(df_props)

            # append to original geojson
            # loop through all features and find the one with the correct UFID
            for origfeature in geoms_orig['features']:
                if origfeature['properties']['UFID'] == dfgr[0]:
                    # then replace the properties with A COPY of the full parent_dict
                    # print(origfeature)
                    origfeature['properties'] = parent_dict.copy()


        # save individually
        # PARENTS
        df_polystem = pd.concat(df_tosave_individuals)
        df_polystem.to_csv((savepath_csv / 'individual_polystems').joinpath('{}_parents_lineaments.csv'.format(polystem)))
        # CHILDREN
        df_polystem_children = pd.concat(df_tosave_individual_children)
        df_polystem_children.to_csv((savepath_csv / 'individual_polystems_children').joinpath('{}_children_lineaments.csv'.format(polystem)))
        # Save appended GEOJSON
        # save and close polygon again
        # TODO and implement in the future: drop 'empty' (zero throughout) features
        with open(savegeojson.joinpath(load_polypath.name), 'w') as f:
            dump(geoms_orig, f)
        with open(savegeojson_CHILDREN.joinpath(load_polypath.name), 'w') as f:
            dump(geoms_children, f)

    # concat
    df_all_parents = pd.concat(df_parents_list)
    df_all_children = pd.concat(df_childs_list)

    # save
    df_all_parents.to_csv(savepath_csv.joinpath('ALL_parents_lineaments.csv'))
    df_all_children.to_csv(savepath_csv.joinpath('ALL_children_lineaments.csv'))


    #%% convert dumped parent and children geojson to shapefile

    os.makedirs(savegeojson / 'shape_files', exist_ok=True)
    os.makedirs(savegeojson_CHILDREN / 'shape_files', exist_ok=True)

    geojfiles = sorted((savegeojson).glob('*.geojson'))
    for geojfile in geojfiles:
        # convert PARENTS to shapefile as well
        command = 'ogr2ogr -nlt POLYGON -overwrite -skipfailures {} {}'.format((savegeojson / 'shape_files').joinpath(geojfile.stem + '.shp'), savegeojson.joinpath(geojfile.stem + '.geojson'))
        print(command)
        os.system(command)
        # convert CHILDREN to shapefile as well
        # note that I've adapted all paths, but the names should be the same!
        command = 'ogr2ogr -nlt POLYGON -overwrite -skipfailures {} {}'.format((savegeojson_CHILDREN / 'shape_files').joinpath(geojfile.stem + '.shp'), savegeojson_CHILDREN.joinpath(geojfile.stem + '.geojson'))
        print(command)
        os.system(command)



    #%%
    # make sub-csvs for categories!

    bands_df = df_all_parents[df_all_parents['id_int']==1]
    doublridge_df = df_all_parents[df_all_parents['id_int']==2]
    ridgecomplex_df = df_all_parents[df_all_parents['id_int']==3]
    undifflinea_df = df_all_parents[df_all_parents['id_int']==4]
    # save!
    bands_df.to_csv(savepath_csv.joinpath('BANDS_parents_lineaments.csv'))
    doublridge_df.to_csv(savepath_csv.joinpath('DOUBLERIDGES_parents_lineaments.csv'))
    ridgecomplex_df.to_csv(savepath_csv.joinpath('RIDGECOMPLEXES_parents_lineaments.csv'))
    undifflinea_df.to_csv(savepath_csv.joinpath('UNDIFFERENTIATEDLINEAE_parents_lineaments.csv'))
    # RC and bands together:
    rc_bands_df = pd.concat([bands_df, ridgecomplex_df])
    rc_bands_df.to_csv(savepath_csv.joinpath('BANDS_and_RIDGECOMPLEXES_parents_lineaments.csv'))
    # UL and DB together
    ul_db_df = pd.concat([doublridge_df, undifflinea_df])
    ul_db_df.to_csv(savepath_csv.joinpath('DOUBLERIDGES_and_UNDIFFERENTIATED_parents_lineaments.csv'))

    # and the same for children:
    bands_df = df_all_children[df_all_children['id_int']==1]
    doublridge_df = df_all_children[df_all_children['id_int']==2]
    ridgecomplex_df = df_all_children[df_all_children['id_int']==3]
    undifflinea_df = df_all_children[df_all_children['id_int']==4]
    # save!
    bands_df.to_csv(savepath_csv.joinpath('BANDS_children_lineaments.csv'))
    doublridge_df.to_csv(savepath_csv.joinpath('DOUBLERIDGES_children_lineaments.csv'))
    ridgecomplex_df.to_csv(savepath_csv.joinpath('RIDGECOMPLEXES_children_lineaments.csv'))
    undifflinea_df.to_csv(savepath_csv.joinpath('UNDIFFERENTIATEDLINEAE_children_lineaments.csv'))

    # now also filter curved features.
    # rc bands
    rc_bands_df_curved = rc_bands_df[rc_bands_df['azimuth_std'] > 15]
    rc_bands_df_straight = rc_bands_df[rc_bands_df['azimuth_std'] <= 15]
    # save
    rc_bands_df_curved.to_csv(savepath_csv.joinpath('BANDS_and_RIDGECOMPLEXES_parents_lineaments_curved.csv'))
    rc_bands_df_straight.to_csv(savepath_csv.joinpath('BANDS_and_RIDGECOMPLEXES_parents_lineaments_straight.csv'))
    # ul db
    ul_db_df_curved = ul_db_df[ul_db_df['azimuth_std'] > 15]
    ul_db_df_straight = ul_db_df[ul_db_df['azimuth_std'] <= 15]
    # save
    ul_db_df_curved.to_csv(savepath_csv.joinpath('DOUBLERIDGES_and_UNDIFFERENTIATED_parents_lineaments_curved.csv'))
    ul_db_df_straight.to_csv(savepath_csv.joinpath('DOUBLERIDGES_and_UNDIFFERENTIATED_parents_lineaments_straight.csv'))

    #%% 

    # # plotting for a first look
    # import matplotlib.pyplot as plt
    # plt.scatter(df_all_parents['length_km'], df_all_parents['fragmentations'])


    #%% Age-binned tables: For the temporal azimuth analysis, I produce tables that contain all child azimuths sorted into a specified number of age bins (rather FPK bins)



