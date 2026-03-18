# 2024-03-21
# Caroline Haslebacher
# simple script to delete all files that contain no SSI ID from the list

#%% import
import os
from pathlib import Path
import argparse
from MAIN_loop1_to_loop4 import get_polygeopaths
import pandas as pd
import numpy as np

import shutil



#%%
current = os.getcwd()
titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

#%% define path

if __name__ == '__main__':
    '''
    The basepath indicates where to find the 'geotiffs/original' and 'polygons/original' folder.
    This then makes a new for_analysis subfolder.
    '''
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Get basepath.')

    parser.add_argument('--basepath', required=True,
                    metavar="where to find //geotiffs/for_analysis and //polygons/for_analysis",
                    help='as posix')
    
    parser.add_argument('--ssi_list', required=True,
                    metavar="name of ssi_list csv file to find in basepath. Can be 0 if no list is available.",
                    help='e.g. list_regmaps_ssi_ids.csv')

    args = parser.parse_args()

    geopoly_basepath = Path(args.basepath) # titaniach / 'lineament_detection/galileo_manual_segmentation/'
    # equirectangular
    geo_orig = geopoly_basepath / 'geotiffs/original'
    poly_orig = geopoly_basepath / 'polygons/original'
    geopaths = sorted(geo_orig.glob('*.tif'))
    polypaths = sorted(poly_orig.glob('*.geojson'))

    # read in ssi list
    if args.ssi_list == '0':
        # implement solution of empty list
        # then, we simply take all SSI
        ssi_list = [gp.stem for gp in geopaths]
    else:
        ssi_list = list(pd.read_csv((geopoly_basepath).joinpath(args.ssi_list)))
    # e.g.
    # ['C0466664152R',
    # 'C0466664153R',
    # 'C0466664165R',
    # 'C0466664166R',
    # 'C0466664178R',
    # 'C0466664179R',
    # 'C0466664200R',
    # 'C0466664201R', ...
    os.makedirs(geopoly_basepath / 'polygons' / 'for_analysis', exist_ok=True)
    os.makedirs(geopoly_basepath / 'geotiffs'  / 'for_analysis', exist_ok=True)

    debug_list_poly = []
    debug_list_geo = []
    # loop through and rename with system commands
    # os.rename(source, destination, *, src_dir_fd = None, dst_dir_fd = None)
    for polyp in polypaths:
        splits = polyp.stem.split('_') # e.g ['2024', '02', '29', '02', '21', 'C0466669578R', '0']
        ssi_id = splits[5] # e.g. C0466669578R
        # we also want to rename
        datestr = ''.join(splits[:5]) # e.g. '202402290221'
        new_ssi_id = '_'.join(splits[5:]) # note: this gets a 0 or higher appended due to cut_size algorithm in the stitching tool. e.g. 'C0466669578R_0'
        if ssi_id in ssi_list:
            # print('I copy {}'.format(((polyp.parents[1] / 'for_analysis').joinpath(new_ssi_id + '_' + datestr + '.geojson') )))
            shutil.copyfile(polyp, ((polyp.parents[1] / 'for_analysis').joinpath(new_ssi_id + '_' + datestr + '.geojson')) )
            debug_list_poly.append(ssi_id)

    print('geopaths')
    # print(geopaths)
    for geop in geopaths:
        ssi_id = geop.stem
        if ssi_id in ssi_list:
            # print('I copy {}'.format((geop.parents[1] / 'for_analysis').joinpath(geop.name)))
            debug_list_geo.append(ssi_id)
            shutil.copyfile(geop, (geop.parents[1] / 'for_analysis').joinpath(geop.name))

    print(len(ssi_list))
    print(len(debug_list_poly))
    print(len(debug_list_geo))
    print(np.unique(debug_list_poly) == np.unique(debug_list_geo))


#%% loop through all files


