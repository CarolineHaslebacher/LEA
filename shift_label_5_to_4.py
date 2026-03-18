# Caroline Haslebacher
# 2024-05-10
# Finally, I wrote a script that I can run only once, to shift labels 5 to 4 for older .geojson files

# input: geojson in folder Z:\Groups\PIG\Caroline\lineament_detection\galileo_manual_segmentation\data\polygons\for_analysis\shift5to4
# output: geojson in folder Z:\Groups\PIG\Caroline\lineament_detection\galileo_manual_segmentation\data\polygons\for_analysis\

#%% import

import os
from pathlib import Path
import json
from geojson import dump
# from osgeo import gdal
import numpy as np
from datetime import datetime
import argparse

#%% define path to titania

current = os.getcwd()
titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

#%%
geopoly_basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/'
# equirectangular
poly_orig = geopoly_basepath / 'data/polygons/for_analysis/shift5to4'
polypaths = sorted(poly_orig.glob('*.geojson'))

savegeojson = geopoly_basepath / 'data/polygons/for_analysis'

###### assign fixed id (date and miliseconds?)

for polypath in polypaths:
    print(polypath)
    del_idcs = []
    # open geojson
    with open(polypath) as jsf:
        geoms = json.load(jsf)
    # loop through features
    for i, feature in enumerate(geoms['features']):
        # check if id_int is equal to 4. if yes, continue with next mask_idx (then, it is a cusp and we delete it)
        if geoms['features'][i]['properties']['id_int'] == 4:
            # delete entry
            del_idcs.append(i)
            # print('deleting this one later')
            continue
        elif geoms['features'][i]['properties']['id_int'] == 5:
            # print('shifting')
            geoms['features'][i]['properties']['id_int'] = 4


    # delete cycloid cusps:
    del_idcs.reverse()
    for delidx in del_idcs:
        del geoms['features'][delidx]

    # save and close polygon again
    with open(savegeojson.joinpath(polypath.name), 'w') as f:
        dump(geoms, f)


# %%
