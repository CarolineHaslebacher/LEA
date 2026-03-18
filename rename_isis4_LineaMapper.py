# 2024-03-20
# Caroline Haslebacher
# This script quickly renames LineaMapper output from 2024_02_28_17_06_C0349875126R_0.geojson to C0349875126R0_2024_02_28_17_06.geojson

#%%

import os
from pathlib import Path


#%%
current = os.getcwd()
titaniach = Path(current.split('Caroline')[0]) / 'Caroline'

#%% define path

polypath =  titaniach / 'lineament_detection/RegionalMaps_CH_NT_EJL_LP/extraction_from_LineaMapper/data/polygons/for_analysis'

# get all geojson files in the directory
polys = sorted(polypath.glob('*.geojson'))

# loop through and rename with system commands
# os.rename(source, destination, *, src_dir_fd = None, dst_dir_fd = None)
for polyp in polys:
    splits = polyp.stem.split('_') # e.g ['2024', '02', '29', '02', '21', 'C0466669578R', '0']
    datestr = ''.join(splits[:5]) # e.g. '202402290221'
    new_ssi_id = '_'.join(splits[5:]) # note: this gets a 0 or higher appended due to cut_size algorithm in the stitching tool. e.g. 'C0466669578R_0'
    new_name = polyp.parents[0].joinpath(new_ssi_id + '_' + datestr + '.geojson')
    os.rename(polyp, new_name)

# %%
