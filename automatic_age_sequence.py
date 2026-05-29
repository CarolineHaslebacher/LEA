# 2026-05-29
# Caroline Haslebacher
# Automatic age sequence generation from lineament polygons

#%% import modules
import os
import argparse
from pathlib import Path
import numpy as np
import cv2
from datetime import datetime
import time

from science_with_manual_segs_algo import * # cut_bounding_box, load_tiff_arr_specified, geojson_to_mask



#%% define functions
def mask_to_MINAREARECT(masktf):
    '''
    input is a mask to fit (masktf), a numpy (!) sparse COO matrix (https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_array.html#scipy.sparse.coo_array )
    output: the minimum area rectangle's coordinates

    test with, e.g.:
    mask_idx = 3
    masktf = sparse_masks[mask_idx]
    beta = fit_line_to_mask(sparse_masks[mask_idx])
    print(beta)
    plt.imshow(masks[:,:,mask_idx])
    '''
    # first, generate a 'contours'. list of numpy arrays
    ftl_idcs = [ np.array([[x, y]]) for x, y in  zip(masktf.row.tolist(), masktf.col.tolist() ) ]
    # stack to array of shape [Npoints, 1, 2] (2 for the x/y pair)
    ftl_idcs = np.stack(ftl_idcs)
    centerpnts, (p01, p03), alpha = cv2.minAreaRect(ftl_idcs)

    return (p01, p03)


#%% main
if __name__ == '__main__': # this allows me to import functions defined above!

    # '''
    # This script reads in the tessellated orthographically projected geotiff and geojson and extracts useful lineament characteristics from it.
    # This is then saved in basepath.parents[0] / 'output' as children (reprojected to equirectangular) and stitched back to parents later.
    # '''
    # # Parse command line arguments
    # parser = argparse.ArgumentParser(
    #     description='Get basepath.')

    # parser.add_argument('--basepath', required=True,
    #                 metavar="where to find /data/geotiffs/for_analysis and /data/polygons/for_analysis",
    #                 help='as posix')

    # args = parser.parse_args()

    # basepath = Path(args.basepath) # titaniach / 'lineament_detection/galileo_manual_segmentation/data'



    area_limit = 20 # masks filled with fewer than XX pixels get ignored
    '''
    test/Debug with:

    '''
    
    #### for  development 
    current = os.getcwd()
    titaniach = Path(current.split('Caroline')[0]) / 'Caroline'
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    outputfoldername = 'output/automatic_age_sequence'
    ######


    geo_orig = basepath / 'geotiffs/for_analysis'
    poly_orig = basepath / 'polygons/for_analysis' # this path gets redefined later!
    geopaths, polypaths = get_polygeopaths(poly_orig, geo_orig) # might be adapted

    # for dev,overwrite polypaths:
    polypaths = [basepath / 'polygons/for_analysis' / 'E4ESMACSTR01_GalileoSSI_Equi-cog.geojson'] # Jesse has manually made an age sequence for this polygon


    # specify the basepath (used for saving)
    savebasepath = basepath.parents[0] / outputfoldername # e.g. titaniach / 'lineament_detection/galileo_manual_segmentation/output/extraction_loop3'
    os.makedirs(savebasepath, exist_ok=True)

    # measure time
    start_all = time.time()

    now = datetime.now()
    dt_string = now.strftime("%Y_%m_%d_%H_%M_")

#%%
    # loop through polygons and geotiffs
    # we take the full polygons into account in this script, as opposed to fragmented snippets.
    # orthographic re-projection needed? --> can happen for each polygon
    for pidx, polypath in enumerate(polypaths): # 17ESREGMAP03: [8:9] if all polypaths
        # to test:
        # polypath = polypaths[214]
        # pidx = 214 # for debugging sunazimuth of 304
        print(polypath)
        # get correct geopath. It should by design be the same index, but we might not start at 0, therefore, search
        # make sure we get really the right geopath index (.tif file from .geojson)
        geopath = 0 # security measure
        for gidx, tempgeopath in enumerate(geopaths):
            if polypath.stem == tempgeopath.stem:
                geopath = tempgeopath
        if geopath == 0:
            print('No geopath found for {}. I CONTINUE WITH NEXT POLYPATH.'.format(polypath))
            continue

        polypaths_posix = [x.as_posix() for x in polypaths]
        poly_idx = polypaths_posix.index(polypath.as_posix()) # not necessarily the same if we do not start at 0 in the list (for testing for example!!) 
        # the poly_idx is needed for the indexing of the different projection pathlists

        print(polypath)

        print('pidx={}'.format(pidx))
        print('poly_index (position in list): {}'.format(poly_idx))

        # measure time of one polypath
        start_pidx = time.time()

        # x is longitude, y is latitude
        # load array for current polygon path. in Orthogonal projection
        try:
            # TODO: think about what do I really need from this?
            arr, geoms, dataset, masks, mask_ids, sparse_masks, holes_counts = load_tiff_arr_specified(polypath, geopath)
            # print('length geoms: {}'.format(len(geoms['features'])))
        except ValueError: # then, there is no saved polygon for this geopath. This might happen if the geotiff is just on the border, or if only an excerpt was mapped (17ESREGMAP02)
            print('WARNING: No polygon data for {}. I continue without this array.'.format(polypath))
            # still add area data?
            # yes! might also just be an 'empty' geojson!
            # linament AREA TO CSV
            save_empty_linea_area(polypath, savebasepath, geopath)            
            continue
        ulx, xres, xskew, uly, yskew, yres, lrx, lry = get_geotransform(dataset)
        meanres = abs(xres)

        # check if array is 3 dimensional. if so, take arr[0] (this is a 'dumb' measure I am taking now to get rid of the alpha channel I've exported accidentally form the Lommel Seeliger...)
        if len(arr.shape) == 3:
            arr =  arr[0]
            print('WARNING! I DETECTED A THREE DIMENSIONAL INPUT ARRAY AND WILL ONLY CONSIDER ARR[0].')

        # mask geotiffs simply by multiplying with masks
        # https://stackoverflow.com/questions/56891206/numpy-multiplication-over-axis
        # use broadcasting by adding a dimension to arr
        extrac = arr[:,:,None] * masks
        extrac[extrac<0] = 0 # get rid of fill values!

        # indices to delete:
        del_idcs = []

        # check if geoms['features'] and masks.shape[-1] are the same
        if len(geoms['features']) != masks.shape[-1]:
            raise Warning('lenght of geoms is {}, but shape of masks is {}'.format(len(geoms['features']), masks.shape))


        # go through geoms, find gaps with minimum area rectangle boxes
        for i, feature in enumerate(geoms['features']):
            print(feature)
            subpoly_list = []
            # if there is only one subpolygon, skip
            if len(feature['geometry']['coordinates']) == 1:
                continue
            for fts_idx in range(len(feature['geometry']['coordinates'])):
                print(fts_idx) # sub-polygon nr.
                # for each sub-polygon, transform to a shapely polygon object
                subpoly = Polygon(feature['geometry']['coordinates'][fts_idx][0])
                # append to list
                subpoly_list.append(subpoly)
            # loop through all subpolygons to search for the closest other subpolygon in the list
            for poly1 in subpoly_list:
                for poly2 in subpoly_list:
                    if poly1 == poly2:
                        continue
                    print(poly1)
                    print(poly2)
                    print('')

            # from that, we can generate a ccbox (cross-cutting box)

from itertools import combinations
from shapely.geometry import box, MultiPolygon
from shapely.ops import nearest_points

def get_closest_subpolygons_box(poly1, poly2):
    """
    Constructs a bounding box between the two most adjacent subpolygons.
    """
    if len(multipoly.geoms) < 2:
        raise ValueError("MultiPolygon must contain at least two subpolygons.")

    # Find the closest pair of subpolygons using all combinations
    min_dist = float('inf')
    closest_pair = None

    for poly1, poly2 in combinations(multipoly.geoms, 2):
        # Calculate the distance between the two polygons
        dist = poly1.distance(poly2)
        if dist < min_dist:
            min_dist = dist
            closest_pair = (poly1, poly2)

    # Extract the actual closest points on the boundaries of these two subpolygons
    poly1, poly2 = closest_pair
    p1, p2 = nearest_points(poly1, poly2)

    # Construct a rectangular box between the two closest point coordinates
    # We use p1.x, p1.y, p2.x, p2.y to create the envelope
    min_x = min(p1.x, p2.x)
    max_x = max(p1.x, p2.x)
    min_y = min(p1.y, p2.y)
    max_y = max(p1.y, p2.y)

    # Add a minimal buffer so the box has a non-zero area (optional)
    if min_x == max_x: max_x += 0.0001
    if min_y == max_y: max_y += 0.0001

    return box(min_x, min_y, max_x, max_y)


from shapely.geometry import Polygon

# Define your polygon coordinates (e.g., L-shape)
coords = [(0, 0), (5, 0), (5, 2), (2, 2), (2, 5), (0, 5)]
poly = Polygon(coords)

# Get the minimum rotated bounding box
min_bbox = poly.minimum_rotated_rectangle

# Output the new bounding box coordinates
print(list(min_bbox.exterior.coords))

            centerpnts, (p01, p03), alpha = cv2.minAreaRect(ftl_idcs)


        # go through masks, every mask contains exactly one feature
        for mask_idx in range(masks.shape[-1]): # range(486,masks.shape[-1]): # masks.shape[-1]

            ######### STEP 1: pre-processing of mask
            # also, we need the current mask and the bounding box
            # extract the current mask
            m = masks[:,:,mask_idx] 
            # get bounding box and cutout (not minimum area bbox, but simply cutoff of black lines to reduce storage)
            p,k,w,h, cut, extrac_cut, p_center, k_center = cut_bounding_box(m, mask_idx, extrac)

            # test if mask is not almost empty. This might happen if a lineament got dissected during tessellation
            # for this, we simply calculate the area in px
            areapx = len(cut[cut>0])
            if areapx < area_limit:
                print('empty mask found. deleting and skipping. except for area calculation.')
                # add to list of indices to delete
                del_idcs.append(mask_idx)
                continue # then, we skip 

            ############ STEP 2: construct bounding box for fragments


        # now, we can savely delete all features that are stored in del_idcs list
        # ATTENTION: you cannot delete iteratively, as this will shift indices!!! We find a solution: delete biggest indices first --> this will not shift the earlier indices
        del_idcs.reverse()
        for delidx in del_idcs:
            del geoms['features'][delidx]



# %%
