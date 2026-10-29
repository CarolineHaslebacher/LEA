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
from shapely.geometry import Polygon, box,  MultiPolygon, LineString
from itertools import combinations
from shapely.ops import nearest_points
from shapely.plotting import plot_polygon

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


        # go through each feature, find gaps with minimum area rectangle boxes
        # 'subpoly' means sub-polygons of one feature. This might be the cross-cutting we are looking for!
        # it is not straight-forward to 'find the gap', because we do not know where the lineament 'is going'  
        for i, feature in enumerate(geoms['features']):
        for i, feature in enumerate(geoms['features'][:1]):
            print(feature)
            ccboxes_poly_list = []
            subpoly_list = []
            del_idcs_ccbox = []
            # if there is only one subpolygon, skip. There are no cross-cutting features recorded (does not mean that there aren't any though).
            if len(feature['geometry']['coordinates']) == 1:
                continue

            for fts_idx in range(len(feature['geometry']['coordinates'])):
                print(fts_idx) # sub-polygon nr.
                # for each sub-polygon, transform to a shapely polygon object
                # Hint: you can plot this polygon by executing 'subpoly'
                # or with: 
                # fig, ax = plt.subplots()
                # plot_polygon(subpoly, ax=ax, color='lightblue', alpha=0.7)
                subpoly = Polygon(feature['geometry']['coordinates'][fts_idx][0])
                # append to list
                subpoly_list.append(subpoly)

            # loop through all subpolygons to search for the closest other subpolygon in the list
            fig, ax = plt.subplots()
            plot_polygon(subpoly, ax=ax, color='darkblue', alpha=0.7)
            for poly1 in subpoly_list:
                ccboxes_poly = get_closest_subpolygons_box(subpoly_list, poly1)
                print(len(ccboxes_poly))

                # plot during development
                plot_polygon(poly1, ax=ax, color='lightblue', alpha=0.7)
                cci = 0
                for ccbox_poly in ccboxes_poly:
                    print(cci)
                    plot_polygon(ccbox_poly, ax=ax, color='green', alpha=0.3)
                    cci += 1
                    ccboxes_poly_list.append(ccbox_poly)

            # delete boxes that overlap to a large extent with (or contain) parts of a polygon
            for ccbox_idx, ccbox in enumerate(ccboxes_poly_list):
                for subpoly in subpoly_list:
                    if ccbox.intersects(subpoly):
                        # Get the overlapping geometry
                        inter_geom = ccbox.intersection(subpoly)
                        # union_geom = ccbox.union(subpoly)

                        # Calculate the intersection-over-feature_area (IoU) area
                        Iofa = inter_geom.area / subpoly.area
                        if Iofa > 0.1: # if the feature and the box overlap for more than 10% of the feature, we discard the box
                            if not ccbox_idx in del_idcs_ccbox:
                                del_idcs_ccbox.append(ccbox_idx)
                for other_ccbox_idx, other_ccbox in enumerate(ccboxes_poly_list):
                    # skip if they are equal
                    if not ccbox_idx == other_ccbox_idx:
                        # check if they overlap too much
                        intersection_ccboxes = ccbox.intersection(other_ccbox).area
                        union_ccboxes = ccbox.union(other_ccbox).area
                        IoU = intersection_ccboxes / union_ccboxes
                        if IoU > 0.5: # I expect such boxes to overlap to a large degree, let's say 50%
                            # if we need to remove one, check one of the two respective indexes is already listed
                            if not (ccbox_idx in del_idcs_ccbox or other_ccbox_idx in del_idcs_ccbox):
                                # if not, append this index
                                del_idcs_ccbox.append(ccbox_idx)

            # now, we can savely delete all features that are stored in del_idcs list
            # ATTENTION: you cannot delete iteratively, as this will shift indices!!! We find a solution: delete biggest indices first --> this will not shift the earlier indices
            # and we only want to have unique values

            del_idcs_ccbox.reverse()
            for delidx in del_idcs_ccbox:
                del ccboxes_poly_list[delidx]

            # plot for debugging:
            for ccbox_poly in ccboxes_poly_list:
                fig, ax = plt.subplots()
                for spoly in subpoly_list:
                    plot_polygon(spoly, ax=ax, color='blue', alpha=0.3)
                print(ccbox_poly)
                plot_polygon(ccbox_poly, ax=ax, color='green', alpha=0.3)
                plt.show()
            print(len(ccboxes_poly_list))


            # # now, let's also remove overlapping 'duplicate' boxes:
            # print(len(ccboxes_poly_list)) # 36
            # unique_ccboxes_poly_list = []
            # for poly in ccboxes_poly_list:
            #     # Check if the polygon topologically matches any already added
            #     if not any(poly.equals(u_poly) for u_poly in unique_ccboxes_poly_list):
            #         unique_ccboxes_poly_list.append(poly)
            #     # also test whether any IoU of ccboxes with themselves is overlapping too much 
            #     intersection = poly.intersection(subpoly)
                

            # print(len(unique_ccboxes_poly_list)) # 9
            # # plot for debugging:
            # for ccbox_poly in unique_ccboxes_poly_list:
            #     fig, ax = plt.subplots()
            #     for spoly in subpoly_list:
            #         plot_polygon(spoly, ax=ax, color='blue', alpha=0.3)
            #     print(ccbox_poly)
            #     plot_polygon(ccbox_poly, ax=ax, color='green', alpha=0.3)
            #     plt.show()

                        


            


TODO: try for i, feature in enumerate(geoms['features'][:1]): and debug missing gaps
and push to github

if this all works, I can continue with the search for overlaps

#%%
            # from that, we can generate a ccbox (cross-cutting box)
def get_min_width(poly):
    # Compute the minimum rotated rectangle
    min_rect = poly.minimum_rotated_rectangle

    # In Shapely 2.0+, use mbr.boundary.coords (or mbr.exterior.coords)
    coords = list(min_rect.boundary.coords)

    # Calculate the lengths of the adjacent edges forming the rectangle
    # The MBR has 5 coordinate pairs (the last one closes the polygon: index 0 to 4)
    width = LineString([coords[0], coords[1]]).length
    height = LineString([coords[1], coords[2]]).length

    # Ensure width is the shorter dimension and height is the longer dimension
    w = min(width, height)
    h = max(width, height)

    return w

def get_closest_subpolygons_box(subpoly_list, poly1, width=None):

    # If I only search for one closest polygon, I might miss the one on the other side of the subpolygon
    # fix: search for the two closest. Eliminate overlapping ccboxes in the end
    pairs = [
    (p, poly1.distance(p))
    for p in subpoly_list
    if not p.equals(poly1)
    ]
    pairs_sorted = sorted(pairs, key=lambda item: item[1])

    gap_boxes = []

    # dev plotting
    fig, ax = plt.subplots()
    plot_polygon(poly1, ax=ax, color='darkblue', alpha=0.7)

    for poly2, pdist in pairs_sorted[:2]:
        # print(poly2)
        # plot_polygon(poly2, ax=ax, color='lightblue', alpha=0.7)
        p1, p2 = nearest_points(poly1, poly2)

        dx, dy = p2.x - p1.x, p2.y - p1.y
        length = math.hypot(dx, dy)
        if length == 0:
            raise ValueError("Polygons touch or overlap; no gap to fill.")

        ux, uy = dx / length, dy / length
        nx, ny = -uy, ux

        if width is None:
            width = np.mean([get_min_width(poly1), get_min_width(poly2)])
            # width = length * 0.5
        half = width / 2

        # Extend slightly past p1/p2 so the difference cleanly trims to the edges
        eps = length * 0.01
        a = (p1.x - ux * eps, p1.y - uy * eps)
        b = (p2.x + ux * eps, p2.y + uy * eps)

        box_poly = Polygon([
            (a[0] + nx * half, a[1] + ny * half),
            (b[0] + nx * half, b[1] + ny * half),
            (b[0] - nx * half, b[1] - ny * half),
            (a[0] - nx * half, a[1] - ny * half),
        ])

        # Trim away anything overlapping the two polygons
        gap_box = box_poly.difference(poly1).difference(poly2)
        oriented = gap_box.minimum_rotated_rectangle
        gap_boxes.append(oriented)

    return gap_boxes

#%%


#%%

# from shapely.geometry import Polygon

# # Define your polygon coordinates (e.g., L-shape)
# coords = [(0, 0), (5, 0), (5, 2), (2, 2), (2, 5), (0, 5)]
# poly = Polygon(coords)

# # Get the minimum rotated bounding box
# min_bbox = poly.minimum_rotated_rectangle

# # Output the new bounding box coordinates
# print(list(min_bbox.exterior.coords))

#             centerpnts, (p01, p03), alpha = cv2.minAreaRect(ftl_idcs)


#         # go through masks, every mask contains exactly one feature
#         for mask_idx in range(masks.shape[-1]): # range(486,masks.shape[-1]): # masks.shape[-1]

#             ######### STEP 1: pre-processing of mask
#             # also, we need the current mask and the bounding box
#             # extract the current mask
#             m = masks[:,:,mask_idx] 
#             # get bounding box and cutout (not minimum area bbox, but simply cutoff of black lines to reduce storage)
#             p,k,w,h, cut, extrac_cut, p_center, k_center = cut_bounding_box(m, mask_idx, extrac)

#             # test if mask is not almost empty. This might happen if a lineament got dissected during tessellation
#             # for this, we simply calculate the area in px
#             areapx = len(cut[cut>0])
#             if areapx < area_limit:
#                 print('empty mask found. deleting and skipping. except for area calculation.')
#                 # add to list of indices to delete
#                 del_idcs.append(mask_idx)
#                 continue # then, we skip 

#             ############ STEP 2: construct bounding box for fragments


#         # now, we can savely delete all features that are stored in del_idcs list
#         # ATTENTION: you cannot delete iteratively, as this will shift indices!!! We find a solution: delete biggest indices first --> this will not shift the earlier indices
#         del_idcs.reverse()
#         for delidx in del_idcs:
#             del geoms['features'][delidx]



# %%
