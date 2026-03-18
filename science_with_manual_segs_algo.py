# 2022-05-11
# Science with manual segmentations
# See "science_with_manual_segmentations" powerPoint
# Caroline Haslebacher

#%% import

import os
# from winreg import HKEY_PERFORMANCE_DATA
from geojson.geometry import Point
from osgeo import gdal, gdalnumeric, ogr
from pathlib import Path
import matplotlib.pyplot as plt

import json
import skimage.draw
import numpy as np
import skimage.io

import cv2
from PIL import Image, ImageChops

import skimage.measure

from geojson import MultiPolygon, Feature, FeatureCollection, dump, Point
import math
import scipy
from scipy.sparse import coo_matrix

import random

from datetime import datetime
import time


import scipy
import scipy.signal as signal
from scipy.spatial.distance import cdist

import xml.etree.ElementTree as ET
import json
import pandas as pd
import skimage.measure
import pickle

import argparse
from MAIN_loop1_to_loop4 import get_polygeopaths

# %matplotlib inline

#%%

#%%
# copy functions from preparing_Galileo_tiles_from_qgis.py
def xcoord_to_px(xarr, ulx, xres):
    return (( xarr - ulx) / xres).astype('int64')

def ycoord_to_px(yarr, uly, yres):
    return (( yarr - uly) / yres).astype('int64')

def px_to_xcoord(xpx, ulx, xres):
    return (xpx*xres + ulx)

def px_to_ycoord(ypx, uly, yres):
    return (ypx*yres + uly)

def coordm_to_lon(coordm): # coordinate in metres from qgis
    radius = 1560800 # metres (corrected on 2023-12-06)
    
    if coordm < 0:
        lon = coordm/(np.pi*radius) * 180
        # For longitudes to the left of the 0 median in qgis:Take absolute value and add 180°​
        lon = abs(lon)+180
    else:
        lon = coordm/(np.pi*radius) * 180
        # then, in QGIS, it goes from the median to the right instead of to the left, like in the usgs map
        lon = 180-lon
    return lon

def coordm_to_lat(coordm):
    radius = 1560800 # metres (corrected on 2023-12-06)
    lat = coordm/(np.pi*radius)*180
    return lat

def mercator_coordm_to_lat(coordm):
    '''
    input: the latitude coordinates in meters from Mercator projection
    output: the latitude coordinates in degrees (0-90°)
    '''
    radius = 1560800 # metres
    lat = np.arctan(np.sinh( coordm / radius))
    return lat*180/np.pi # to degrees, not radians

def mercator_coordm_to_lon(coordm): # coordinate in metres from qgis
    '''
    input: the longitude coordinates in meters from Mercator projection
    output: the longitude coordinates in degrees (0-90°)
    '''
    radius = 1560800 # metres

    # I think I can use the same if else conditional statement than for equirectangular to account for east/west-positive
    # NO: in contrary to equirectangular, the mercator projection I selected here is 
    if coordm < 0:
        lon = coordm/(np.pi*radius) * 180
        # In contrary to 
        lon = abs(lon)
    else:
        lon = coordm/(np.pi*radius) * 180
        lon = 360-lon
    return lon

def ortho_coordm_to_lon_lat(x, y, lambda_0, phi_1): 
    '''
    input: 
    - longitude (x) in meters, from Orthographic projection
    - latitude (y) in meters, from Orthographic projection
    - lambda_0: center longitude of projection, in radians
    - phi_1: center latitutde of projection, in radians
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

def equi_coord_deg_to_coordm(lond, latd, phi0=0):
    '''
    input:
    - lond: longitude in degrees (0-360E)
    - latd: latitude in degrees (-90 to 90 (N to S))
    - phi0 is the central meridian in degrees (0 - 360deg)

    output:
    - longitude in meters
    - latitude in meters

    fixed:
    - radius
    - center_lon = 0 deg
    '''
    radius = 1560800 # metres
    latm = radius* np.pi/180*latd
    lonm = radius* np.pi/180*(lond - phi0)
    return(lonm, latm)



def get_lon_lat(coordlon, coordlat):
    return coordm_to_lon(coordlon), coordm_to_lat(coordlat)

# 180def calc_stress_field(lon, lat):
    # this function calculates the stress in Pa for a given longitude and latitude


def get_geotransform(dataset):
    # from here: https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
    # ulx, uly is the upper left corner, lrx, lry is the lower right corner
    ulx, xres, xskew, uly, yskew, yres  = dataset.GetGeoTransform()
    lrx = ulx + (dataset.RasterXSize * xres)
    lry = uly + (dataset.RasterYSize * yres)

    return ulx, xres, xskew, uly, yskew, yres, lrx, lry

def geojson_to_mask(geoms, dataset):
    # get coordinates
    ulx, xres, xskew, uly, yskew, yres, lrx, lry = get_geotransform(dataset)

    # intialize numpy array which gets filled
    npolyg = len(geoms['features'])
    mask = np.zeros([dataset.RasterYSize, dataset.RasterXSize, npolyg], dtype=np.uint8)
    mask_ids = []
    holes_counts_list = []
    del_colls = []

    for i, feature in enumerate(geoms['features']):
        # sometimes, a feature consists of multiple polygons! We loop through with fts_idx
        # print(i)
        # print(feature['properties']['id_int'])
        # check if this is really a Multipolygon and not a LineString!
        if feature['geometry']['type'] == 'LineString':
            print('LINESTRING FOUND. I IGNORE.')
            continue
        elif feature['geometry']['type'] == 'MultiLineString':
            print('MultiLINESTRING FOUND. I IGNORE.')
            continue
        elif feature['geometry']['type'] == 'Point':
            print('POINT FOUND. I IGNORE.')
            continue
        elif feature['geometry']['type'] == 'GeometryCollection':
            print('GeometryCollection FOUND. I investigate what is inside.')
            count_newfeat = 0
            for collectionfeature in feature['geometry']['geometries']:
                # print(collectionfeature)
                if collectionfeature['type'] != 'LineString': # other unknown types get caught later
                    # append to the list
                    newfeat = Feature()
                    # inherit properties:
                    newfeat['properties'] = feature['properties'].copy()
                    newfeat['geometry'] = collectionfeature 
                    geoms['features'].append(newfeat)
                    # print(newfeat)
                    # increase count
                    count_newfeat = count_newfeat + 1
            # delete GeometryCollection now:
            # check if we are really deleting the geometrycollection
            if feature['geometry']['type'] == 'GeometryCollection':
                # print('adding Geometrycollection to the to-be-deleted list: {}'.format(geoms['features'][i]))
                del_colls.append(i) # I do not delete here, as this would lead to chaos in the loop. I delete further down. I don't need to worry because I append and then CONTINUE for each geometrycollection found
            # We need to add entries to the mask. For each added feature, we add one. (The geometrycollection is already counted for in npolyg)
            # since empty masks get filtered out later anyway, we don't care if there are empty masks
            # we need to make the mask of shape (size, size, npoly) bigger by count_newfeat-1 (minus one, because we delete the GeometryCollection)
            added_mask = np.zeros([dataset.RasterYSize, dataset.RasterXSize, count_newfeat], dtype=np.uint8)  
            # print(mask.shape)
            # concatenate the existing and the new mask along the npolyg dimension
            mask = np.concatenate([mask, added_mask], axis=2)
            # print(mask.shape)
            continue # now, we can continue and loop through newly added features (right? --> yes! tested. If you append in a loop through the list, you can have an infinite loop!)

        elif feature['geometry']['type'] != 'Polygon' and feature['geometry']['type'] != 'MultiPolygon':
            raise Warning('unknown type found: {}'.format(feature['geometry']['type']))

        # print(feature['geometry']['coordinates'])
        # if i in [30, 35, 36, 41, 86]:
        #     continue
        holes_count = 0 # count number of holes!
        for fts_idx in range(len(feature['geometry']['coordinates'])):
            holes = False # set holes index to False for each part
            
            p = feature['geometry']['coordinates'][fts_idx]
            if len(p) == 1:
                p = p[0] # then, we need to reduce the dimension
                # now, len(p) are the number of points
            # Handling holes (from conversion of shapefiles to geojson)
            # if len(p) still 2 or 3, then, this means it has 1 or 2 holes.
            # we can better check the shape of parr

            # convert to array
            parr = np.asarray(p, dtype=object)
            # HOLES: in the case of holes, parr.shape is for example (3,) instead of e.g. (5,2).
            # also, numpy throws a warning: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray.
            if len(parr.shape) == 1:
                # change index
                holes = True
                # consider the holes count
                holes_count += parr.shape[0] - 1 # e.g. 3 parts, means 2 holes (1 is the outer part)
                # the first entry is 
                parr = np.asarray(p[0])               

            # convert coordinates to pixels
            # p_conv = [(p[idx][0] - uly)/yres for idx in range(0, 1)]

            parrx = xcoord_to_px(parr[:,0], ulx, xres)
            parry = ycoord_to_px(parr[:,1], uly, yres)

            # Get indexes of pixels inside the polygon and set them to 1
            # mask_pol = skimage.draw.polygon2mask(mask.shape, parrxy)
            rr, cc = skimage.draw.polygon(parry, parrx, mask.shape)
            mask[rr, cc, i] = 1

            if holes:
                # then, we subtract all holes by setting the mask to 0 in these points
                for pi in range(1, len(p)):
                    # loop through each hole
                    parr = np.asarray(p[pi])               

                    parrx = xcoord_to_px(parr[:,0], ulx, xres)
                    parry = ycoord_to_px(parr[:,1], uly, yres)

                    # Get indexes of pixels belonging to the hole and set them to 0$
                    # ATTENTION: this way, holes to not count as a disruption, as they are not counted as a 'part'
                    rr, cc = skimage.draw.polygon(parry, parrx, mask.shape)
                    mask[rr, cc, i] = 0

        # we need to save the class_id for every mask number
        mask_ids.append(int(feature['properties']['id_int']))
        # per feature, we add the total number of holes to the list
        holes_counts_list.append(holes_count)

    # delete now geometryCollections        
    del_colls.reverse()
    # print(mask.shape)
    for delidx in del_colls:
        del geoms['features'][delidx]
        mask = np.delete(mask, delidx, axis=-1)
    # print(mask.shape)

    # and delete corresponding masks!

        # # tested with 'G7ESTYRMAC02'
        # fig = plt.figure(frameon=False, figsize=(10,10))

        # plt.imshow(arr+100, cmap=plt.cm.gray, interpolation='nearest')
        # for maskt in range(mask.shape[2]):
        #     plt.imshow(mask[:,:,maskt], cmap=plt.cm.gray, alpha=.08, interpolation='bilinear')

        # plt.show()

    return mask, mask_ids, holes_counts_list, geoms # added return of geoms, because we might have changed this here!

def load_tiff_arr(polypath, geopaths):
    # this function loads an array for a certain polygon path

    # make sure we get really the right geopath index (.tif file from .geojson)
    for gidx, geopath in enumerate(geopaths):
        if '_'.join(polypath.stem.split('_')[:-1]) in str(geopath):
            pidx = gidx

    # open geojson
    with open(polypath) as jsf:
        geoms = json.load(jsf)

    # open geotiff
    dataset = gdal.Open(geopaths[pidx].as_posix(), gdal.GA_ReadOnly)  #, gdal.GA_ReadOnly not to change the file!

    # get projection information
    # proj = dataset.GetProjection()
    # transform to numpy array
    arr = dataset.ReadAsArray()
    arr[arr<0] = 0 # get rid of fill values here
    # plt.imshow(arr, vmin=0)
    # plt.show()

    # fill mask
    masks, mask_ids, holes_counts, geoms = geojson_to_mask(geoms, dataset)
    # save memory:
    sparse_masks = []
    for mask_idx in range(masks.shape[-1]): # for mask_idx in mask layers (e.g. 18 layers)
        # since most of the mask array elements are zero, we can make use of sparse (https://pythonspeed.com/articles/numpy-memory-footprint/)
        sparse_masks.append(scipy.sparse.coo_matrix(masks[:,:,mask_idx])) ## .row and .col give the indices as an array
    # stack to array
    sparse_masks = np.stack(sparse_masks)

    return arr, geoms, dataset, masks, mask_ids, sparse_masks, holes_counts

def load_tiff_arr_specified(polypath, geopath):
    '''
    # this function loads an array for a certain polygon path AND A GIVEN GEOPATH, IN CONTRARY TO LOAD_TIFF_aRR
    
    '''
    # open geojson
    with open(polypath) as jsf:
        geoms = json.load(jsf)

    # open geotiff
    print(geopath)
    dataset = gdal.Open(geopath.as_posix(), gdal.GA_ReadOnly)  #, gdal.GA_ReadOnly not to change the file!

    # get projection information
    # proj = dataset.GetProjection()
    # transform to numpy array
    arr = dataset.ReadAsArray()
    arr[arr<1e-10] = 0 # get rid of fill values here
    # plt.imshow(arr, vmin=0)
    # plt.show()

    # fill mask
    masks, mask_ids, holes_counts, geoms = geojson_to_mask(geoms, dataset)
    # save memory:
    sparse_masks = []
    for mask_idx in range(masks.shape[-1]): # for mask_idx in mask layers (e.g. 18 layers)
        # since most of the mask array elements are zero, we can make use of sparse (https://pythonspeed.com/articles/numpy-memory-footprint/)
        sparse_masks.append(scipy.sparse.coo_matrix(masks[:,:,mask_idx])) ## .row and .col give the indices as an array
    # stack to array
    sparse_masks = np.stack(sparse_masks)

    return arr, geoms, dataset, masks, mask_ids, sparse_masks, holes_counts


def aux_to_df(auxfile):
    tupleinfo = [] # list of tuples

    # metadata_file = ssi_path.joinpath(ssi_name + '.aux.xml')

    tree = ET.parse(auxfile)
    root = tree.getroot()
    metadata = root.find('Metadata')
    metadata_as_json = json.loads(metadata.text)

    # print(metadata_as_json['IsisCube']['Mapping'])

    # get Min Longitude
    minlon = metadata_as_json['IsisCube']['Mapping']['MinimumLongitude']
    # get Max Longitude
    maxlon = metadata_as_json['IsisCube']['Mapping']['MaximumLongitude']
    # get Min Latitude
    minlat = metadata_as_json['IsisCube']['Mapping']['MinimumLatitude']
    # get Max Latitude
    maxlat = metadata_as_json['IsisCube']['Mapping']['MaximumLatitude']

    # get resolution
    pixres = metadata_as_json['IsisCube']['Mapping']['PixelResolution']['value']

    tupleinfo.append((auxfile, auxfile.stem,  pixres, minlon, maxlon, minlat, maxlat))
    # print(len(tupleinfo))

        
    # construct pandas dataframe with column headings
    colnames = ['image_file', 'filename', 'pixel_resolution [m/px]', 'MinLon [deg]', 'MaxLon [deg]', 'MinLat [deg]', 'MaxLat [deg]']
    df_auxinfo = pd.DataFrame(tupleinfo, columns = colnames)

    return df_auxinfo


def mesh(shape0, shape1):
    # generates meshgrid
    xline = np.arange(0, shape0)
    yline = np.arange(0, shape1)
    x, y = np.meshgrid(xline,yline)
    return x,y

def mask0to90(alpha, shape0, shape1, width): # generates mask for angles 0 to 90 deg
    x,y = mesh(shape0, shape1)
    # take int() of shape1/2, otherwise we get an empty mask for 0°
    mask = (y - int(shape1/2) < np.tan(alpha*np.pi/180)*(x- int(shape1/2)) + width/(2*np.cos(alpha*np.pi/180))) & (y- int(shape1/2) > np.tan(alpha*np.pi/180)*(x- int(shape1/2)) - width/(2*np.cos(alpha*np.pi/180)))
    return mask

def maskov90(alpha, shape0, shape1, width): # generates mask for angles above 90 deg
    x,y = mesh(shape0, shape1)
    mask = (y - int(shape1/2) > -np.tan((180-alpha)*np.pi/180)*(x- int(shape1/2)) + width/(2*np.cos(alpha*np.pi/180))) & (y- int(shape1/2) < -np.tan((180-alpha)*np.pi/180)*(x- int(shape1/2)) - width/(2*np.cos(alpha*np.pi/180)))
    return mask

def mask_dict(angles, shape0, shape1, width):
    # this function returns a dictionary with one kernel for each angle (180 degrees divided by the number of angles)
    # with these kernel, the main direction of one lineament via convolution can be determined

    # initialise dictionary
    di_angle_idc = {angle: None for angle in angles}

    for alpha in angles:
        mask = mask0to90(alpha, shape0, shape1, width)
        if alpha > 90:
            mask = maskov90(alpha, shape0, shape1, width)
        # append mask to dict
        di_angle_idc[alpha] = mask

    return di_angle_idc

def maindir(img, angles, shape0, shape1, width):
    # shape1 is the shape[1] of the kernel!
    # detection of maximum output of convolution of kernel with image
    # the idea is that the maximum indicates the main direction
    # initialise minimum so that first value is for sure below min
    minimun = img.shape[0]*img.shape[1]

    # prepare dictionary with output values for further investigation
    outdict = {angle: {'conv': None, 'nonzero': None} for angle in angles}

    di_angle_idc = mask_dict(angles, shape0, shape1, width)

    maxangle=None
    for angle in angles:
        output = signal.convolve2d(img, di_angle_idc[angle], mode='same')
        # append to output dictionary
        outdict[angle]['conv'] = output
        # for debugging:
        # plt.imshow(output)
        # plt.show()
        # angle of main direction is where image is the least blurred --> count pixels above 0
        nonz = len(output[output>0])
        outdict[angle]['nonzero'] = nonz
        if nonz < minimun:
            minimun = nonz
            maxangle = angle
    if maxangle == None:
        print('maxangle is None, investigate.')

    return maxangle

# example use:
# test image 3
# mask_idx = 9
# tstimg = masks[:,:,mask_idx]
# plt.imshow(tstimg)
# size= 10
# width=int(size/3)
# di = np.zeros((size, size))

# nangles = 12 # with 18, we have one image per 10° for example. So this is our sensitivity
# angles = np.arange(0,180,180/nangles)

# maxangle = maindir(tstimg, angles, size, size)

def get_empty_class_dict(num_classes=6):
    return {i: [] for i in range(1, num_classes)}

def get_mean_std_of_class(extrac, geoms, save=True):
    # we can average the non.zero pixels 
    # or better, show the histogram
    means = []
    stds = []
    hists = [] 
    # bins are always the same
    for idx in range(extrac.shape[-1]):
        # exclude 0
        hist, bins = np.histogram(extrac[:,:,idx], range=(0.0001, 1), bins=20, density=True)
        mean = extrac[:,:,idx].mean()
        std = extrac[:,:,idx].std()
        means.append(mean)
        stds.append(std)
        hists.append(hist)

    # now sort into categories
    # len(geoms['features']) = len(means) and so
    num_classes = 6
    means_dict = get_empty_class_dict()
    stds_dict = get_empty_class_dict()
    hists_dict = get_empty_class_dict()
    # {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    # now we have a dict where we can sort in our means and hists, depending on which category it belongs to
    # produce a list of indices with corresponding feature category
    ls_cats = [feature['properties']['id_int'] for feature in geoms['features']]

    for idx in range(len(ls_cats)):
        cat = ls_cats[idx]
        means_dict[cat].append(means[idx])
        stds_dict[cat].append(stds[idx])
        hists_dict[cat].append(hists[idx])

    # average
    means_avg = {i: np.mean(means_dict[i]) for i in range(1, num_classes)}
    stds_avg = {i: np.mean(stds_dict[i]) for i in range(1, num_classes)}
    hists_avg = {i: np.mean(hists_dict[i], axis=0) for i in range(1, num_classes)}

    # save
    if save:
        # directory
        # current = os.getcwd()
        # titaniach = Path(current.split('Caroline')[0]) / 'Caroline' --> outdated. see lines below.
        current = Path.cwd()
        titaniach = current.parents[0]
        # basesavep = titaniach / 'lineament_detection/galileo_manual_segmentation/output' / 'diagrams' --> outdated. see lines below.
        basesavep = titaniach / 'output' / 'diagrams'
        os.makedirs(basesavep, exist_ok=True)

        # plot averaged histograms and means
        plt.errorbar(np.arange(1,num_classes), means_avg.values(), stds_avg.values(), marker='o', linestyle='')
        for a,b in enumerate(means_avg.values()): 
            plt.text(a+1.05, b, '{:.4f}'.format(b))
        plt.title('means of ' + '_'.join(polypath.stem.split('_')[:-1]))
        plt.savefig(basesavep.joinpath('_'.join(polypath.stem.split('_')[:-1]) + '_means.png'))
        plt.show()

        for cidx in range(1,num_classes): # leave out 0 category
            try:
                plt.plot(bins[1:], hists_avg[cidx], label=f'category {cidx}')
            except ValueError:
                print('Category is not in this image: {}'.format(cidx))

        plt.legend()
        plt.title('histograms of ' + '_'.join(polypath.stem.split('_')[:-1]))
        plt.savefig(basesavep.joinpath('_'.join(polypath.stem.split('_')[:-1]) + '_histograms.png'))
        plt.show()

    return means_dict, stds_dict, hists_dict, hists_avg

def gaussian_kernel(size, sigma=1):
    size = int(size) // 2 # // floor division (rounds down to nearest wholest number after dividing by 2)
    x, y = np.mgrid[-size:size+1, -size:size+1]
    if sigma > 0.00:
        normal = 1 / (2.0 * np.pi * sigma**2)
    else:
        print('sigma is zero! I will set it to 1.')
        sigma = 1
        normal = 1 / (2.0 * np.pi * sigma**2)
    g =  np.exp(-((x**2 + y**2) / (2.0*sigma**2))) * normal
    return g


# NOTE: this is currently a Baustelle! I can finish it, but don't need it until double ridge profiles..
# def get_profile(azimuth, cut, dn_masked_cut, screen_width=5, width_only=False): # dn_masked_cut is extrac_cut
#     '''
#         This function scans a lineament by defining a line perpendicular to its azimuth
#         It calculates a profile and calculates a width in pixels.
#         For this, it first rotates the image back to have an azimuth of 0 degrees and then scans through.
#         This implies a small error source from the interporlation during rotation.

#         # screen_width = 1 # this gives us a 1-line profile. we 'screen' the image
#         maxangle = azimuth
#         dn_masked_cut = extrac_cut
#         screen_width = round(lengthpx/20)
#     '''

#     # convert array to pillow image
#     # Rotate the image by azimuth counter clockwise
#     # center of rotation is center of image.
#     # we use 'expand' to expand the array so that everything fits in
#     im = Image.fromarray(255*cut) #np.uint8(cm.gist_earth(myarray)*255))
#     im_masked = Image.fromarray(dn_masked_cut) # dn_masked_cut is already uint8
#     # rotate
#     im_rotated = im.rotate(angle=azimuth, resample=Image.Resampling.BICUBIC, expand=True)
#     im_masked_rotated = im_masked.rotate(angle=azimuth, resample=Image.Resampling.BICUBIC, expand=True)
#     # back to numpy array:
#     cut_rotated = 1*np.array(im_rotated).astype('bool') # recast to boolean
#     dn_masked_cut_rotated = np.array(im_masked_rotated)
#     # cut unnecessary pixels:
#     _, _, _, _, cut_rotated, _, _, _ = cut_bounding_box(cut_rotated.astype('uint8'))
#     _, _, _, _, dn_masked_cut_rotated, _, _, _ = cut_bounding_box(dn_masked_cut_rotated.astype('uint8'))

#     # check if screen_width is non-zero:
#     if screen_width == 0:
#         screen_width = 1

#     # x_read is the indicator of the x-line that we read out (we take the profile of)
#     while x_read <= cut_rotated.shape:







#     # beta = alpha - 90°, to get direction of maximal change
#     # 'kernel' image for extracting the profile
#     beta = 90 #  (90 + azimuth)%180 # The modulo 180 helps to keep beta within 0 and 180 degrees.
#     # NOTE: if maxangle == 0 or == 180 deg, the tangens of 90 and 270 degrees is infinite. 
#     # Therefore, our y_lim would be infinite!
#     # we take care of this further down

#     di = np.zeros( cut.shape )
#     # construct meshgrid
#     yline = np.arange(0, di.shape[0])
#     xline = np.arange(0, di.shape[1])
#     x, y = np.meshgrid(xline,yline)
#     # x-shift for angles >0, <90,
#     # y-shift for angles >90 <=180, =0
#     # go through x-axis/ y-axis
#     # x- and y-shifts tell me how I have to shift the mask for each scan so that I go 'forwards'
#     # # NOTE: 2022-02: I commented the 'if beta >0...' lines below, because I realized that I only need 
#     # if beta > 0 and beta < 90:
#         # x_lim = abs(di.shape[1]/2+abs(np.tan(beta*np.pi/180)*di.shape[1]/2) + screen_width/(2*np.cos(beta*np.pi/180)))
#         # x_shifts = np.arange(-math.ceil(x_lim), math.floor(x_lim), screen_width) # or + math.ceil ??
#         # y_shifts = np.zeros(len(x_shifts))8
#     if round(beta,0) == 90 or round(beta,0) == 270: # NOTE: the round(beta,0) makes sure that decimal points are neglected. If for example, the azimuth was 4.3e-6, beta would not exactly be 90, but the round(4.3e-6,0) is 0.0.
#         # print('beta is 90 or 270')
#         y_shifts = np.arange(0, di.shape[1], screen_width)
#         x_shifts = np.zeros(len(y_shifts))
#     else:
#         y_lim = 35 # abs(di.shape[0]/2 + abs(np.tan(beta*np.pi/180)*di.shape[0]/2) + screen_width/(2*np.cos(beta*np.pi/180)))
#         y_shifts = np.arange(-math.ceil(y_lim), math.floor(y_lim), screen_width) # start, stop, space  # 
#         x_shifts = np.zeros(len(y_shifts))
#     # list of profiles along direction of biggest change (maxangle)
#     ls_profile = []
#     # plt.imshow(cut)
#     # plt.show()

#     for idx, (x_shift, y_shift) in enumerate(zip(x_shifts, y_shifts)):
#         # if idx%100 == 0:
#         #     print(idx)
#         # NOTE: I substitute screen_width/(2*np.cos(beta*np.pi/180))) with 1/(2*np.cos(beta*np.pi/180))) in order to only make slices of 1
#         # special case for beta=90
#         if beta > 90:
#             mask = (y - di.shape[0]/2 - y_shift  > -np.tan((180-beta)*np.pi/180)*(x- di.shape[1]/2 - x_shift) + 1/(2*np.cos(beta*np.pi/180))) & (y- di.shape[0]/2 - y_shift  < -np.tan((180-beta)*np.pi/180)*(x- di.shape[1]/2  - x_shift) - width/(2*np.cos(beta*np.pi/180)))
#         else:
#             mask = (y - di.shape[0]/2 - y_shift < np.tan(beta*np.pi/180)*(x- di.shape[1]/2 - x_shift) + 1/(2*np.cos(beta*np.pi/180))) & (y- di.shape[0]/2 - y_shift  > np.tan(beta*np.pi/180)*(x- di.shape[1]/2 - x_shift) - width/(2*np.cos(beta*np.pi/180)))
#         plt.imshow(mask + cut)
#         plt.show()
#         # append to profile
#         # NOTE: in order to account for a screen_width > 1, we must find a way to take the mean along the opposite mask direction...
#         # for now, I use only slices of width = 1
#         if width_only == False: # then, we also want the profile
#             profile = dn_masked_cut[mask][dn_masked_cut[mask]>0]

#         else: # else, we do the same, but with the mask instead of th extracted dn_masked_cut (because sometimes we don't have it)
#             profile = cut[mask][cut[mask]>0]
#         # for both cases, append only if the profile is not empty. This happens often because the x and y-shift are not limited correctly.
#         if len(profile) > 0:
#             ls_profile.append(profile)
        

#     # plot
#     # color development:
#     # from matplotlib.colors import to_hex
#     # collength = len(ls_profile)
#     # cols = [to_hex(plt.cm.viridis(i / collength)) for i in range(collength)] 

#     # calculate width by taking the median and ignoring empty values!
#     lwidth = np.median([len(arr) for arr in ls_profile])
#     if width_only:
#         return lwidth
#     # np.median([len(arr) for arr in ls_profile]) # --> this only works if screen_width is 1!
#     # get max of profile length
#     maxplen = max([len(arr) for arr in ls_profile])

#     # first, add nan to account for different widthspx
#     # then, average by ignoring nan values
#     # initialise new list
#     ls_cprofile = [] # 'c' stands for centered
#     for pidx, prof in enumerate(ls_profile):
#         # check how much is missing to 'full length'
#         difftomax = maxplen - len(prof)
#         if difftomax > 0:
#             # in this case, we need to add nan on both sides
#             # we initialise a new array filled with nan values
#             newarr = np.full(maxplen, np.nan)
#             addi = int(difftomax/2)
#             # we center (with a bias for even numbers) the profile
#             if difftomax%2 == 0: # if number is even
#                 newarr[addi:-addi] = prof
#             else: # for odd numbers of addi, addi rounds down, so we need to add on one side one place
#                 newarr[addi:-(addi+1)] = prof

#             # add newarr to list
#             ls_cprofile.append(newarr)
#         else:
#             # just add profile to list like it is
#             ls_cprofile.append(prof)

#     # stack to new array of shape (number of profiles, maxplen)
#     profarr = np.stack(ls_cprofile) # e.g. of shape (177, 10)

#     # now we can average over axis=0, ignoring nan
#     avgprof = np.nanmean(profarr, axis=0)
#     stdprof = np.nanstd(profarr, axis=0)

#     # plt.errorbar(np.arange(maxplen), avgprof, yerr=stdprof)
#     return avgprof, stdprof, maxplen, lwidth

# def get_profile_with_minAreaRect(maxangle, ftl_idcs ,cut, dn_masked_cut, screen_width=5):
#     '''
#         2023-11: new function based on shifting the points from the rotated rectangle to the right
#     '''
#     p0, p1, p2, p3 = cv2.boxPoints(cv2.minAreaRect(ftl_idcs))
#        # box points: array([[1351.0784 ,  943.363  ], # P0
#     #    [1995.6895 ,  727.33124], # P1
#     #    [2017.5178 ,  792.4644 ], # P2
#     #    [1372.9067 , 1008.49615]], dtype=float32) # P3 

#         # append to profile
#         profile = dn_masked_cut[mask][dn_masked_cut[mask]>0]
#         ls_profile.append(profile)

def get_segmental_lineament(azimuth, cut, screen_width=25, plot=False):
    '''
        This function provides a solution for the azimuth calculation of curved lineaments.
        As an output, you get masked versions of 'cut', with segmented lineaments. These can then be fed into an azimuth, length, width, and center_lon/lat analysis function.

        We need an approximate azimuth (=maxangle) for cutting, the cut (bounding box around the lineament mask, with the boolean lineament mask), and a screen_width (which says after how many pixels to tile up. This can be input dependend on the length in pixels, for example.)

        This function scans a lineament by defining a line perpendicular to its azimuth
        # screen_width = 25 # every 25 pixels, we make a mask
        maxangle = azimuth
        screen_width = round(lengthpx/4)
    '''
    
    # check if screen_width is non-zero:
    if screen_width == 0:
        screen_width = 1

    # convert array to pillow image
    # Rotate the image by azimuth counter clockwise
    # center of rotation is center of image.
    # we use 'expand' to expand the array so that everything fits in
    im = Image.fromarray(255*cut) #np.uint8(cm.gist_earth(myarray)*255))
    # toggle this line for extrac_cut usage: im_masked = Image.fromarray(dn_masked_cut) # dn_masked_cut is already uint8
    # rotate
    im_rotated = im.rotate(angle=azimuth, resample=Image.Resampling.BICUBIC, expand=True)
    # toggle this line for extrac_cut usage: im_masked_rotated = im_masked.rotate(angle=azimuth, resample=Image.Resampling.BICUBIC, expand=True)
    # back to numpy array:
    cut_rotated = 1*np.array(im_rotated).astype('bool') # recast to boolean
    # toggle this line for extrac_cut usage: dn_masked_cut_rotated = np.array(im_masked_rotated)
    # cut unnecessary pixels:
    _, _, _, _, cut_rotated, _, _, _ = cut_bounding_box(cut_rotated.astype('uint8'))
    # toggle this line for extrac_cut usage: _, _, _, _, dn_masked_cut_rotated, _, _, _ = cut_bounding_box(dn_masked_cut_rotated.astype('uint8'))

    # segment the lineament
    segments = []
    # x_read is the indicator of the x-line that we read out (we take the profile of)
    x_read = 0 # we take row 0 first
    while (x_read+screen_width) <= cut_rotated.shape[0]:
        # mask segment (by cutting mask) and append to list
        cutting_mask = np.zeros(cut_rotated.shape)
        cutting_mask[x_read:(x_read+screen_width), :] = 1
        segment = cut_rotated*cutting_mask
        segments.append(segment)
        # increment x_read by screen_width ('scan down')
        x_read = x_read + screen_width
    # for the last iteration:
    cutting_mask = np.zeros(cut_rotated.shape)
    cutting_mask[x_read:, :] = 1
    segment = cut_rotated*cutting_mask
    segments.append(segment)

    # NOTE: there will be a small overlap for most segments when backrotating.
    # rotate back!
    back_azi = 360-azimuth # the azimuth for backrotation (counter-clock-wise) is 360 deg minus the original azimuth
    backrots = []
    for segment in segments:
        if len(segment[segment>0]) == 0: 
            continue
        if plot:
            plt.imshow(segment)
            plt.show()
        # same rotational transformation as above
        im = Image.fromarray(255*segment)
        im_rotated = im.rotate(angle=back_azi, resample=Image.Resampling.BICUBIC, expand=True)
        segment_backrot =  1*np.array(im_rotated).astype('bool') # recast to boolean
        # check shape, adjust
        if segment_backrot.shape != cut.shape:
            # resize (again, a small negligible error source)
            # If the new array is larger than the original array, then the new array is filled with zeros.
            segment_backrot.resize(cut.shape)
        # append to list
        backrots.append(segment_backrot)
        if plot:
            plt.imshow(segment_backrot)
            plt.show()
        
        # _, _, _, _, segment_backrot, _, _, _ = cut_bounding_box(segment_backrot.astype('uint8'))
    # stack
    segments = np.stack(segments) # (num_tiles, size1, size2), e.g. shape (5, 204, 53)

    return segments


# BELOW IS THE OLD GET_SEGMENTAL_LINEAMENT FUNCTION SNIPPET 
    # # beta = alpha - 90°, to get direction of maximal change
    # # 'kernel' image for extracting the profile
    # beta = (90 + maxangle)%180 # The modulo 180 helps to keep beta within 0 and 180 degrees.
    # # NOTE: if maxangle == 0 or == 180 deg, the tangens of 90 and 270 degrees is infinite. 
    # # Therefore, our y_lim would be infinite!
    # # we take care of this further down

    # di = np.zeros( cut_merc.shape )
    # # construct meshgrid
    # yline = np.arange(0, di.shape[0])
    # xline = np.arange(0, di.shape[1])
    # x, y = np.meshgrid(xline,yline)
    # # x-shift for angles >0, <90,
    # # y-shift for angles >90 <=180, =0
    # # go through x-axis/ y-axis
    # # x- and y-shifts tell me how I have to shift the mask for each scan so that I go 'forwards'
    # # # NOTE: 2022-02: I commented the 'if beta >0...' lines below, because I realized that I only need 
    # # if beta > 0 and beta < 90:
    # #     x_lim = abs(di.shape[1]/2+abs(np.tan(beta*np.pi/180)*di.shape[1]/2) + screen_width/(2*np.cos(beta*np.pi/180)))
    # #     x_shifts = np.arange(-math.ceil(x_lim), math.floor(x_lim), screen_width) # or + math.ceil ??
    # #     y_shifts = np.zeros(len(x_shifts))
    # if round(beta,0) == 90 or round(beta,0) == 270: # NOTE: the round(beta,0) makes sure that decimal points are neglected. If for example, the azimuth was 4.3e-6, beta would not exactly be 90, but the round(4.3e-6,0) is 0.0.
    #     y_shifts = np.arange(0, di.shape[1], screen_width)
    #     x_shifts = np.zeros(len(y_shifts))
    # else:
    #     y_lim = abs(di.shape[0]/2 + abs(np.tan(beta*np.pi/180)*di.shape[0]/2) + screen_width/(2*np.cos(beta*np.pi/180)))
    #     y_shifts = np.arange(-math.ceil(y_lim), math.floor(y_lim), screen_width) # start, stop, space  # 
    #     x_shifts = np.zeros(len(y_shifts))

    # # list of profiles along direction of biggest change (maxangle)
    # ls_profile = []
    # if plot:
    #     plt.imshow(cut_merc)
    #     plt.show()

    # num_tiles = len(x_shifts) # this might lead to one totally empty tile, which can easily be sorted out!
    # for idx, (x_shift, y_shift) in enumerate(zip(x_shifts, y_shifts)):
    #     print(idx)
    #     if plot:
    #         plt.imshow(mask1)
    #         plt.show()
    #     # cast the 'old' mask from the previous loop to 'mask_past'
    #     # unless there is no previous loop, in which case we fill a mask of ones
    #     if idx > 0:
    #         mask_past = mask1
    #     else: # 
    #         mask_past = np.ones(cut_merc.shape, dtype='int')
    #     if plot:
    #         plt.imshow(~mask_past)
    #         plt.show()
    #     # if idx%100 == 0:
    #     #     print(idx)
    #     # NOTE: I substitute screen_width/(2*np.cos(beta*np.pi/180))) with 1/(2*np.cos(beta*np.pi/180))) in order to only make slices of 1
    #     # special case for beta=90
    #     if beta > 90:
    #         mask1 = (y - di.shape[0]/2 - y_shift  > -np.tan((180-beta)*np.pi/180)*(x- di.shape[1]/2 - x_shift) + 1/(2*np.cos(beta*np.pi/180)))
    #     else:
    #         # NOTE: I changed the sign from '<' to '>' in the next line! This way, the segments can be calculated the same as for beta > 90
    #         mask1 = (y - di.shape[0]/2 - y_shift > np.tan(beta*np.pi/180)*(x- di.shape[1]/2 - x_shift) + 1/(2*np.cos(beta*np.pi/180)))
    #     # calculate segment now, by multiplying the cut_merc (our initial mask) with the inverted mask1 times the mask_past, to only take one single segment
    #     segment = cut_merc*~mask1*mask_past
    #     if plot:
    #         plt.imshow(~mask1)
    #         plt.show()
    #         plt.imshow(mask_past)
    #         plt.show()
    #         plt.imshow(segment)
    #         plt.show()
    #     # multiply the cut_merc by the inverted mask1 and the past mask. Except for the last index, for which we also take the real mask1
        
    #     # only append if the segment is not totally empty
    #     if len(segment[segment>0]) > 0: 
    #         segments.append(segment)
    #     if idx == num_tiles-1:
    #         segment = cut_merc*mask1
    #         # again, filter if empty
    #         if len(segment[segment>0]) > 0: 
    #             segments.append(segment)
    #         if plot:
    #             plt.imshow(segment)
    #             plt.show()

    # # stack segments
    # segments = np.stack(segments) # (num_tiles, size1, size2), e.g. shape (5, 204, 53)

    # return segments



#%% 



def get_sub_polygeopaths(polypath_orig, geopath_orig):
    '''
        this function is slightly adapted from get_polygeopaths.
        It searches the correct subdirs where the tessellated files are.
        debug with:
        polypath_orig = poly_orig
        geopath_orig = geo_orig
    '''

    # os.makedirs(basepath / 'data/polygons/for_analysis', exist_ok=True)
    polypaths = sorted((polypath_orig).glob('*/*.geojson'))

    geopaths = []
    # construct geotiff from polypath
    for polypath in polypaths:
        polypath.parents[0]
        # get subdirectory, for example '17ESREGMAP01EXCERPT2_GalileoSSI_Equi-cog'
        subdir = polypath.parts[-2] # parts gives back a tuple, for examle: ('z:\\','Groups', ... 'Orthographic', 'geojson', '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed', '11ESCOLORS01-01_GalileoSSI_Equi-cog_fixed_19.geojson')

        # For this to work, the geotiff and the polygon need to have the exact same name, except for the line ending
        geopath = (geopath_orig / polypath.parts[-2]).joinpath('{}.tif'.format(polypath.stem))
        # append to list
        geopaths.append(geopath)

    return geopaths, polypaths

# 2024-01-26, from LineaMapper_to_img_sparsemasks_experimental.py,
# tested well!
def fit_line_to_mask(masktf, plot=False):
    '''
    input is a mask to fit (masktf), a NUMPY sparse tensor
    masktf = sparse_masks[mask_idx]
    '''
    # first, generate a 'contours'. list of numpy arrays
    # WATCH OUT: we have to swap x and y here! I checked/verified this function for many lineaments.
    ftl_idcs = [ np.array([[y, x]]) for x, y in  zip(masktf.row.tolist(), masktf.col.tolist() ) ]
    # stack to array of shape [Npoints, 1, 2] (2 for the x/y pair)
    ftl_idcs = np.stack(ftl_idcs)
    # fit the line with the off-the-shelve openCV function fitLine. L2 is simple euclidean distance.
    [vx,vy,x,y] = cv2.fitLine(ftl_idcs, cv2.DIST_L2,0,0.01,0.01)
    # calculate the angle
    alpha = np.arctan(vy/vx) # alpha is in radians
    # array([1.3424875], dtype=float32)
    # extract values!
    alpha = alpha[0]
    # calculate rms (developed in LineaMapper_to_img_sparsemasks_experimental.py)
    # 2024-02: assess curviness by goodness of fit. calculation of simple 2D geometry (distance of point to line)
    dist = (1 + vy**2 / vx**2 )**(-1) * ( (y - ftl_idcs[:,:,1])/vx + vy/(vx**2) * ( ftl_idcs[:,:,0] - x ) ) # w1 = (1+yv2/xv2)-1 * ((y0-py)/xv + yv(px-x0)/xv2)
    rms_dist = np.sqrt(np.sum(dist**2)/len(ftl_idcs)) # RMS = sqrt( sum(xi - xi_mean)**2 / N). len(ftl_idcs) only takes the first shape entry, which is the Number of Points
    if plot:
        tfimg = masktf.todense()
        rows,cols = tfimg.shape
        # [vx,vy,x,y] = cv2.fitLine(contours[cnt_idx], cv2.DIST_L2,0,0.01,0.01)
        lefty = int((-x*vy/vx) + y)
        righty = int(((cols-x)*vy/vx)+y)
        lineimg = cv2.line(tfimg*255,(cols-1,righty),(0,lefty),(255,255,0),2) # attention: I have to multiply tfimg by 255 to make it uint8!
        plt.imshow(lineimg)
        plt.show()
    return 90 + alpha*180/np.pi, rms_dist # change radians to degrees. Add 90 degrees to orient it with north = 0°

# the below was commented on 2024-01-26. I found that it was not reliable.
# # openCV fit a line
# # 2023-10-06, from LineaMapper_to_img_sparsemasks_experimental.py
def fit_line_to_mask_MINAREARECT(masktf):
    '''
    input is a mask to fit (masktf), a numpy (!) sparse COO matrix (https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_array.html#scipy.sparse.coo_array )
    output: the degree of rotation in degrees

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
    # attention: widht and length are not ordered by width and length, but rather P0-P3 comes second!
    # So: if P0-P3 is smaller than P0-P1, then it is the width, else the length.
    # Therefore, the below code is not necessary anymore!
    # test if P0-P3 is smaller than P0-P1
    if p03 <= p01: # then, P0-P3 is the width
        # If they are equal, We cannot really say if it is alpha or alpha + 90.
        # then, lw-width is minimal, meaning P0-P3 is the width
        beta = abs(alpha) + 90
    else:
        # then, P0-P3 is the length
        beta = 90 - abs(alpha)    

    # check directly here if P0-P3 is length, I calculate Beta = 90 - abs(alpha),
    # if P0-P3 is width, I calculate Beta = 90 + abs(alpha)
    # This then follows the Rhoden+2013 convention that north is 0°, East is 90°, south is 180°
    # Here, I give one example:
    # minAreaRect: ((1684.298095703125, 867.9136962890625), # Center x, y coords
    #         (68.693603515625, 679.8478393554688), # width, length
    #         71.4721908569336) # alpha
    # box points: array([[1351.0784 ,  943.363  ], # P0
    #    [1995.6895 ,  727.33124], # P1
    #    [2017.5178 ,  792.4644 ], # P2
    #    [1372.9067 , 1008.49615]], dtype=float32) # P3
    # The trick is to know how the points are numbered. See here: https://theailearner.com/tag/cv2-minarearect/ 
    # We calculate the length of the vector P0-P3. s = sqrt((-21.8)**2 + (-65.)**2) = 68. 
    # We compare this to our width and length and see that it matches the width much better.
    # we calculate the Beta angle from the alpha angle through Beta = abs(alpha) + 90.

    # # get points
    # p0, p1, p2, p3 = cv2.boxPoints(cv2.minAreaRect(ftl_idcs))
    # # calculate absolute value of vector in between points P0 and P3
    # p03_abs = p0 - p3
    # lw = np.sqrt( p03_abs[0]**2 + p03_abs[1]**2)
    # # compare to output of minAreaRect
    # if np.argmin((abs(lw-width), abs(lw-length))) == 0:
    #     # then, lw-width is minimal, meaning P0-P3 is the width
    #     beta = abs(alpha) + 90
    # else:
    #     # then, lw-length is minimal, meaning P0-P3 is the length
    #     beta = 90 - abs(alpha)

    # OTher option: fit a line
    # # fit the line with the off-the-shelve openCV function fitLine
    # [vx,vy,x,y] = cv2.fitLine(ftl_idcs, cv2.DIST_L2,0,0.01,0.01)
    # # calculate the angle
    # alpha = np.arctan(vy/vx) # alpha is in radians
    # Beta (might be flawed!), length, width, indices
    return beta, np.max((p01, p03)), np.min((p01, p03)), ftl_idcs # we also return the length and the indices, which we might want to use later

def nan_to_num(x, substitute=-1):
    if math.isnan(x) == True:
        x = substitute
    elif x == np.inf:
        x = substitute
    return x

def cut_bounding_box(m, mask_idx=None, extrac=None):
    '''
        This function cuts all empty black space away from a given boolean mask m.
        and returns a rectangular array containing precisely the mask.
        it also returns the bounding box parameters (ulx, uly, width, height).

        mask_idx and extrac are optional, but if extrac is given, mask_idx must also be given!
        This extracts
    '''
    # get bounding box
    p,k,w,h = cv2.boundingRect(m)
    # cut excerpt
    cut = m[k:k+h, p:p+w]
    if extrac is not None:
        extrac_cut = extrac[k:k+h, p:p+w, mask_idx]
    else: # cast to -1 if empty, so that it can be returned
        extrac_cut = -1
    # for lon/lat, just take center of bounding box as a proxy:
    # this can be referenced later in the script to the whole big image!
    p_center = p+int(w/2) # longitude
    k_center = k+int(h/2) # latitude

    return p,k,w,h, cut, extrac_cut, p_center, k_center

# def from_meta(auxfile, param):
#     '''
#     the parameter must be in the 'Mapping' group of the .aux file
#     From \\titania.unibe.ch\Space\Groups\PIG\Caroline\lineament_detection\galileo_manual_segmentation\code\manual_segmentation_image_info_table.py

#     example for auxfile path: ((titaniach / 'isis/data/galileo/usgs_photogrammetrically' / crs_folder).glob('*.aux.xml'))
#     '''
#     # auxfile = Z:\Groups\PIG\Caroline\isis\data\galileo\Europa_isis4\lbls\ssi_id.lbl

#     tree = ET.parse(auxfile)
#     root = tree.getroot()
#     metadata = root.find('Metadata')
#     metadata_as_json = json.loads(metadata.text)

#     # print(metadata_as_json['IsisCube']['Mapping'])

#     # get the parameter
#     metaparam = metadata_as_json['IsisCube']['Mapping'][param]

#     return metaparam

def load_metacsv(ssi_id):
    dfmeta = pd.read_csv(Path.cwd().joinpath('Galileo_SSI_Europa_labels.csv'))

    return dfmeta

def save_empty_linea_area(polypath, savebasepath, geopath):
    '''
    input:
    polypath, savebasepath, geopath, 
    '''
    save_csv_path_dt = savebasepath / '_'.join(polypath.stem.split('_')[:-1]) / 'csv_area_files'
    os.makedirs(save_csv_path_dt, exist_ok=True)
    # calculate the lineament area per class
    area_semantic_dict = {1: {'fraction': 0, 'km2': 0}, 
                    2: {'fraction': 0, 'km2': 0}, 
                    3: {'fraction': 0, 'km2': 0}, 
                    4: {'fraction': 0, 'km2': 0}, 
                    'total_linea': {'fraction': 0, 'km2': 0}, 
                    'total_area': {'fraction': 0, 'km2': 0}, # NOTE: this is the total area as a safeguard. the fraction should in the end add up to 1 (after merging, 'Merge_polygon_calculations.py')
                }
    # NO! this would not account for tessels that are half black, in border regions!! total_frame_area_px = masks.shape[0] * masks.shape[1]
    # load the array here!
    dataset = gdal.Open(geopath.as_posix(), gdal.GA_ReadOnly)  #, gdal.GA_ReadOnly not to change the file!
    ulx, xres, xskew, uly, yskew, yres, lrx, lry = get_geotransform(dataset) # we need xres and yres for total area calculation
    # transform to numpy array
    arr = dataset.ReadAsArray()
    arr[arr<1e-10] = 0 # get rid of fill values here
    total_frame_area_px = len(arr[arr>0])
    total_linea_area = 0
    # append total area and fraction
    total_frame_area_km2 = abs(total_frame_area_px*xres*yres*1e-6)
    area_semantic_dict['total_area']['km2'] = total_frame_area_km2
    area_semantic_dict['total_area']['fraction'] = 1.0

    # dict to pandas dataframe, to csv
    area_semantic_df = pd.DataFrame(area_semantic_dict)
    area_semantic_df.to_csv(save_csv_path_dt.joinpath(polypath.stem + '_linea_area.csv')) # Here, we need the full path, including the children _1 added!

    return

class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.float32):
            return float(obj)
        return json.JSONEncoder.default(self, obj)

#%%

if __name__ == '__main__': # this allows me to import functions defined above!

    '''
    This script reads in the tessellated orthographically projected geotiff and geojson and extracts useful lineament characteristics from it.
    This is then saved in basepath.parents[0] / 'output' as children (reprojected to equirectangular) and stitched back to parents later.
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
    parser.add_argument('--output_foldername', required=False,
                        default='Equirectangular',
                    metavar="the name of sub-folder to save",
                    help='e.g. Equirectangular')  
    parser.add_argument('--SSI_id', required=False,
                        default='false',
                    metavar="simply omit if not needed. Defines whether the input id (the first string when split by _) is an SSI ID (opposed to a mosaic)",
                    help='is a string, not a bool!')
    # instead of bool: action='store_true': This sets the argument to True if the flag is present on the command line and False if it is omitted. This is suitable for flags that enable a feature, where the default state is False.
    parser.add_argument('--illumination_info', action='store_true', 
                    help='use this flag if you want to attach illumination information from a Galileo SSI image of Europa.') 
    # somehow, the below lines were not used. Maybe I wanted to implement this and then did not.
    # parser.add_argument('--find_source_crs', required=False,
    #                     default='false',
    #                 metavar="simply omit if not needed. Defines whether we search for the original crs. If not, I assume it's Equirectangular.",
    #                 help='is a string, not a bool!')  

    args = parser.parse_args()

    basepath = Path(args.basepath) # titaniach / 'lineament_detection/galileo_manual_segmentation/data'

    equifile = args.crs_file
    savefolder = args.crsfoldername
    SSI_ID_toggle = args.SSI_id
    illum_info = args.illumination_info
    outputfoldername = args.output_foldername 
    '''
    test/Debug with:
    current = os.getcwd()
    titaniach = Path(current.split('Caroline')[0]) / 'Caroline'
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data_phoc_sunAzi'
    illum_info = True
    SSI_ID_toggle = 'false'
    equifile = 'Equirectangular_EUROPA.prj'
    savefolder = 'Equirectangular'
    outputfoldername = 'output/testing_scwmanualsegs_algo'


    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    basepath = titaniach / 'lineament_detection/global_maps_for_book/Europa_globalextraction/data'
    basepath = titaniach / 'lineament_detection/RegionalMaps_CH_NT_EJL_LP/mapping/data'
    basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    equifile = 'Equirectangular_EUROPA.prj'
    savefolder = 'Equirectangular'
    outputfoldername = 'output/testing_scwmanualsegs_algo'
    SSI_ID_toggle = 'false'
    # polypath = Path('//titania.unibe.ch/Space/Groups/PIG/Caroline/lineament_detection/RegionalMaps_CH_NT_EJL_LP/mapping/data/polygons/for_analysis/Orthographic/geojson/17ESREGMAP02chaos_GaliloSSI_Equi-cog/17ESREGMAP02chaos_GaliloSSI_Equi-cog_19.geojson')    
    polypath = Path('//titania.unibe.ch/Space/Groups/PIG/Caroline/lineament_detection/RegionalMaps_CH_NT_EJL_LP/mapping/data/polygons/for_analysis/Orthographic/geojson/17ESREGMAP02chaos_GaliloSSI_Equi-cog/17ESREGMAP02chaos_GaliloSSI_Equi-cog_20.geojson')

    # LineaMapper
    basepath = titaniach / 'publications_presentations/RegMaps_specialIssue_PSJ/full_RegMaps_predictions/do_analysis'
    crsfoldername = 'Equirectangular'
    crs_file = 'Equirectangular_EUROPA.prj'
    output_foldername = 'output/extraction_loop3'
    SSI_id = 'true'


    '''
    
    #%%
    # e.g.
    # basepath = titaniach / 'lineament_detection/galileo_manual_segmentation/data'
    # jsonpath = basepath / 'polygons/for_analysis/Equirectangular'
    # equifile = 'Equirectangular_EUROPA.prj'

    # equirectangular
    geo_orig = basepath / 'geotiffs/for_analysis'
    poly_orig = basepath / 'polygons/for_analysis' # this path gets redefined later!
    #####

    # specify the basepath (used for saving)
    savebasepath = basepath.parents[0] / outputfoldername # e.g. titaniach / 'lineament_detection/galileo_manual_segmentation/output/extraction_loop3'
    os.makedirs(savebasepath, exist_ok=True)
    equipath = poly_orig / savefolder
    ortho_prj_path = basepath / 'orthographic_projection'

    
    # define the correct class dict
    class_dict = {1: "band", 2: "double_ridge", 3: "ridge_complex", 4: "undifferentiated_linea"}

    # define parameters for main direction
    ######
    size= 10
    width=int(size/3)
    di = np.zeros((size, size))
    rms_limit_widthpx = 0.35
    area_limit = 20 # masks filled with fewer than XX pixels get ignored

    # nangles = 36 # with 18, we have one image per 10° for example. So this is our sensitivity
    # angles = np.arange(2,178,180/nangles)

    #####
    # connectivity=2 # 2-jump neighborhoood if connectivity=2

    ### NOTE: maybe not used any longer due to implementation of phocube
    # dict for north azimuth by mosaic (this is not automated because for each mosaic, I have to find the best matching individual SSI)
    # but can be automated for individual images
    north_azimuth = {
        'G7ESLOWFOT01': 84.08,
        'G7ESTYRMAC01': 285.068	,
        'G7ESTYRMAC02': 280.565,
        'G7ESAPEXCR01': 274.27,
        'G7ESAPEXCR04': 273.756,
        'G7ESAPEXCR05': 273.862,
        'E4ESMACSTR01': 90.437,
        'E6ESDRKLIN01': 93.386,
        '11ESCOLORS01-01': 91.129,
        '12ESFRTPLT01EXCERPT1': 91.342,
        '15ESREGMAP01EXCERPT1': 93.747,
        '17ESREGMAP03': 268.765,
        '25ESDARKBP01': 287.91,
        'G7ESAPEXCR02': 274.115,
        '17ESREGMAP01EXCERPT2': 95.438,
        '17ESREGMAP01EXCERPT1': 94.497,
        '11ESREGMAP01EXCERPT1': 89.254,
        '14ESWEDGES01EXCERPT1': 269.757
    }

    # the mosaics below are chosen manually and with the help of Z:\Groups\PIG\Caroline\isis\Europa_code\python_code\get_observation_with_id_clock_id.py
    mosaic_to_ssi_id = {
        'G7ESLOWFOT01': 'C0389767100R',
        'G7ESTYRMAC01': 'C0389772500R',
        'G7ESTYRMAC02': 'C0389773000R',
        'G7ESAPEXCR01': 'C0389778842R',
        'G7ESAPEXCR04': 'C0389780135R',
        'G7ESAPEXCR05': 'C0389780563R',
        'E4ESMACSTR01': 'C0374667300R',
        'E6ESDRKLIN01': 'C0426267400R',
        '11ESCOLORS01-01': 'C0420617239R',
        '12ESFRTPLT01EXCERPT1': 'C0426267200R',
        '15ESREGMAP01EXCERPT1': 'C0449961914R',
        '17ESREGMAP03': 'C0466677052R',
        '25ESDARKBP01': 'C0527275700R',
        'G7ESAPEXCR02': 'C0389779270R',
        '17ESREGMAP01EXCERPT2': 'C0466664352R',
        '17ESREGMAP01EXCERPT1': 'C0466664366R',
        '11ESREGMAP01EXCERPT1': 'C0420619278R',
        '14ESWEDGES01EXCERPT1': 'C0440955239R',
        '17ESREGMAP02testregion': 'C0466676865R',
        '17ESREGMAP02trainingregion': 'C0466676801R',
        '17ESREGMAP02fordb': 'C0466676801R',
        '17ESREGMAP02pr': 'C0466676801R',
        '17ESREGMAP02chaos': 'C0466676801R',
        '17ESREGMAP02forpub2024': 'C0466676801R',
        # '2025': 'C0484888726R', # this is a stereo. the other image is: 'C0449961826R'
        # 'manual': 'C0484888726R', # this is a stereo. the other image is: 'C0449961826R'
        # doing a better workaround by renaming with C0484888726R
    }


    # once and for all, we have only one projected file to import: the oblique orthographic projected geotiff and corresponding feature layer
    geo_orig = basepath / 'geotiffs/for_analysis/Orthographic'
    poly_orig = basepath / 'polygons/for_analysis/Orthographic/geojson'
    geopaths, polypaths = get_sub_polygeopaths(poly_orig, geo_orig) # THIS IS A NEW FUNCTION!

    ortho_proj_path = basepath / 'orthographic_projection'

    # geodict_keys = [polypath.stem for polypath in polypaths]
    # geodict = dict.fromkeys(geodict_keys) # dictionary with information for all mosaics!

    # measure time
    start_all = time.time()

    now = datetime.now()
    dt_string = now.strftime("%Y_%m_%d_%H_%M_")

    # # # or select just one polypath
    # for polypath in [polypaths[0]]: 
    #     print('you have selected only one specific polypath!')
    #%%
    # loop through polygons and geotiffs
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

        # NOTE: the below is not necessary anymore with the new script I have run on 2024-05-10 (\\titania.unibe.ch\Space\Groups\PIG\Caroline\lineament_detection\galileo_manual_segmentation\code\extraction_loops_pub1\shift_label_5_to_4.py) 
        # NOTE: for run with labels 1 - 5, including cusps, I implemented to shift id_int 5 to 4. and to delete (=don't save) id_int 4
        # if ('testregion' in polypath.as_posix()) or ('trainingregion' in polypath.as_posix()) or ('fordb' in polypath.as_posix()) or ('chaos' in polypath.as_posix()) or ('pr' in polypath.as_posix()): # for 17ESREGMAP02 testregion and trainingregion, I labeled 1-4 (no cusps)
        #     shift5to4 = False
        # else:
        #     shift5to4 = True

        if illum_info:
            # NOTE: galileopath is a helper variable that prevents loading a large array multiple times for the same parent polygon
            # therefore, we check whether we have a new polygon here:
            # but first, we need to check whether galileo_name has even been defined yet
            if 'galileo_name' in locals():
                if galileo_name != polypath.stem.split('_')[0]:
                    del galileopath # if it exists, we have to erase it for each new galileo_name
            # note that the sunazimuth extraction only works correctly in equirectangular projection
            galileo_name = polypath.stem.split('_')[0]
            # extract north azimuth and sun_illumination
            # get a corresponding SSI ID, from the mosaic name.
            if SSI_ID_toggle != 'false':
                ssi_id = galileo_name # e.g. 'C0449961800R'
            else: # default
                ssi_id = mosaic_to_ssi_id[galileo_name] # e.g. 17ESREGMAP02trainingregion
                # take 12 first characters for mosaic name (from 17ESREGMAP02trainingregion to 17ESREGMAP02)
                # however, this does not work for a mosaic such as 11ESCOLORS01-01 / 11ESCOLORS01-02
                # therefore, I ignore my special extents off the galileo name: ‘EXCERPT’ and ‘forpub2024’
                # and otherwise I leave it, because it was split already correctly with the line above (galileo_name = ..)
                if 'EXCERPT' in galileo_name:
                    galileo_name = galileo_name.split('EXCERPT')[0]
                if 'forpub2024' in galileo_name:
                    galileo_name = galileo_name.split('forpub2024')[0]
            dfmeta = load_metacsv(ssi_id)
            # we take the subsolar azimuth from the meta dict
            subsolarAZ = dfmeta[dfmeta['SSI_clock_ID'] == ssi_id]['SUB_SOLAR_AZIMUTH'].values[0]
            # also take incidence, phase and emission angle
            incidence = dfmeta[dfmeta['SSI_clock_ID'] == ssi_id]['INCIDENCE_ANGLE'].values[0]
            phase = dfmeta[dfmeta['SSI_clock_ID'] == ssi_id]['PHASE_ANGLE'].values[0]
            emission = dfmeta[dfmeta['SSI_clock_ID'] == ssi_id]['EMISSION_ANGLE'].values[0]
            # NOTE: I think I could improve the code below...
            try:
                northAZ = north_azimuth[galileo_name]
            except KeyError:
                # in this case, no northAZ was found in the dictionary
                # print('no north azimuth found for {}. We take it from the meta file.'.format(polypath.stem))
                northAZ = dfmeta[dfmeta['SSI_clock_ID'] == ssi_id]['NORTH_AZIMUTH'].values[0]
            ### NEW 2025-11-21
            # take pixel-wise subsolar sun azimuth values from phocubes
            # I have individual cubes as well as mosaics.
            # Caveats: 
            # hardcoded cub file path, 
            # Equirectangular projection only, but tessellated files are in orthographic projection
            # I retrieve pixelwise values by transforming lon/lat in meters into pixels (with coordm_to_lon, coordm_to_lat, which are for EUROPA equirectangular)
            current = os.getcwd()
            # I know this will only work for me... could be fixed in the future
            titaniach = Path(current.split('Caroline')[0]) / 'Caroline'
            phocubs = titaniach / 'isis/data/galileo/Europa_isis4/'
            # print('phocubs path is {}, galileo_name is {}'.format(phocubs, galileo_name))
            if 'galileopath' not in locals(): # this checks if the galileopath variable is assigned and therefore stored as a local variable
                # we only load in the pho_arr if not already done (remember, we are in the for loop of the tessellated tiles)
                # this saves a significant amount of computational power
                if SSI_ID_toggle != 'false':
                    # either, galileo_name is the individual id, if SSI_ID_toggle == True
                    # we search the individual cubes
                    pho_paths = sorted((phocubs / 'mapped_equirectangular_sunazimuth_fornoseam').glob('*.cub'))
                else:
                    # else, galileo_name contains the mosaic name
                    # we search the mosaics
                    pho_paths = sorted((phocubs / 'automos_sunazimuth_mosaics' / 'Europa_Mosaics_Equirectangular').glob('*.cub'))
                galileopath = 0
                for indivpath in pho_paths:
                    if galileo_name in indivpath.as_posix():
                        galileopath = indivpath
                        continue
                # print(galileopath)
                # if not pd.isnull(galileopath):
                # open the cube
                pho_dataset = gdal.Open(galileopath, gdal.GA_ReadOnly)
                # this line will simply throw an error if galileopath is still 0
                pho_arr = pho_dataset.ReadAsArray()
                pho_ulx, pho_xres, pho_xskew, pho_uly, pho_yskew, pho_yres, pho_lrx, pho_lry = get_geotransform(pho_dataset)
                # pho_meanres = abs(pho_xres)


        # print(northAZ)
        # print(subsolarAZ)

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
            arr, geoms, dataset, masks, mask_ids, sparse_masks, holes_counts = load_tiff_arr_specified(polypath, geopath)
            # print('length geoms: {}'.format(len(geoms['features'])))
        except ValueError: # then, there is no saved polygon for this geopath. This might happen if the geotiff is just on the border, or if only an excerpt was mapped (17ESREGMAP02)
            print('WARNING: No polygon data for {}. I continue without this array. still, we add area data.'.format(polypath))
            # still add area data?
            # yes! might also just be an 'empty' geojson!
            # linament AREA TO CSV
            save_empty_linea_area(polypath, savebasepath, geopath)            
            continue
        ulx, xres, xskew, uly, yskew, yres, lrx, lry = get_geotransform(dataset)
        meanres = abs(xres)

        # testing coordinates
        # read in projection file and extract lambda_0 and phi_1?
        # NOTE: lambda and phi are in radians
        with open(ortho_proj_path.joinpath(polypath.stem + '.prj')) as orth:
            ortho_proj = orth.read()
            phi_1 = np.pi/180 *float(ortho_proj.split('"Latitude_Of_Center",')[1].split(']')[0]) # latitude, e.g. '-28.3222' degrees, to radians
            lambda_0 = np.pi/180 *float(ortho_proj.split('"Longitude_Of_Center",')[1].split(']')[0]) # center longitude (east-positive) defined in projection, e.g. '288.8334' degrees, to radians
        # for testing: print(ortho_coordm_to_lon_lat(lrx, lry, lambda_0, phi_1)) # e.g. (-24151.124936296364, 27191.029295620425, 5.04109381972979, -0.494315641408338)

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

        # update on 2024-07-06: calculate cumulative area (semantic segmentation) of each category,
        # in pixel and km^2
        # define dicts for cumulative count
        area_id_masks = {1: [], 2: [], 3: [], 4: []}
        # area_cumulative_mask = np.zeros((masks.shape[0], masks.shape[1]))

        # we simply sum along an axis for all masks with the same id_int. This could be done without a for loop, but with index selection.
        # however, since I anyway loop through the masks, and since this is the less error-prone way to do it, I add the masks during the loop and sum along the mask axis after the loop

        # go through masks, every mask contains exactly one feature
        for mask_idx in range(masks.shape[-1]): # range(486,masks.shape[-1]): # masks.shape[-1]

            ######### STEP 1: calculate azimuth by fitting a line. 
            # also, we need the current mask and the bounding box
            # extract the current mask
            m = masks[:,:,mask_idx] 
            # get bounding box and cutout
            p,k,w,h, cut, extrac_cut, p_center, k_center = cut_bounding_box(m, mask_idx, extrac)
            
            # append mask to list by indexing the dict with the id_int
            area_id_masks[geoms['features'][mask_idx]['properties']['id_int']].append(m)
            # we do this step before we filter out masks that are too small for the other calculations.

            # test if mask is not almost empty. This might happen if a lineament got dissected during tessellation
            # for this, we simply calculate the area in px
            areapx = len(cut[cut>0])
            if areapx < area_limit:
                print('empty mask found. deleting and skipping. except for area calculation.')
                # do area calculation!
                # still ad lineament density of zero
                # linament AREA TO CSV
                save_empty_linea_area(polypath, savebasepath, geopath)
                # add to list of indices to delete
                del_idcs.append(mask_idx)
                continue # then, we skip 

            
            azimuth, RMS = fit_line_to_mask(sparse_masks[mask_idx], plot=False) # NOTE (2022-02): I could simply have taken the 'cut' instead of the sparse_masks... but it doesn't matter too much, since I take only the indices!!
            
            # define screen width:
            # for this, we need the approximate length in pixels:
            # length
            _, lengthpx, widthpx, _ = fit_line_to_mask_MINAREARECT(sparse_masks[mask_idx]) # we take the minimum area rectangle for the length
            if lengthpx > 400:
                screen_width = round(lengthpx/20)
            else:
                # else, to get a confident width, we need a low number of screening width
                screen_width = 2

            # run 'get_profile' only to get the approximate width
            _, lengthpx, widthpx, _ = fit_line_to_mask_MINAREARECT(sparse_masks[mask_idx]) 

            # Security measure: make sure that widthpx is not zero. # did not happen!
            if widthpx == 0 or lengthpx == 0:
                print('WARNING: widthpx is {}. Lengthpx is {}\n \n'.format(widthpx, lengthpx))
                # the below does not work for command line...
                continue
                # plt.imshow(cut)
                # plt.show()
                # bool_skip = bool(input("Do you want to skip this lineament? type anything (yes), 0 (leave empty): " ))
                # if bool_skip:
                #     continue
                # else:
                #     bool_terminate = bool(input("Do you want to stop and investigate? (otherwise, it gets segmented!) 1 (yes), 0 (leave empty): "))
                #     if bool_terminate:
                #         break

            ########## STEP 2: add new features.
            # get center lon and lat
            lon = px_to_xcoord(p_center, ulx, xres)
            lat = px_to_ycoord(k_center , uly, yres)
            # to degrees
            # with a dedicated function for Orthographic coordinate transformation
            lond, latd = ortho_coordm_to_lon_lat(lon, lat, lambda_0, phi_1)

            # get length and width in meters
            length_m = meanres*lengthpx # length in meters. meanres must come from the same projection as lengthpx
            width_m = int(widthpx) *meanres
            if math.isnan(lengthpx) == True:
                print('Problem for the length: {}'.format(lengthpx))
            # area to meter squared
            # 1 px is equal to xres*yres m^2
            aream = abs(areapx*xres*yres)
            # append all properties
            # add a column with the category name
            geoms['features'][mask_idx]['properties']['unit'] = class_dict[geoms['features'][mask_idx]['properties']['id_int']]
            # add azimuth, RMS, length, width, area, lon, lat
            geoms['features'][mask_idx]['properties']['azimuth'] = round(azimuth, 2)
            geoms['features'][mask_idx]['properties']['RMS'] = round(RMS, 2)
            geoms['features'][mask_idx]['properties']['RMS_pwpx'] = round(RMS/widthpx, 2)
            geoms['features'][mask_idx]['properties']['length_km'] = round(1e-3*length_m, 4) # to km
            geoms['features'][mask_idx]['properties']['width_km'] = round(1e-3*width_m, 4) # to km
            geoms['features'][mask_idx]['properties']['area_km2'] = round(aream*1e-6, 5) # to km2, from m2 to km2 you multiply by (10^-3)^2 = 10^-6
            geoms['features'][mask_idx]['properties']['lon_d'] = round(360-lond, 3) # to make it west_positive on the paper (but not in the projection)
            geoms['features'][mask_idx]['properties']['lat_d'] = round(latd, 3)
            # mean IFs and stds
            meanIF = np.nanmean(extrac_cut[extrac_cut>0])
            stdIF = np.nanstd(extrac_cut[extrac_cut>0])
            geoms['features'][mask_idx]['properties']['meanIF'] = float(nan_to_num(meanIF)) # '{:.3f}'.format(meanIF)
            geoms['features'][mask_idx]['properties']['stdIF'] = float(nan_to_num(stdIF)) # '{:.5f}'.format(stdIF)

            # check length again
            if geoms['features'][mask_idx]['properties']['length_km'] == 0.0:
                print('WARNGING AGAIN: I FOUND AN EMPTY LENGTH FEATURE\n \n')

            if illum_info:
                try:
                    ## sunazimuth pixelwise extraction
                    # I have to solve projection issues. polygons are loaded as orthographic.
                    # lon/lat from orthographic lineament to coordm (coordinates in meters) with new implemented function equi_coord_deg_to_coordm
                    # extract the central meridian of the current phocube: (e.g. 180)
                    pho_phi0 = float(pho_dataset.GetProjection().split('"central_meridian",')[1].split(']')[0])
                    equi_lonm, equi_latm = equi_coord_deg_to_coordm(lond, latd, pho_phi0)
                    # coordm to pixels with function coordm_to_lat and coordm_to_lon
                    pho_x_px = xcoord_to_px(np.array(equi_lonm), pho_ulx, pho_xres)
                    pho_y_px = ycoord_to_px(np.array(equi_latm), pho_uly, pho_yres)
                    # extract sunazimuth:
                    sunazi = pho_arr[pho_y_px, pho_x_px] # NOTE: x and y are changed to extract the right pixels (Tested!)
                    if sunazi < 0: 
                        # then, this is unfortunately a no-data array
                        # we can as a very quick-and-dirty hack just take the mean of the subimg.
                        #but better, we look for the nearest neighbor:
                        # https://gis.stackexchange.com/questions/207700/nearest-numpy-array-element-whose-value-is-less-than-the-current-element 
                        a,b = np.where(pho_arr > 0)
                        xli = np.stack([a,b]).T # create a 2-dim array of (row,column)
                        scidist = np.argmin(cdist(np.array([[pho_y_px, pho_x_px]]), xli)) # first distances are calculated between (row, col) of your input value, and the nearest index value is selected
                        nearest_pho_y_px, nearest_pho_x_px = xli[scidist] # this is result (row,col) of the nearest element above 0
                        # extract again, with nearest neighbor coordinates:
                        sunazi = pho_arr[nearest_pho_y_px, nearest_pho_x_px]
                    # NOTE: don't forget that the extracted sunazimuth is in a wrong coordinate system:
                    # (east is 0, south is 90deg)
                    # change it to a geological crs (North is 0, East is 90 deg)
                    sunazi = (sunazi + 90)%360 # tested
                    # assign it:
                    geoms['features'][mask_idx]['properties']['sunazimuth'] = sunazi
                except ValueError:
                    print('Sunazimuth from phocube failed. falling back to hardcoded values.')
                    # add sun azimuth and other illumination conditions from the hardcoded values (although it is the same for all features)
                    # there are two cases. solved geometrically in Book5 and tested on three examples in science_with_manual_segmentations.pptx
                    if subsolarAZ >= northAZ:
                        sunAZ = subsolarAZ - northAZ
                    else:
                        sunAZ = 360 - northAZ + subsolarAZ # (From) = 360° - (northAZ + subsolarAZ)
                    geoms['features'][mask_idx]['properties']['sunazimuth'] = sunAZ
                geoms['features'][mask_idx]['properties']['incidence'] = incidence
                geoms['features'][mask_idx]['properties']['phase'] = phase
                geoms['features'][mask_idx]['properties']['emission'] = emission

            # add resolution (this can have small discrepancies between tessells, I think due to gdal functions. But in the decimal order (e.g. 0.2))
            geoms['features'][mask_idx]['properties']['res_mpx'] = round(meanres, 2)

        # linament AREA TO CSV
        save_csv_path_dt = savebasepath / '_'.join(polypath.stem.split('_')[:-1]) / 'csv_area_files'
        os.makedirs(save_csv_path_dt, exist_ok=True)
        # calculate the lineament area per class
        area_semantic_dict = {1: {'fraction': 0, 'km2': 0}, 
                        2: {'fraction': 0, 'km2': 0}, 
                        3: {'fraction': 0, 'km2': 0}, 
                        4: {'fraction': 0, 'km2': 0}, 
                        'total_linea': {'fraction': 0, 'km2': 0}, 
                        'total_area': {'fraction': 0, 'km2': 0}, # NOTE: this is the total area as a safeguard. the fraction should in the end add up to 1 (after merging, 'Merge_polygon_calculations.py')
                       }
        # NO! this would not account for tessels that are half black, in border regions!! total_frame_area_px = masks.shape[0] * masks.shape[1]
        total_frame_area_px = len(arr[arr>0])
        total_linea_area_masks = []
        for key_id, id_mask in area_id_masks.items():
            if len(id_mask) == 0:
                continue
            id_mask_stacked = np.stack(id_mask) # shape e.g. (4, 105, 85) --> (num_masks, shape0, shape1)
            # sum along the first axis, and convert to bool (this essentially converts to semantic segmentation, but not totally, because one pixel can be cateogrised as different lineament categories)
            id_mask_summed = np.sum(id_mask_stacked, axis=0).astype(dtype=bool)*1
            # add summed mask to total area masks
            total_linea_area_masks.append(id_mask_summed)
            # calculate area in px
            id_areapx = len(id_mask_summed[id_mask_summed>0])
            # calculate fraction of lineament area over total area
            area_frac = id_areapx / total_frame_area_px
            # write to dict
            area_semantic_dict[key_id]['fraction'] = round(area_frac, 4)
            # and area in km
            # 1 px is equal to xres*yres m^2
            id_aream = abs(id_areapx*xres*yres)
            id_area_km2 = round(id_aream*1e-6, 5)
            # append to dict
            area_semantic_dict[key_id]['km2'] = id_area_km2
        # WATCH OUT HERE: We need to calculate the total lineament area from all masks stacked at the same time!!
        # because otherwise, the total_linea_fraction can be higher than 1, for example if the full image is covered by band predictions, and there are double ridge predictions ON TOP!
        # stack masks for fully-semantic total linea area
        total_linea_msum = np.sum(np.stack(total_linea_area_masks), axis=0).astype(dtype=bool)*1
        total_linea_area_px = len(total_linea_msum[total_linea_msum>0])
        total_linea_area_km2 = abs(total_linea_area_px*xres*yres*1e-6) # tested! correct
        # append total area and fraction
        total_frame_area_km2 = abs(total_frame_area_px*xres*yres*1e-6)
        area_semantic_dict['total_linea']['km2'] = total_linea_area_km2
        area_semantic_dict['total_linea']['fraction'] = round(total_linea_area_km2/total_frame_area_km2, 5)
        area_semantic_dict['total_area']['km2'] = total_frame_area_km2
        area_semantic_dict['total_area']['fraction'] = 1.0
        # NOTE: the total linea area fraction is NOT the sum of the individual fractions due to overlaps!

        # dict to pandas dataframe, to csv
        area_semantic_df = pd.DataFrame(area_semantic_dict)
        area_semantic_df.to_csv(save_csv_path_dt.joinpath(polypath.stem + '_linea_area.csv')) # Here, we need the full path, including the children _1 added!


        # now, we can savely delete all features that are stored in del_idcs list
        # ATTENTION: you cannot delete iteratively, as this will shift indices!!! We find a solution: delete biggest indices first --> this will not shift the earlier indices
        del_idcs.reverse()
        for delidx in del_idcs:
            del geoms['features'][delidx]

        # add the coordinate reference system
        geoms['crs'] ={'type': 'name',
        'properties': {'name': ortho_proj}} # we take the imported projection! It can happen that 276 degrees gets converted to 84 degrees somehow.

        # save to geojson!
        # use date for saving to prevent writing errors
        # dt_ymd = now.strftime("%Y_%m_%d")
        save_json_path_dt = savebasepath / '_'.join(polypath.stem.split('_')[:-1]) / 'json_files'
        save_crsjson_path = savebasepath / '_'.join(polypath.stem.split('_')[:-1]) / 'json_files_crs'
        geojson_filepath = save_json_path_dt.joinpath(polypath.stem + '_calc.geojson')
        # make an individual folder for each day (to prevent a huge chaos!)
        os.makedirs(save_json_path_dt, exist_ok=True)
        os.makedirs(save_crsjson_path, exist_ok=True)
        with open(geojson_filepath, 'w') as f:
            dump(geoms, f, cls=CustomEncoder) # the custom encoder makes sure that float32 gets converted to float

        # to shapefile
        shapefiles_path = savebasepath / '_'.join(polypath.stem.split('_')[:-1]) 
        os.makedirs(shapefiles_path / 'shape_files', exist_ok=True)
        geojfile = geojson_filepath
        # convert to shapefile
        command = 'ogr2ogr -nlt POLYGON -overwrite -makevalid -skipfailures {} {}'.format(shapefiles_path.as_posix() + '/shape_files/' + geojfile.stem + '.shp', shapefiles_path.as_posix() + '/json_files/' + geojfile.name)
        # print(command)
        os.system(command)


        # edit on 2024-03-01: project individual tiles back to equirectangular (or other projection).
        # because we want to merge the children feature layers later
        os.makedirs(shapefiles_path / 'back_to_prj/shape_files', exist_ok=True)
        os.makedirs(shapefiles_path / 'back_to_prj/json_files', exist_ok=True)
        os.makedirs(shapefiles_path / 'back_to_prj/json_files_crs', exist_ok=True)
        command = 'ogr2ogr -nlt MULTIPOLYGON -overwrite -skipfailures -makevalid -s_srs {} -t_srs {} {} {}'.format(ortho_prj_path.joinpath("{}.prj".format(polypath.stem)).as_posix(), equipath.joinpath(equifile).as_posix(), shapefiles_path.as_posix() + '/back_to_prj/json_files/' + geojfile.name, shapefiles_path.as_posix() + '/json_files/' + geojfile.name)
        # print(command)
        os.system(command)
        # unfortunately, the below does not add the crs at all to the geojson. no errors thrown
        # # somehow, crs is not added in 'back_to_prj'! therefore, we simply add it here
        # command = 'ogr2ogr -nlt MULTIPOLYGON -overwrite -a_srs {} {} {}'.format(equipath.joinpath(equifile).as_posix(), shapefiles_path.as_posix() + '/back_to_prj/json_files_crs/' + geojfile.name, shapefiles_path.as_posix() + '/back_to_prj/json_files/' + geojfile.name) # output input
        # print(command)
        # os.system(command)
        # # convert to shapefile NOTE: this does not work flawless yet! (thinks that source SRS is WGS) --> new: also does not work with -s_srs added.
        # command = 'ogr2ogr -nlt MULTIPOLYGON -overwrite -skipfailures -makevalid -s_srs {} -t_srs {} {} {}'.format(ortho_prj_path.joinpath("{}.prj".format(polypath.stem)).as_posix(), equipath.joinpath(equifile).as_posix(), shapefiles_path.as_posix() + '/back_to_prj/shape_files/' + geojfile.stem + '.shp', shapefiles_path.as_posix() + '/back_to_prj/json_files_crs/' + geojfile.name)
        # # print(command)
        # os.system(command)


        # print time
        end = time.time()
        print('needed {:.2f} seconds'.format(end-start_pidx))

    # measure time for all polypaths
    end = time.time()
    print('needed {:.2f} seconds'.format(end-start_all))


#%%

