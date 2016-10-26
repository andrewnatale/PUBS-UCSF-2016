#!/usr/bin/env python2
import numpy as np
import math
#import matplotlib.pyplot as plt
from scipy import ndimage as ndi
#from sys import argv
import nd2reader
from skimage import color, filters, measure, exposure, io
from scipy import stats
import os
import pickle

# normalize
def normalize(image):
    #norm_image = image.astype(np.double)
    norm_image = image - np.min(image)
    norm_image /= np.max(norm_image)
    return norm_image

# copy/convert an nd2 image object to a float64 array
def nd2tofloat(nd2):
    new_img = np.zeros(nd2.shape)
    for i,v in np.ndenumerate(nd2):
        new_img[i] = float(v)
    return normalize(new_img)

# segment a bright field image
def segment_bf(image):
    # make a copy
    segment_target = np.copy(image)
    # should already be normalized at load
    #segment_target = normalize(segment_target)

    # subtract background
    segment_target = segment_target - filters.gaussian(segment_target, sigma=15)
    segment_target = segment_target + abs(np.min(segment_target))
    # enhance contrast
    plow, phigh = np.percentile(segment_target, (8, 50))
    segment_target = exposure.rescale_intensity(segment_target, in_range=(plow, phigh))
    # smooth
    segment_target = filters.gaussian(segment_target, sigma=2)
    # otsu threshold
    thresh = filters.threshold_otsu(segment_target)

    # use threshold to create a binary mask
    binary = np.zeros_like(segment_target)
    binary[segment_target < thresh] = 1
    # fill holes
    binary = ndi.morphology.binary_fill_holes(binary)

    # make initial labels
    label_1 = measure.label(binary)
    # filter for size and shape of regions
    filter_mask = np.ones_like(np.bincount(label_1.ravel()), dtype=bool)
    filter_mask[0] = 0
    ratio_thresh_mult = 1.4
    for region in measure.regionprops(label_1):
        if region.perimeter > ratio_thresh_mult \
          * math.sqrt(4*math.pi*region.area):
            filter_mask[region.label] = 0
        if region.area < 100 or region.area > 2000:
            filter_mask[region.label] = 0
    cleaned_binary = filter_mask[label_1]
    label_2 = measure.label(cleaned_binary)

    return image, segment_target, label_2

def pixel_correlation(labels, reference, gfp):
    # use labels to set up list of properties
    properties_reference = measure.regionprops(labels, reference)
    properties_gfp = measure.regionprops(labels, gfp)

    # iterate over all labeled regions (cells) to get image wide coorelation
    slopes_r_vals = []
    for i in range(len(measure.regionprops(labels))):
        # construct a mask identifying which regions have 0 intesity in both channels
        intensity_reference = normalize(properties_reference[i].intensity_image)
        intensity_gfp = normalize(properties_gfp[i].intensity_image)
        mask = np.ones_like(intensity_reference, dtype=bool)
        # cut out the bottom 10% of pixels
        ref_thres = np.percentile(intensity_reference[np.nonzero(intensity_reference)], 10)
        gfp_thres = np.percentile(intensity_gfp[np.nonzero(intensity_gfp)], 10)
        for index, value in np.ndenumerate(intensity_reference):
            if value <= ref_thres and intensity_gfp[index] <= gfp_thres:
                mask[index] = False
        # pair up the pixel intesity in each channel if at least one is nonzero
        gfp_vals = []
        reference_vals = []
        for index, value in np.ndenumerate(mask):
            if value:
                gfp_vals.append(intensity_gfp[index])
                reference_vals.append(intensity_reference[index])
        # linear regression on intesity values, then save slope and r value
        line = stats.linregress(gfp_vals,reference_vals)
        slopes_r_vals.append((line[0],line[2]))
    return slopes_r_vals

def bright_feature_distance(labels, reference, gfp):
    properties_reference = measure.regionprops(labels, reference)
    properties_gfp = measure.regionprops(labels, gfp)
    distance_list = []
    for i in range(len(measure.regionprops(labels))):
        intensity_reference = normalize(properties_reference[i].intensity_image)
        intensity_gfp = normalize(properties_gfp[i].intensity_image)
        # use the brightest 10% of pixels for distance calculations
        ref_thres = np.percentile(intensity_reference[np.nonzero(intensity_reference)], 90)
        gfp_thres = np.percentile(intensity_gfp[np.nonzero(intensity_gfp)], 90)
        mask_reference = intensity_reference <= ref_thres
        mask_gfp = intensity_gfp > gfp_thres
        # plt.imshow(intensity_reference)
        # plt.show()
        # plt.imshow(mask_reference)
        # plt.show()
        # plt.imshow(intensity_gfp)
        # plt.show()
        # plt.imshow(mask_gfp)
        # plt.show()
        d_to_ref = ndi.morphology.distance_transform_edt(mask_reference)
        d = d_to_ref[mask_gfp]
        distance_list.append(d)
    return np.concatenate(distance_list).ravel()

# def examine_segmentation(image1, image2, image3):
#     fig = plt.figure()
#     a=fig.add_subplot(1,3,1)
#     imgplot = plt.imshow(image1, cmap=plt.cm.gray)
#     a.set_title('Original')
#     a=fig.add_subplot(1,3,2)
#     imgplot = plt.imshow(image2, cmap=plt.cm.gray)
#     a.set_title('Modified')
#     a=fig.add_subplot(1,3,3)
#     imgplot = plt.imshow(image3)
#     a.set_title('Segmentation')
#     plt.show()

# def color_channels_with_labels(labels, gfp_in, dapi_in, cy3_in):
#     gfp_colors = color.label2rgb(labels, gfp_in, bg_label=0)
#     dapi_colors = color.label2rgb(labels, dapi_in, bg_label=0)
#     cy3_colors = color.label2rgb(labels, cy3_in, bg_label=0)
#     fig = plt.figure()
#     a=fig.add_subplot(1,3,1)
#     imgplot = plt.imshow(gfp_colors)
#     a.set_title('gfp')
#     a=fig.add_subplot(1,3,2)
#     imgplot = plt.imshow(dapi_colors)
#     a.set_title('dapi')
#     a=fig.add_subplot(1,3,3)
#     imgplot = plt.imshow(cy3_colors)
#     a.set_title('cy3')
#     plt.show()

nd2_files = []
for filename in os.listdir(os.getcwd()):
    if filename.endswith(".nd2"):
        nd2_files.append(filename)

for filename in nd2_files:
    well_id = filename.split('_')[1]
    compound_img = nd2reader.Nd2(filename)
    dapi_to_gfp_distance = []
    cy3_to_gfp_distance = []
    cy3_to_dapi_distance = []
    dapi_to_gfp_pixel_coor = []
    cy3_to_gfp_pixel_coor = []
    cy3_to_dapi_pixel_coor = []
    well_data = {}
    well_segmentation = []
    for i in range(0,6,1):
        for img in compound_img.select(fields_of_view=i):
            if img.channel == 'FITClong':
                gfp = nd2tofloat(img)
            elif img.channel == 'DAPI':
                dapi = nd2tofloat(img)
            elif img.channel == 'CY3':
                cy3 = nd2tofloat(img)
            elif img.channel == 'BF-Cy3':
                bf = nd2tofloat(img)
        h, i, j = segment_bf(bf)
        dapi_to_gfp_distance.extend(bright_feature_distance(j, dapi, gfp))
        cy3_to_gfp_distance.extend(bright_feature_distance(j, cy3, gfp))
        cy3_to_dapi_distance.extend(bright_feature_distance(j, dapi, cy3))
        dapi_to_gfp_pixel_coor.extend(pixel_correlation(j, dapi, gfp))
        cy3_to_gfp_pixel_coor.extend(pixel_correlation(j, cy3, gfp))
        cy3_to_dapi_pixel_coor.extend(pixel_correlation(j, dapi, cy3))
        well_segmentation.append(j)
    well_data['DAPI_GFP_distance'] = dapi_to_gfp_distance
    well_data['CY3_GFP_distance'] = cy3_to_gfp_distance
    well_data['CY3_DAPI_distance'] = cy3_to_dapi_distance
    well_data['DAPI_GFP_pixel_coor'] = dapi_to_gfp_pixel_coor
    well_data['CY3_GFP_pixel_coor'] = cy3_to_gfp_pixel_coor
    well_data['CY3_DAPI_pixel_coor'] = cy3_to_dapi_pixel_coor
    well_data['segmentation_labels'] = well_segmentation
    pickle.dump(well_data, open('%s.pic' % well_id, 'wb'))
