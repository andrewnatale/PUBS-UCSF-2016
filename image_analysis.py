#!/usr/bin/env python2

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from skimage import data, io, filters, img_as_float, exposure, transform, measure, feature
from skimage.morphology import disk, watershed, opening, dilation, erosion, closing, white_tophat, black_tophat
#from skimage.segmentation import slic, join_segmentations
#from skimage.feature import canny
from scipy import fftpack
from scipy import ndimage as ndi
from sys import argv

# imput and various iterations of image
script, filename = argv
image = io.imread(filename)

# segment a bright field image
def segment_bf(image):
    # convert to float
    image = img_as_float(image)
    original = image

    # enhance contrast
    p2, p98 = np.percentile(image, (2, 98))
    image = exposure.rescale_intensity(image, in_range=(p2, p98))

    # subtract a broad gaussian from the background
    image = image - filters.gaussian(image, sigma=20)
    image = image + abs(np.min(image))

    # smooth and set a threshold
    smoothed = filters.gaussian(image, sigma=3)
    thresh = filters.threshold_otsu(smoothed)

    # use threshold to create a binary mask
    binary = np.zeros_like(image)
    binary[smoothed < thresh*1.1] = 1
    # fill holes
    binary = ndi.morphology.binary_fill_holes(binary)


    label_objects = measure.label(binary)
    print 
    """
    ratio_thresh_mult = 1.5
    for region in measure.regionprops(label_objects):
        if region.perimeter > ratio_thresh_mult \
          * math.sqrt(4*math.pi*region.area):
            print region.label
    """
    #plt.scatter(area, perimeter)
    #plt.show()
    #plt.hist(ratio)
    #plt.show()
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = np.zeros_like(sizes, dtype=bool)
    for index, elem in enumerate(sizes):
        if elem > 100 and elem < 5000:
            mask_sizes[index] = 1
    mask_sizes[0] = 0
    cleaned = mask_sizes[label_objects]

    #print type(mask_sizes)
    #print mask_sizes.shape
    #print mask_sizes

    #print mask_sizes.shape
    #print label_objects.ravel()
    """

    """
    return (original, cleaned)

def normalize(img):
	high = np.amax(img)
	low = np.amin(img)
	return (img - low) / (high - low)


def find_cells(img):
	strong_blur = filters.gaussian(img, 20)
	no_back = img - strong_blur
	no_back = normalize(no_back)
	equalized_no_back = exposure.equalize_hist(no_back)
	equalized_no_back = normalize(equalized_no_back)
	edges_nb = feature.canny(equalized_no_back, sigma=5)
	close_nb = ndi.binary_closing(edges_nb, structure=np.ones((3, 3)), iterations=1)
	fill_close_nb = ndi.binary_fill_holes(close_nb)
	open_fcnb = ndi.binary_opening(fill_close_nb, structure=np.ones((10, 10)))
	open_bigger = ndi.morphology.binary_dilation(open_fcnb, iterations=5)
	#border = morph.binary_dilation(open_fcnb) - open_fcnb
	#dist = ndi.distance_transform_edt(open_fcnb)
	#local_peaks = feature.peak_local_max(dist, min_distance=12, threshold_abs=4,
	#										labels=open_fcnb, indices=False)
	#markers = ndi.label(local_peaks)[0]
	#labels = morph.watershed(-dist, markers, mask=open_bigger)
	#find_boundaries = seg.find_boundaries(labels)
	return equalized_no_back, edges_nb


#mask = np.zeros_like(image)
# mask[image < 4000] = 1
#mask[image < 3500] = 1
# watershed(image_sobel, mask)
#closed = ndi.morphology.binary_closing(wall, iterations=4)
#filled = ndi.morphology.binary_fill_holes(closed)
#cells = np.zeros_like(image)
#cells[filled == 1] = 1
#cells[wall == 1] = 2
def image_comparison(image1, image2):
    fig = plt.figure()
    a=fig.add_subplot(1,2,1)
    imgplot = plt.imshow(image1, cmap=plt.cm.gray)
    a.set_title('Before')
    a=fig.add_subplot(1,2,2)
    imgplot = plt.imshow(image2, cmap=plt.cm.gray)
    a.set_title('After1')
    #a=fig.add_subplot(1,3,3)
    #imgplot = plt.imshow(image3)
    #a.set_title('After2')
    plt.show()

# image should be floating point
def image_to_hist(image):
    plt.hist(image.ravel(), bins=256, range=(0.0, 1.0), fc='k', ec='k')
    plt.show()

h, i = segment_bf(image)
#h, i = find_cells(image)
#image_comparison(h, i)
#image_to_hist(h)
