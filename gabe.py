#!/usr/bin/env python2
from skimage.morphology import watershed, disk
from skimage.filters import sobel
from skimage import exposure
from skimage import transform as tf
from skimage.segmentation import slic, join_segmentations
from skimage import data, img_as_float
from scipy import fftpack
from skimage.filters.rank import median
from scipy import ndimage as ndi
from skimage import data
# ----------------------------------------------------------------------
import skimage.feature
import skimage.filters
import os
import sys
import skimage
#import nd2reader
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.misc
import numpy as np
#GIT_DIR = os.path.normpath(os.path.dirname(os.getcwd()))
#PUBS_DIR = os.path.normpath(os.path.dirname(GIT_DIR))
#IMAGE_DIR = os.path.normpath(PUBS_DIR + '/images_day_1/')
# ----------------------------------------------------------------------

def segment_image(photo_matrix, background_threshold, foreground_threshold):
    edges = sobel(photo_matrix)
    markers = np.zeros_like(photo_matrix)
    foreground, background = 1, 2
    markers[photo_matrix < background_threshold] = background
    markers[photo_matrix > foreground_threshold] = foreground
    ws = watershed(edges, markers)
    segmentation_matrix = ndi.label(ws == foreground)[0]
    return segmentation_matrix


def count_DAPI(img):
	marker_mat = np.zeros_like(img)

	# block_size = 7
	# binary_adaptive = skimage.filters.threshold_adaptive(img, block_size, offset=10)
	global_thresh = skimage.filters.threshold_otsu(img)
	binary_global = img > global_thresh



	# edges = sobel(binary_global)
	edges = sobel(img)

	background, foreground = 1, 2
	background_threshold = 30
	foreground_threshold = 150
	markers = np.zeros_like(binary_global)
	# print binary_global
	markers[img < global_thresh] = background
	markers[img > global_thresh] = foreground
	# print np.amax(markers)
	markers = binary_global + 1

	segmentation = watershed(edges, markers)
	segmentation = ndi.binary_fill_holes(segmentation - 1)
	labeled_coins, _ = ndi.label(segmentation)

	image_label_overlay = skimage.color.label2rgb(labeled_coins, image=img)



	# print labeled_coins
	# # print markers
	# # print global_thresh
	plt.imshow(edges, cmap = plt.cm.gray)
	# # plt.imshow(image_label_overlay, cmap = plt.cm.gray)
	plt.imshow(img, cmap=plt.cm.gray, interpolation='nearest')
	plt.contour(segmentation, [0.5], linewidths=1.2, colors='y')
	plt.show()

im = scipy.misc.imread(sys.argv[1])
# print im
count_DAPI(im)
