#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from skimage import data, io, filters, img_as_float, exposure, transform, measure
from skimage.morphology import disk, watershed, opening, dilation, erosion, closing, white_tophat, black_tophat
#from skimage.segmentation import slic, join_segmentations
#from skimage.feature import canny
from scipy import fftpack
from scipy import ndimage as ndi
from sys import argv

# imput and various iterations of image
script, filename = argv
image = io.imread(filename)
image_flt = img_as_float(image)
image_invert = 65535 - image

def segment_image(image):
    #z = closing(image, selem=np.ones((300,300)))
    gaussian_filtered = filters.gaussian(image, sigma=3)
    thresh = filters.threshold_otsu(gaussian_filtered)
    binary = np.zeros_like(gaussian_filtered)
    binary[gaussian_filtered < thresh] = 1
    filled = ndi.morphology.binary_fill_holes(binary)
    label_objects = measure.label(filled)
    #sizes = np.bincount(label_objects.ravel())
    #mask_sizes = sizes > 100
    #mask_sizes[0] = 0
    #cleaned = mask_sizes[label_objects]
    for region in measure.regionprops(label_objects):
        print region.area
        print region.perimeter
    return (image, binary, label_objects)


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

h, i, j = segment_image(image_flt)
image_comparison(i, j)
