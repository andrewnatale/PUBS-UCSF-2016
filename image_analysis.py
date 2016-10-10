#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from skimage import data, io, filters, img_as_float, exposure, transform
from skimage.morphology import disk, watershed, opening, dilation, erosion
from skimage.segmentation import slic, join_segmentations
from skimage.feature import canny
from scipy import fftpack
from scipy import ndimage as ndi
from sys import argv

script, filename = argv
image = io.imread(filename)
image_flt = img_as_float(image)

#median_filtered = filters.median(image, disk(3))
#edges = canny(median_filtered/1., sigma=6)

#image_sobel = filters.sobel(image)

markers = np.zeros_like(image)
#markers[image < 3600] = 2
#markers[image > 3600] = 1
markers[image > 5200] = 1
fill = opening(markers)
fill = opening(fill)
fill = opening(fill)
fill = opening(fill)
fill = opening(fill)
fill = opening(fill)
#print image_sobel.min()
#print image_sobel.max()

#seg = watershed(image_sobel, markers)
#print image.median()
plt.imshow(fill, cmap=plt.cm.gray)
#plt.imshow(image, cmap=plt.cm.gray)
plt.show()
