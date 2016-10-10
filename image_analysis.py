#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from skimage import data, io, filters, img_as_float, exposure, transform
from skimage.morphology import disk, watershed
from skimage.segmentation import slic, join_segmentations
from scipy import fftpack
from scipy import ndimage as ndi
from sys import argv

script, filename = argv
image = io.imread(filename)
image_flt = img_as_float(image)
image_sobel = filters.sobel(image)

markers = np.zeros_like(image)
markers[image < 3500] = 1
markers[image > 5300] = 2

#print image_sobel.min()
#print image_sobel.max()

seg = watershed(image_sobel, markers)
#print image.median()
plt.imshow(seg, cmap=plt.cm.gray)
#plt.imshow(image, cmap=plt.cm.gray)
plt.show()
