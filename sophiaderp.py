import numpy as np
import skimage
from skimage import filters, io as io
from skimage.morphology import disk, watershed
from skimage.filters.rank import median
from sys import argv 
import matplotlib.pyplot as plt, matplotlib.cm as cm
from scipy import ndimage as ndi

script, file_name = argv

photo_file = (io.imread(file_name, as_grey=True))
photo_matrix = skimage.img_as_float(photo_file)

#median = median(photo_matrix, selem = disk(2))

def segment_image(photo_matrix, background_threshold, foreground_threshold):
    edges = skimage.filters.sobel(photo_matrix)
    markers = np.zeros_like(photo_matrix)
    foreground, background = 1 , 2
    markers[photo_matrix < background_threshold] = background
    markers[photo_matrix > foreground_threshold] = foreground
    ws = watershed(edges, markers)
    segmentation_matrix = ndi.label(ws == foreground)[0]
    return segmentation_matrix

segment_image(photo_matrix, 0.2,0.7 )

#np.ma.masked_where(photo_matrix > 0.2,photo_matrix)


plt.imshow(photo_file, cmap = plt.cm.gray)
plt.show()
