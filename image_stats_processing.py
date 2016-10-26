#!/usr/bin/env python2
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from sys import argv
import nd2reader
from skimage import color, filters, measure, exposure
from scipy import stats
import os
import pickle

data_file_list = ['WellD11.pic', 'WellE07.pic']

with open(data_file_list[0], 'rb') as data_file:
    data = pickle.load(data_file)

access = {
'A':'DAPI_GFP_distance',
'B':'CY3_GFP_distance',
'C':'CY3_DAPI_distance',
'D':'DAPI_GFP_pixel_coor',
'E':'CY3_GFP_pixel_coor',
'F':'CY3_DAPI_pixel_coor',
'G':'segmentation_labels'}

plt.hist(data[access['A']], bins=50)
plt.show()
plt.hist(data[access['B']], bins=50)
plt.show()
plt.hist(data[access['C']], bins=50)
plt.show()

plt.hist([z[0] for z in data[access['D']]], bins=50)
plt.show()
plt.hist([z[1] for z in data[access['D']]], bins=50)
plt.show()
plt.hist([z[0] for z in data[access['E']]], bins=50)
plt.show()
plt.hist([z[1] for z in data[access['E']]], bins=50)
plt.show()
plt.hist([z[0] for z in data[access['F']]], bins=50)
plt.show()
plt.hist([z[1] for z in data[access['F']]], bins=50)
plt.show()

for img in data[access['G']]:
    plt.imshow(img)
    plt.show()
