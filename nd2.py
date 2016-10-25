#!/usr/bin/env python2
import nd2reader
from sys import argv

script, nd2 = argv

nd2 = nd2reader.Nd2(str(nd2))

dapi = []
fitc = []
fitclong = []
cy3 = []
bfcy3 = []

for image in nd2.select(channels = 'DAPI'):
    dapi.append(image)

for image in nd2.select(channels = 'FITC'):
    fitc.append(image)

for image in nd2.select(channels = 'FITClong'):
    fitclong.append(image)

for image in nd2.select(channels = 'CY3'):
    cy3.append(image)

for image in nd2.select(channels = 'BF-Cy3'):
    bfcy3.append(image)

test = cy3[0]
print type(test)
