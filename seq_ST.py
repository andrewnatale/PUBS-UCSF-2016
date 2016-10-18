#from Bio import *
#from Bio import SeqIO
#import scipy
from scipy import stats
import numpy as np
#import numpy
import math
import seaborn as sns
import matplotlib.pyplot as plt
import pickle as pkl
import sys
from collections import Counter

# open up pickle file with data
with open(sys.argv[1], 'rb') as data_file:
    data = pkl.load(data_file)

# open pickle file containing map of barcodes to alleles
with open('allele_dic_with_WT.pkl', 'r') as allele_map_file:
    allele_map = pkl.load(allele_map_file)

# open pickle file with codon to amino acid translations
with open('translate.pkl', 'r') as translation_file:
    translation = pkl.load(translation_file)

# get a list containing each (allele, pos#) tuple once
#allele_list = sorted(Counter(allele_map.values()).keys())
# use the list to create 2 new empty dictionaries
#allele_counts = dict(zip(allele_list, np.zeros((len(allele_list), 9))))

# sum up all barcodes that point to the same allele
#for barcode in allele_map:
#    allele_counts[allele_map[barcode]] += data[barcode]

# normalize to total reads per sample
start_count = len(data.keys())
for barcode in data.keys():
    if 0.0 in data[barcode]:
        del data[barcode]
end_count = len(data.keys())

#print start_count, end_count

for n in range(0,9):
    sample_sum = np.sum(np.array([data[i][n] for i in data]))
    for barcode in data:
        data[barcode][n] = data[barcode][n] / sample_sum

for barcode in data:
    data[barcode] = np.log2(data[barcode])

# x vals (# of generations) for slope calculations
time_points = [0.0, 1.74, 3.71, 0, 1.75, 3.37, 0., 2., 4.]

# calculate slope for every allele for each experiment
barcode_slopes = dict(zip(data.keys(), np.zeros((len(data.keys()), 3))))

for barcode in data:
    for n in range(0,3):
        yvals = data[barcode][(n*3):(n*3+3)]
        xvals = time_points[(n*3):(n*3+3)]
        line = stats.linregress(xvals, yvals)
        barcode_slopes[barcode][n] = line[0]

#print barcode_slopes

def get_fitness(barcode_slopes):
    # get avg wt slope for each experiment
    wt_slopes = np.zeros((1,3))
    wt_count = 0
    for barcode in barcode_slopes:
        if allele_map[barcode][1] == 'WT':
            wt_slopes += barcode_slopes[barcode]
            wt_count += 1

    wt_slopes = wt_slopes / float(wt_count)

    relative_fitness = {}
    for barcode in barcode_slopes:
        relative_fitness[barcode] = barcode_slopes[barcode] / wt_slopes

# compare fitness

#for allele in relative_fitness:
#    print allele, relative_fitness[allele]
"""
print 'sum of counts per sample:'
print np.sum(np.array(data.values()), axis=0)
print 'after normalization:'
print np.sum(np.array(allele_counts.values()), axis=0)
print 'max and min slopes per sample:'
print np.amax(np.array(allele_slopes.values()), axis=0)
print np.amin(np.array(allele_slopes.values()), axis=0)
print 'max and min after normalization to avg wt slope'
print np.amax(np.array(relative_fitness.values()), axis=0)
print np.amin(np.array(relative_fitness.values()), axis=0)


for key, value in dict.iteritems():
    #day 1
    a = value[0]
    b = value[1]
    c = value[2]
    #day 2
    d = value[3]
    e = value[4]
    f = value[5]

    score1 = [a, b, c]
    score2 = [d, e, f]
    lists = np.array([score1, score2])
    # y value is list of averages
    avg = np.mean(lists, axis=0)
    x = [0, 1, 2]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, avg)
    rsquared = r_value**2
#
#handle = open('fastq_test', 'rU')
#records = list(SeqIO.parse(handle, 'fastq'))
#
##each record has 6 lines, the sequence is in the 6th line
#for record in records:
#        print repr(record.seq)



# barcodes = [day1, day1, day1, day2, day2, day2, wt, wt, wt]
#wt_T0=barcodes[barcode[6]]
#wt_T1=barcodes[barcode[7]]
#wt_T2=barcodes[barcode[8]]

#wt_slope = 0.00001

#def fitness_score(barcode):
#    score1 = [a, b, c]
#    score2 = [d, e, f]
#    lists = np.array([score1, score2])
#    # y value is list of averages
#    avg = np.mean(lists, axis=0)
#    x = [0, 1, 2]
#    slope, intercept, r_value, p_value, std_err = stats.linregress(x, avg)
#    rsquared = r_value**2
#
#

"""
