#!/usr/bin/env python2
from scipy import stats
import numpy as np
import math
import pickle as pkl
import sys
from collections import Counter
from matplotlib import pyplot as plt
from copy import deepcopy

# data outputs:
# - array of per aa/pos raw fitness scores for every experiment + masks
# - counts, frequency, r-squared for every codon

# open up pickle file with data dictionary in the form:
# 'barcode sequence' : [array of counts per barcode in each sample]
with open(sys.argv[1], 'rb') as data_file:
    data_set = pkl.load(data_file)

# data and dictionaries needed to process sequencing counts

# open pickle file containing map of barcodes to alleles
with open('allele_dic_with_WT.pkl', 'r') as allele_map_file:
    allele_map = pkl.load(allele_map_file)
# open pickle file with codon to amino acid translations
with open('translate.pkl', 'r') as translation_file:
    translation = pkl.load(translation_file)
# 76 aa long ubqitin sequence
ub_seq = 'MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG*'
# x vals (# of generations) for slope calculations from our experiments
time_points = [0.0, 1.74, 3.71, 0.0, 1.75, 3.37, 0., 2., 4.]

# build a dictionary to help convert barcodes to mutants
codon_to_aa = {}
for codon in sorted(Counter(allele_map.values()).keys()):
    if codon[1] == 'WT':
        if ub_seq[codon[0]-1] == '*':
            codon_to_aa[codon] = (codon[0], 'STOP')
        else:
            codon_to_aa[codon] = (codon[0], ub_seq[codon[0]-1])
    else:
        codon_to_aa[codon] = (codon[0], translation[codon[1].replace('T', 'U')])

# aggregate counts to the codon level
#codon_list = sorted(Counter(allele_map.values()).keys())
#codon_counts = dict(zip(codon_list, np.zeros((len(codon_list), 9))))

# calculate slope for every codon for each experiment
#codon_slopes = dict(zip(codon_counts.keys(), np.zeros((len(codon_counts.keys()), 6))))

barcode_slopes = dict(zip(data_set.keys(), np.zeros((len(data_set.keys()), 6))))

start_count = len(data_set.keys())
for barcode in data_set.keys():
    if 0.0 in data_set[barcode]:
        del data_set[barcode]
end_count = len(data_set.keys())

for n in range(0,9):
    sample_sum = np.sum(np.array([data_set[i][n] for i in data_set]))
    for barcode in data_set:
        data_set[barcode][n] = data_set[barcode][n] / sample_sum

# log2 of the frequency
for barcode in data_set:
     data_set[barcode] = np.log2(data_set[barcode])

for barcode in data_set:
    for n in range(0,3):
        yvals = data_set[barcode][(n*3):(n*3+3)]
        xvals = time_points[(n*3):(n*3+3)]
        line = stats.linregress(xvals, yvals)
        # store the slope and r-squared
        barcode_slopes[barcode][n] = line[0]
        barcode_slopes[barcode][n+3] = line[4]
        #codon_stats[codon][n*2+1] = line[2]**2

mutant_list = sorted(Counter(codon_to_aa.values()).keys())
component_slopes = {}
component_slopes_ctrl = {}
for mutant in mutant_list:
    component_slopes[mutant] = []
    component_slopes_ctrl[mutant] = []

for barcode in data_set:
    component_slopes[codon_to_aa[allele_map[barcode]]].append(float(barcode_slopes[barcode][0]))
    component_slopes[codon_to_aa[allele_map[barcode]]].append(float(barcode_slopes[barcode][1]))
    component_slopes_ctrl[codon_to_aa[allele_map[barcode]]].append(float(barcode_slopes[barcode][2]))

#for mutant in component_slopes:
#    print mutant, component_slopes[mutant]

pkl.dump(component_slopes, open("component_slopes.pkl", 'wb'))
pkl.dump(component_slopes_ctrl, open("component_slopes_ctrl.pkl", 'wb'))
"""
for barcode in data_set:
    codon_counts[allele_map[barcode]] += data_set[barcode]

# delete any codons that have no counts at any position
start_count = len(codon_counts.keys())
for codon in codon_counts.keys():
    if 0.0 in codon_counts[codon]:
        del codon_counts[codon]
end_count = len(codon_counts.keys())

# codon stats will hold (pos, codon):[]
#codon_stats = dict(zip(codon_counts.keys(), np.zeros((len(codon_counts.keys()), 6))))

# add log2 of raw codon counts to the stats dictionary
#for codon in codon_counts:
#        codon_stats[codon][0] = math.log(np.sum(codon_counts[codon][0:3]), 2)
#        codon_stats[codon][2] = math.log(np.sum(codon_counts[codon][3:6]), 2)
#        codon_stats[codon][4] = math.log(np.sum(codon_counts[codon][6:]), 2)

# convert counts to frequency of each codon in each sample
for n in range(0,9):
    sample_sum = np.sum(np.array([codon_counts[i][n] for i in codon_counts]))
    for codon in codon_counts:
        codon_counts[codon][n] = codon_counts[codon][n] / sample_sum

# log2 of the frequency
for codon in codon_counts:
     codon_counts[codon] = np.log2(codon_counts[codon])

# calculate slope for every codon for each experiment
codon_slopes = dict(zip(codon_counts.keys(), np.zeros((len(codon_counts.keys()), 6))))
for codon in codon_counts:
    for n in range(0,3):
        yvals = codon_counts[codon][(n*3):(n*3+3)]
        xvals = time_points[(n*3):(n*3+3)]
        line = stats.linregress(xvals, yvals)
        # store the slope and r-squared
        codon_slopes[codon][n] = line[0]
        codon_slopes[codon][n+3] = line[2]**2
        #codon_stats[codon][n*2+1] = line[2]**2

# get the average wt slope (only exact wt codons, not wt synonym)
wt_slopes = np.zeros((1,3))
wt_count = 0
for codon in codon_slopes:
    if codon[1] == 'WT':
        wt_slopes += codon_slopes[codon][0:3]
        wt_count += 1
wt_slopes = wt_slopes / float(wt_count)

# sum up all the codon slopes into mutant positions
mutant_list = sorted(Counter(codon_to_aa.values()).keys())
mutant_slopes = dict(zip(mutant_list, np.zeros((len(mutant_list), 4))))
for codon in codon_slopes:
    mutant_slopes[codon_to_aa[codon]][0:3] \
      += codon_slopes[codon][0:3]
    mutant_slopes[codon_to_aa[codon]][3] += 1

# use the count to get the average
for pos in sorted(mutant_slopes):
    if mutant_slopes[pos][3] != 0:
        mutant_slopes[pos] = mutant_slopes[pos] / mutant_slopes[pos][3]

# subtract avg wt from every point
fitness_scores = {}
for pos in mutant_slopes:
    fitness_scores[pos] = mutant_slopes[pos][0:3] - wt_slopes

#for codon in codon_stats:
#    print codon, codon_stats[codon]

def compare_hist(flat_array1, flat_array2):
    fig = plt.figure()
    a=fig.add_subplot(1,2,1)
    hist = plt.hist(flat_array1)
    a.set_title('barcode')
    a=fig.add_subplot(1,2,2)
    hist = plt.hist(flat_array2)
    a.set_title('codon')
    plt.show()

def single_hist(flat_array):
    plt.hist(flat_array, bins=500, range=(0.0, 1000.))
    plt.show()

#print np.sum(np.array(codon_counts.values()), axis=0)
#print np.sum(np.array(codon_counts.values()))
#print start_count, end_count

pkl.dump(fitness_scores, open("fitness_scores1.pkl", 'wb'))

#log_freqs = []
#r_squareds = []
#for codon in codon_stats:
#    log_freqs.append(codon_stats[codon][0])
#    log_freqs.append(codon_stats[codon][2])
#    log_freqs.append(codon_stats[codon][4])
#    r_squareds.append(codon_stats[codon][1])
#    r_squareds.append(codon_stats[codon][3])
#    r_squareds.append(codon_stats[codon][5])
"""
