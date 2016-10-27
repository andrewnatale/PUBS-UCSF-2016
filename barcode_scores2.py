#!/usr/bin/env python2
from scipy import stats
import numpy as np
import math
import pickle as pkl
import sys
from collections import Counter
from matplotlib import pyplot as plt
from copy import deepcopy
from statsmodels import robust

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

# possibly dangerous filtering and pseudocounting
start_count = len(data_set.keys())
strange_bcodes = []
for barcode in data_set.keys():
    if 0 in data_set[barcode][[0,3,6]]:
        print data_set[barcode]
        del data_set[barcode]
    else:
        data_set[barcode] += 0.1
end_count = len(data_set.keys())
print start_count, end_count
 # or np.sum(data_set[barcode]) <= 30:
# convert counts to frequency
for n in range(0,9):
    sample_sum = np.sum(np.array([data_set[i][n] for i in data_set]))
    for barcode in data_set:
        data_set[barcode][n] = data_set[barcode][n] / sample_sum

# log2 of the frequency
for barcode in data_set:
     data_set[barcode] = np.log2(data_set[barcode])

# calculate slopes for every barcode for each experiment (also s.e.)
barcode_slopes = dict(zip(data_set.keys(), np.zeros((len(data_set.keys()), 3))))
barcode_se = dict(zip(data_set.keys(), np.zeros((len(data_set.keys()), 3))))

for barcode in data_set:
    for n in range(0,3):
        yvals = data_set[barcode][(n*3):(n*3+3)]
        xvals = time_points[(n*3):(n*3+3)]
        line = stats.linregress(xvals, yvals)
        # store the slope and s.e.
        barcode_slopes[barcode][n] = line[0]
        barcode_se[barcode][n] = line[4]

# get the average wt slope for each experiment (only exact wt codons, not wt synonym)
wt_slopes = np.zeros((3))
wt_count = 0
for barcode in barcode_slopes:
    if allele_map[barcode][1] == 'WT':
        wt_slopes += barcode_slopes[barcode][0:3]
        wt_count += 1
wt_slopes = wt_slopes / float(wt_count)

# build a empty dictionaries which will hold lists of all the individual barcode fitness scores for each position
mutant_list = sorted(Counter(codon_to_aa.values()).keys())
component_slopes_day1 = {}
component_slopes_day2 = {}
component_slopes_ctrl = {}
for mutant in mutant_list:
    component_slopes_day1[mutant] = []
    component_slopes_day2[mutant] = []
    component_slopes_ctrl[mutant] = []

# put the scores into the dictionary and subtract the avg wt value for that experiment
for barcode in data_set:
    component_slopes_day1[codon_to_aa[allele_map[barcode]]].append(barcode_slopes[barcode][0] - wt_slopes[0])
    component_slopes_day2[codon_to_aa[allele_map[barcode]]].append(barcode_slopes[barcode][1] - wt_slopes[1])
    component_slopes_ctrl[codon_to_aa[allele_map[barcode]]].append(barcode_slopes[barcode][2] - wt_slopes[2])

#pkl.dump(component_slopes_day1, open("component_slopes_day1.pkl", 'wb'))
#pkl.dump(component_slopes_day2, open("component_slopes_day2.pkl", 'wb'))
#pkl.dump(component_slopes_ctrl, open("component_slopes_ctrl.pkl", 'wb'))

def position_stats(score_list):
    arr_of_list = np.array(score_list)
    list_len = len(score_list)
    list_mean = np.mean(arr_of_list)
    list_median = np.median(arr_of_list)
    # standard deviation
    list_std = np.std(arr_of_list)
    # standard error
    list_se = list_std / math.sqrt(list_len)
    # median absolute deviation
    list_mad = robust.scale.mad(arr_of_list)
    # trim out values outside the range median +/- mad
    trimmed_list = [i for i in arr_of_list if i >= list_median - list_mad and i <= list_median + list_mad]
    trimmed_len = len(trimmed_list)
    trimmed_mean = np.mean(np.array(trimmed_list))
    return trimmed_mean

fitness_scores = dict(zip(mutant_list, np.zeros((len(mutant_list), 3))))
for pos in fitness_scores:
    fitness_scores[pos][0] = position_stats(component_slopes_day1[pos])
    fitness_scores[pos][1] = position_stats(component_slopes_day2[pos])
    fitness_scores[pos][2] = position_stats(component_slopes_ctrl[pos])

pkl.dump(fitness_scores, open("fitness_scores_test.pkl", 'wb'))
#pkl.dump(component_slopes_ctrl, open("component_slopes_ctrl.pkl", 'wb'))
