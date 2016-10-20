from scipy import stats
import numpy as np
import math
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

# 76 aa long ubqitin sequence
ub_seq = 'MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG*'

# delete any barcodes that have no counts
start_count = len(data.keys())
for barcode in data.keys():
    if 0.0 in data[barcode]:
        del data[barcode]
end_count = len(data.keys())

#print start_count, end_count

# convert counts to frequency of each barcode in each sample
for n in range(0,9):
    sample_sum = np.sum(np.array([data[i][n] for i in data]))
    for barcode in data:
        data[barcode][n] = data[barcode][n] / sample_sum

# log2 of the frequency
for barcode in data:
    data[barcode] = np.log2(data[barcode])

# x vals (# of generations) for slope calculations
time_points = [0.0, 1.74, 3.71, 0.0, 1.75, 3.37, 0., 2., 4.]

# calculate slope for every barcode for each experiment
barcode_slopes = dict(zip(data.keys(), np.zeros((len(data.keys()), 3))))
for barcode in data:
    for n in range(0,3):
        yvals = data[barcode][(n*3):(n*3+3)]
        xvals = time_points[(n*3):(n*3+3)]
        line = stats.linregress(xvals, yvals)
        barcode_slopes[barcode][n] = line[0]

# build a dictionary to help convert barcodes to mutants
codon_list = sorted(Counter(allele_map.values()).keys())
codon_to_aa = {}
for codon in codon_list:
    if codon[1] == 'WT':
        if ub_seq[codon[0]-1] == '*':
            codon_to_aa[codon] = (codon[0], 'STOP')
        else:
            codon_to_aa[codon] = (codon[0], ub_seq[codon[0]-1])
    else:
        codon_to_aa[codon] = (codon[0], translation[codon[1].replace('T', 'U')])

mutant_list = sorted(Counter(codon_to_aa.values()).keys())
mutant_slopes = dict(zip(mutant_list, np.zeros((len(mutant_list), 4))))

# get the average wt slope (only exact wt codons, not wt synonym)
wt_slopes = np.zeros((1,3))
wt_count = 0
for barcode in barcode_slopes:
    if allele_map[barcode][1] == 'WT':
        wt_slopes += barcode_slopes[barcode][0:3]
        wt_count += 1
wt_slopes = wt_slopes / float(wt_count)

#print wt_count
#print wt_slopes

# sum up all the barcode slopes into mutant positions
for barcode in barcode_slopes:
    mutant_slopes[codon_to_aa[allele_map[barcode]]][0:3] \
      += barcode_slopes[barcode][0:3]
    mutant_slopes[codon_to_aa[allele_map[barcode]]][3] += 1

# use the count to get the average
for pos in sorted(mutant_slopes):
    if mutant_slopes[pos][3] != 0:
        mutant_slopes[pos] = mutant_slopes[pos] / mutant_slopes[pos][3]

# subtract avg wt from every point
fitness_scores = {}
for pos in mutant_slopes:
    fitness_scores[pos] = mutant_slopes[pos][0:3] - wt_slopes

pkl.dump(fitness_scores, open("fitness_scores.pkl", 'wb'))
