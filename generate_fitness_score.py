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
from plot_function import heatmap_function

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

    # sum up all the barcode slopes into mutant positions
    for barcode in barcode_slopes:
        mutant_slopes[codon_to_aa[allele_map[barcode]]][0:3] \
          += barcode_slopes[barcode][0:3]
        mutant_slopes[codon_to_aa[allele_map[barcode]]][3] += 1

    # use the count to get the average
    for pos in sorted(mutant_slopes):
        if mutant_slopes[pos][3] != 0:
            mutant_slopes[pos] = mutant_slopes[pos] / mutant_slopes[pos][3]

    wt_slopes = np.zeros((1,3))
    wt_count = 0
    for pos in mutant_slopes:
        if pos[1] == ub_seq[pos[0]-1]:
            wt_slopes += mutant_slopes[pos][0:3]
            wt_count += 1
    wt_slopes = wt_slopes / float(wt_count)
    print wt_count
    print wt_slopes

    relative_fitness = {}
    for pos in mutant_slopes:
        relative_fitness[pos] = mutant_slopes[pos][0:3] - wt_slopes

    return relative_fitness

def process_fitness(fitness_data, native_seq_file):

    # mapping of aa types to numbers that loosely
    # orders them from hydrophobic to polar with 0 as 'stop'
    aminotonumber = {'STOP': 0, 'W': 1,
                     'F': 2, 'Y': 3,
                     'L': 4, 'I': 5,
                     'M': 6, 'V': 7,
                     'C': 8, 'A': 9,
                     'G': 10, 'P': 11,
                     'S': 12, 'T': 13,
                     'N': 14, 'Q': 15,
                     'H': 16, 'R': 17,
                     'K': 18, 'D': 19,
                     'E': 20}

    # load a reference sequence - this will be used to normalize the data
    # as well as to provide labels. this should be a fasta file and only
    # the first record in the file will be used
    native_seq = native_seq_file

    # open fitness data
    #data = open(fitness_data_file, 'r')
    #datalines = data.readlines()
    #data.close()
    #datalines = [i.strip('\r\n') for i in datalines]
    # remove column labels (but save it just in case)
    #fitness_header = datalines.pop(0)

    # count the number of replicates in the input file
    #expreps = len(datalines[0].split(',')) - 2
    expreps = 3
    # read data into a dictionary, replacing the aa letter code with a number
    # which will be used to arrange them (see 'aminotonumber' dictionary)
    fitness_dict = fitness_data
    #for line in datalines:
    #    elem = line.split(',')
    #    fitness_dict[(elem[0], aminotonumber[elem[1]])] = \
    #      [float(i) for i in elem[2:]]

    # get length of sequence and first resi position
    seq_length = len(Counter([i[0] for i in fitness_dict])) +1
    #for key, value in Counter([i[0] for i in fitness_dict]).iteritems():
    #    print key, value
    first_resi = int(sorted([i[0] for i in fitness_dict])[0])
    #print first_resi
    # create an empty 3d array of the correct size for all the data
    fitness_array = np.zeros((expreps,21,seq_length))

    for elem in fitness_dict:
        #print fitness_dict[elem][i]
        fitness_array[:, int(aminotonumber[elem[1]]), int(elem[0]) - first_resi] \
          = fitness_dict[elem]

    # look for and remember bad data points
    #bad_data = []
    #print 'Bad data points:'
    #for index, j in np.ndenumerate(fitness_array):
    #    if j == 999.9:
    #        bad_data.append(index)
    #        print index
    """
    # normalization types
    # None: don't do any normalization, results in
    # a mix of pos and neg values
    if normalization == None:
        for index, j in np.ndenumerate(fitness_array):
            if index in bad_data:
                fitness_array[index] = 0.0
        # average all the expreps
        mean_fitness_array = np.mean(fitness_array, axis=0)
    # 'stop': use the highest stop codon fitness as a cutoff, set it to
    # zero and shift all the data, then discard values below zero
    elif normalization == 'stop':
        for i in range(expreps):
            # get max stop codon values for this replicate
            max_stop = np.amax(fitness_array[i,aminotonumber['*'],:])
            # shift everything by the max stop codon value
            fitness_array[i,:,:] = fitness_array[i,:,:] - max_stop
        # set bad data points and values below zero to zero
        for index, j in np.ndenumerate(fitness_array):
            if j < 0:
                fitness_array[index] = 0.0
            if index in bad_data:
                fitness_array[index] = 0.0
        # average all the expreps
        mean_fitness_array = np.mean(fitness_array, axis=0)
    # 'prob': same as stop, but then normalize each position column
    # so that it sums to 1
    elif normalization == 'prob':
        for i in range(expreps):
            # get max stop codon values for this replicate
            max_stop = np.amax(fitness_array[i,aminotonumber['*'],:])
            # shift everything by the max stop codon value
            fitness_array[i,:,:] = fitness_array[i,:,:] - max_stop
        # set any values below zero to zero
        for index, j in np.ndenumerate(fitness_array):
            if j < 0:
                fitness_array[index] = 0.0
            if index in bad_data:
                fitness_array[index] = 0.0
        # average all the expreps
        mean_fitness_array = np.mean(fitness_array, axis=0)
        # normalize each position independently
        for i in range(seq_length):
            mean_fitness_array[:,i] = mean_fitness_array[:,i] \
              / np.sum(mean_fitness_array[:,i])
    # end of the normalization loop
    """

    """
    # compose data labels
    sequence_labels = []
    for n in range(first_resi, first_resi + seq_length, 1):
        sequence_labels.append(('%s%d') % (native_seq[0].seq[n-1], n))

    mutation_labels = []
    for value in sorted(aminotonumber.values()):
        for key in aminotonumber:
            if value == aminotonumber[key]:
                mutation_labels.append(key)

    # return a 2d array which is the mean of the replicates
    # stacked in the 3d array, plus some useful labels for plots
    return (mean_fitness_array, sequence_labels, mutation_labels)
    """
    return np.mean(fitness_array[0:2,:,:], axis=0)

fitness = get_fitness(barcode_slopes)

#for i in range(0,3):
#    for elem in fitness:
#        print fitness[elem][i,]

heatmap_function(process_fitness(fitness, ub_seq))
plt.show()
