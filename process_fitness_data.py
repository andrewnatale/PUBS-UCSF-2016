import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import pickle as pkl
import sys
from collections import Counter
#from plot_function import heatmap_function

with open('fitness_scores.pkl', 'rb') as fitness_scores_file:
    fitness_scores = pkl.load(fitness_scores_file)

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

# 76 aa long ubqitin sequence
ub_seq = 'MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG*'

# we don't want anything from position 0 or 1
for pos in fitness_scores.keys():
    if pos[0] == 0 or pos[0] == 1 or pos[0] == 77:
        del fitness_scores[pos]

# how many data points per position
expreps = len(fitness_scores.values()[0][0])
# get length of sequence and first resi position
seq_length = len(Counter([i[0] for i in fitness_scores]))
first_resi = int(sorted([i[0] for i in fitness_scores])[0])

# create an empty 3d array of the correct size for all the data
fitness_array = np.zeros((expreps+2,21,seq_length))

for elem in fitness_scores:
    coord = (int(aminotonumber[elem[1]]), int(elem[0]) - first_resi)
    fitness_array[2:, coord[0], coord[1]] = fitness_scores[elem]
    # this slice will be a mask that tells us which positions have no data
    fitness_array[0,  coord[0], coord[1]] = 1
    # this slice indicates wt or not
    if elem[1] == ub_seq[elem[0]-1]:
        fitness_array[1,  coord[0], coord[1]] = 1

# compose data labels
sequence_labels = []
for n in range(first_resi, first_resi + seq_length, 1):
    sequence_labels.append(('%s%d') % (ub_seq[n-1], n))
#print sequence_labels

mutation_labels = []
for value in sorted(aminotonumber.values()):
    for key in aminotonumber:
        if value == aminotonumber[key]:
            mutation_labels.append(key)
#print mutation_labels

def plot_hmap(data, row_labels, column_labels):
    sns.set(font_scale=1.2)
    grid_kws = {"height_ratios": (.9, .05), "hspace": .2}
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    ax = sns.heatmap(data, ax=ax,
                     linewidths=0.5,
                     #cmap='Blues',
                     cbar_ax=cbar_ax,
                     cbar_kws={"orientation": "horizontal"},
                     xticklabels=row_labels,
                     yticklabels=column_labels)
    plt.sca(ax)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.show()

# takes a 2D array
def normalize_to_mean_stop(fitness_array):
    # get max stop codon values for this replicate
    #max_stop = np.amax(fitness_array[aminotonumber['STOP'],:])
    avg_stop = np.mean(fitness_array[aminotonumber['STOP'],:])
    # shift everything by the max stop codon value
    fitness_array = fitness_array - avg_stop
    # neg values go to zero
    fitness_array[fitness_array < 0] = 0
    return fitness_array

# takes a 2D array
def normalize_to_max_stop(fitness_array):
    # get max stop codon values for this replicate
    max_stop = np.amax(fitness_array[aminotonumber['STOP'],:])
    #avg_stop = np.mean(fitness_array[aminotonumber['STOP'],:])
    # shift everything by the max stop codon value
    fitness_array = fitness_array - max_stop
    # neg values go to zero
    fitness_array[fitness_array < 0] = 0
    return fitness_array

def normalize_info_content(fitness_array):
    # normalize each position independently
    for i in range(seq_length):
        fitness_array[:,i] = fitness_array[:,i] \
          / np.sum(fitness_array[:,i])
    return fitness_array

# average of two experiments
#plot_hmap(np.mean(fitness_array[2:4,:,:], axis=0),sequence_labels,mutation_labels)

# fitness difference between experimental and control conditions
#plot_hmap(np.mean(fitness_array[2:4,:,:], axis=0) - fitness_array[4,:,:],sequence_labels,mutation_labels)

# control
#plot_hmap(fitness_array[3,:,:],sequence_labels,mutation_labels)

# mask which shows where we have no data
#plot_hmap(fitness_array[0,:,:],sequence_labels,mutation_labels)

# mask which shows which positions are wt
#plot_hmap(fitness_array[1,:,:],sequence_labels,mutation_labels)

# try normalization
norm = normalize_info_content(normalize_to_max_stop(np.mean(fitness_array[2:4,:,:], axis=0)))
plot_hmap(norm,sequence_labels,mutation_labels)


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
