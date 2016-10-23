import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import pickle as pkl
import sys
from collections import Counter
import pandas as pd
#from plot_function import heatmap_function

with open('fitness_scores1.pkl', 'rb') as fitness_scores_file:
    fitness_scores = pkl.load(fitness_scores_file)

ddg_data = pd.read_table('uby_1ubq.tsv')
ddg_dict = {}
for index, row in ddg_data.iterrows():
    mutant = row['Mutation']
    ddg_dict[(int(''.join(i for i in mutant[2:] if i.isdigit())), mutant[-1:])] = \
      row['ddg_stability']

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
fitness_array = np.zeros((expreps+3,21,seq_length))

for elem in fitness_scores:
    coord = (int(aminotonumber[elem[1]]), int(elem[0]) - first_resi)
    fitness_array[0:3, coord[0], coord[1]] = fitness_scores[elem]
    # this slice will be a mask that tells us which positions have no data
    fitness_array[3, coord[0], coord[1]] = 1
    # this slice indicates wt or not
    if elem[1] == ub_seq[elem[0]-1]:
        fitness_array[4,  coord[0], coord[1]] = 1
    # this slice contains some rosetta ddg data
    try:
        ddg_score = ddg_dict[elem]
    # if it can't be looked up, means its wt with score of 0
    except KeyError:
        ddg_score = 0.0
    fitness_array[5, coord[0], coord[1]] = ddg_score

# dictionary to help index the array
data_sets = {'day1':0,
             'day2':1,
             'ctrl':2,
             'data?':3,
             'wt?':4,
             'ddg':5}

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

def plot_hmap(data, row_labels, column_labels, _min=None, _max=None):
    sns.set(font_scale=1.2)
    grid_kws = {"height_ratios": (.9, .05), "hspace": .2}
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    ax = sns.heatmap(data, ax=ax,
                     linewidths=0.5,
                     #cmap='Blues',
                     cbar_ax=cbar_ax,
                     cbar_kws={"orientation": "horizontal"},
                     xticklabels=row_labels,
                     yticklabels=column_labels,
                     vmin=_min,
                     vmax=_max)
    plt.sca(ax)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.show()

# takes a 2D array
def threshold_to_mean_stop(fitness_array):
    # get max stop codon values for this replicate
    #max_stop = np.amax(fitness_array[aminotonumber['STOP'],:])
    avg_stop = np.mean(fitness_array[aminotonumber['STOP'],:])
    # shift everything by the max stop codon value
    fitness_array = fitness_array - avg_stop
    # neg values go to zero
    fitness_array[fitness_array < 0] = 0
    return fitness_array

# takes a 2D array
def threshold_to_max_stop(fitness_array):
    # get max stop codon values for this replicate
    max_stop = np.amax(fitness_array[aminotonumber['STOP'],:])
    #avg_stop = np.mean(fitness_array[aminotonumber['STOP'],:])
    # shift everything by the max stop codon value
    fitness_array = fitness_array - max_stop
    # neg values go to zero
    fitness_array[fitness_array < 0] = 0
    return fitness_array

def arbitrary_threshold(fitness_array, threshold):
    fitness_array = fitness_array - threshold
    fitness_array[fitness_array < 0] = 0
    return fitness_array

def normalize_info_content(fitness_array):
    # normalize each position independently
    for i in range(seq_length):
        fitness_array[:,i] = fitness_array[:,i] \
          / np.sum(fitness_array[:,i])
    return fitness_array

# histogram of stop codon fitness scores, day1, day2, ctrl
#plt.hist(fitness_array[[data_sets['day1'],data_sets['day2'],data_sets['ctrl']], aminotonumber['STOP'], :].flatten(), bins=20)
#plt.show()

# heatmap of ddg monomer data
ddg_array = fitness_array[data_sets['ddg'],1:,:]
print np.amin(ddg_array), np.amax(ddg_array)
plt.hist(ddg_array.flatten(), bins=20)
plt.show()
#plot_hmap(ddg_array, sequence_labels, mutation_labels[1:], -10.0, 10.0)

# average of two experiments
#avg_perturb = np.mean(fitness_array[2:4,:,:], axis=0)
#plot_hmap(arbitrary_threshold(avg_perturb, -0.2), sequence_labels, mutation_labels)

# fitness difference between experimental and control conditions
#plot_hmap(np.mean(fitness_array[2:4,:,:], axis=0) - fitness_array[4,:,:],sequence_labels,mutation_labels)

# control
#plot_hmap(fitness_array[3,:,:],sequence_labels,mutation_labels)

# mask which shows where we have no data
#plot_hmap(fitness_array[0,:,:],sequence_labels,mutation_labels)

# mask which shows which positions are wt
#plot_hmap(fitness_array[1,:,:],sequence_labels,mutation_labels)

# try normalization
#norm = normalize_info_content(normalize_to_mean_stop(np.mean(fitness_array[2:4,:,:], axis=0)))
#plot_hmap(norm,sequence_labels,mutation_labels)

#for index, val in np.ndenumerate(fitness_array[1,:,:]):
#    if val == 1.:
#        if fitness_array[0,index[0],index[1]] == 1.:
#            print fitness_array[4,index[0],index[1]]
