import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import pickle as pkl
import sys
from collections import Counter
import pandas as pd
from scipy import stats
from statsmodels import robust
#from plot_function import heatmap_function

with open('fitness_scores_v4.pkl', 'rb') as fitness_scores_file:
    fitness_scores = pkl.load(fitness_scores_file)

# load up ubiqitin monomer rosetta ddg data
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

# opposite dictionary for reverse lookups
numbertoamino = {}
for key, value in aminotonumber.iteritems():
    numbertoamino[value] = key

# 76 aa long ubqitin sequence
ub_seq = 'MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG*'

# we don't want anything from position 0, 1, or 77
for pos in fitness_scores.keys():
    if pos[0] == 0 or pos[0] == 1 or pos[0] == 77:
        del fitness_scores[pos]

# how many data points per position
expreps = len(fitness_scores.values()[0])
# get length of sequence and first resi position
seq_length = len(Counter([i[0] for i in fitness_scores]))
first_resi = int(sorted([i[0] for i in fitness_scores])[0])

# create an empty 3d array of the correct size for all the data
fitness_array = np.zeros((expreps+2,21,seq_length))

for elem in fitness_scores:
    coord = (int(aminotonumber[elem[1]]), int(elem[0]) - first_resi)
    fitness_array[0:3, coord[0], coord[1]] = fitness_scores[elem]
    # this slice will be a mask that tells us which positions have no data
    #fitness_array[3, coord[0], coord[1]] = 1
    # this slice indicates wt or not
    if elem[1] == ub_seq[elem[0]-1]:
        fitness_array[3,  coord[0], coord[1]] = 1
    # this slice contains some rosetta ddg data
    try:
        ddg_score = ddg_dict[elem]
    # if it can't be looked up, means its wt with score of 0
    except KeyError:
        ddg_score = 0.0
    fitness_array[4, coord[0], coord[1]] = ddg_score

# dictionary to help index the array
data_sets = {'day1':0,
             'day2':1,
             'ctrl':2,
             'wt?':3,
             'ddg_m':4}

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

def plot_hmap(data, row_labels, column_labels, min_=None, max_=None, cmap_=None):
    sns.set(font_scale=1.2)
    grid_kws = {"height_ratios": (.9, .05), "hspace": .2}
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    ax = sns.heatmap(data, ax=ax,
                     linewidths=0.2,
                     cmap=cmap_,
                     cbar_ax=cbar_ax,
                     cbar_kws={"orientation": "horizontal"},
                     xticklabels=row_labels,
                     yticklabels=column_labels,
                     vmin=min_,
                     vmax=max_)
    plt.sca(ax)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.show()

def threshold(fitness_array, threshold):
    fitness_array = fitness_array - threshold
    fitness_array[fitness_array < 0] = 0
    return fitness_array

def rescale(fitness_array, scale_factor=None):
    if not scale_factor:
        # get avg wt
        wt_sum = 0.0
        wt_count = 0.0
        for index, val in np.ndenumerate(fitness_array):
            if ub_seq[index[1] + first_resi] == numbertoamino[index[0]]:
                wt_sum += fitness_array[index]
                wt_count += 1
        scale_factor = wt_sum / wt_count
    return fitness_array / scale_factor

def normalize_info_content(fitness_array):
    # normalize each position independently
    for i in range(seq_length):
        fitness_array[:,i] = fitness_array[:,i] \
          / np.sum(fitness_array[:,i])
    return fitness_array

def compare_hist(flat_array1, flat_array2):
    fig = plt.figure()
    a=fig.add_subplot(1,2,1)
    hist = plt.hist(flat_array1, bins = 15)
    a.set_title('p-fluorophenylalanine')
    a=fig.add_subplot(1,2,2)
    hist = plt.hist(flat_array2, bins = 15)
    a.set_title('control')
    plt.show()

def pssm_out(data, sequence_labels, mutation_labels, identifier, out_file):
    # open a file for writing
    with open(out_file, 'w') as pssm_file:
        # put an identifier on the first line
        pssm_file.write('ID %s\n' % identifier)
        # write position values
        pssm_file.write('P0 ')
        for elem in mutation_labels:
            pssm_file.write('%s ' % elem)
        pssm_file.write('\n')
        # write data for each position on a line
        # start by writing the position labels
        position_count = 0
        for i in sequence_labels:
            pssm_file.write('%s ' % i[1:])
            for j in data[:,position_count]:
                pssm_file.write('%s ' % str(j*20))
            pssm_file.write('\n')
            position_count += 1
        # write footer
        pssm_file.write('XX\n')
        pssm_file.write('//\n')



# break out avg experimental fitness and control fitness
#avg_exp = np.mean(fitness_array[[data_sets['day1'],data_sets['day2']],:,:], axis=0)
day1_exp = fitness_array[data_sets['day1'],:,:]
day2_exp = fitness_array[data_sets['day2'],:,:]
ctrl = fitness_array[data_sets['ctrl'],:,:]

# histogram of stops vs non_stops
day1_scores = day1_exp[1:,:].flatten()[~np.isnan(day1_exp[1:,:].flatten())]
day2_scores = day2_exp[1:,:].flatten()[~np.isnan(day2_exp[1:,:].flatten())]
ctrl_scores = ctrl[1:,:].flatten()[~np.isnan(ctrl[1:,:].flatten())]
day1_stops = day1_exp[0,:].flatten()[~np.isnan(day1_exp[0,:].flatten())]
day2_stops = day2_exp[0,:].flatten()[~np.isnan(day2_exp[0,:].flatten())]
ctrl_stops = ctrl[0,:].flatten()[~np.isnan(ctrl[0,:].flatten())]

# scatter plot, day1 vs day2
line = stats.linregress(day1_scores, day2_scores)
slope, intercept, rvalue, pvalue, stderr = line
print line
fig, ax = plt.subplots()
plt.scatter(day1_scores, day2_scores)
xlim = ax.get_xlim()
x = np.linspace(xlim[0], xlim[1], 100)
best_fit_line = float(slope) * x + intercept
plt.plot(x, best_fit_line)
plt.show()

# histogram of stops and non-stops
non_stops = np.concatenate([day1_scores, day2_scores, ctrl_scores])
stops = np.concatenate([day1_stops, day2_stops, ctrl_stops])
plt.hist(stops, bins=50)
plt.hist(non_stops, bins=50, alpha=0.5)
plt.show()

# get the median of all stop codon fitness values
median_stop = np.median(stops)
stops_mad = robust.scale.mad(stops)

thresh = float(median_stop) + float(stops_mad)

# threshold exp and control data using median control stop codon
day1_thresh = threshold(day1_exp, thresh)
day2_thresh = threshold(day2_exp, thresh)
ctrl_thresh = threshold(ctrl, thresh)
# rescale each so wt=1
day1_rescale = rescale(day1_thresh, -1.0*thresh)
day2_rescale = rescale(day2_thresh, -1.0*thresh)
ctrl_rescale = rescale(ctrl_thresh, -1.0*thresh)

average = (day1_rescale + day2_rescale) / 2
difference_map = average - ctrl_rescale



plot_hmap(day1_thresh,sequence_labels,mutation_labels)
plot_hmap(day2_thresh,sequence_labels,mutation_labels)
plot_hmap(average,sequence_labels,mutation_labels)
plot_hmap(ctrl_thresh,sequence_labels,mutation_labels)
plot_hmap(difference_map,sequence_labels,mutation_labels)


# empty array for z scores
# z_scores = np.zeros((3,21,seq_length), dtype=np.double)
# z_scores.fill(np.nan)
# for index, val in np.ndenumerate(day1_exp):
#     z_scores[0,index[0],index[1]] = (val - np.mean(day1_scores)) / np.std(day1_scores)
# for index, val in np.ndenumerate(day2_exp):
#     z_scores[1,index[0],index[1]] = (val - np.mean(day2_scores)) / np.std(day2_scores)
# for index, val in np.ndenumerate(ctrl):
#     z_scores[2,index[0],index[1]] = (val - np.mean(ctrl_scores)) / np.std(ctrl_scores)
#
# scatter_list = []
# for index, val in np.ndenumerate(z_scores):
#     if not math.isnan(val):
#         scatter_list.append((float(index[2]+2), val, index[1], index[0]))
# print scatter_list
#
#
#
# use_colors = {0: "red", 1: "green", 2: "blue"}
# #print [use_colors[i[3]] for i in scatter_list]
# plt.scatter([i[0] for i in scatter_list], [i[1] for i in scatter_list], color=[use_colors[i[3]] for i in scatter_list])
# plt.show()



# ******************
# raw fitness scores
# ******************

# average of perturbation experiments
#avg_perturb = np.mean(fitness_array[[data_sets['day1'],data_sets['day2']],:,:], axis=0)

#plot_hmap(avg_perturb, sequence_labels, mutation_labels)

# control experiment
#plot_hmap(fitness_array[data_sets['ctrl'],:,:], sequence_labels, mutation_labels)



# **************************
# thresholded fitness scores
# **************************

# histograms of stop codon fitness values
#perturb_stops = fitness_array[[data_sets['day1'],data_sets['day2']], aminotonumber['STOP'], :].flatten()
#ctrl_stops = fitness_array[data_sets['ctrl'], aminotonumber['STOP'], :].flatten()
#print np.median(perturb_stops)
#print np.median(ctrl_stops)
#compare_hist(perturb_stops, ctrl_stops)

# average perturb data sets thresholded based on median stop codon value
#avg_perturb = np.mean(fitness_array[[data_sets['day1'],data_sets['day2']],:,:], axis=0)
#median_stop_perturb = np.median(fitness_array[[data_sets['day1'],data_sets['day2']], aminotonumber['STOP'], :].flatten())
#rescale_perturb = rescale(arbitrary_threshold(avg_perturb, median_stop_perturb), -1.0*median_stop_perturb)
#print median_stop
#plot_hmap(rescale(arbitrary_threshold(avg_perturb, median_stop), -1.0*median_stop)[1:,:], sequence_labels, mutation_labels[1:], max_=2.5, cmap_='Greens')

# control data set thresholded based on median stop codon value
#ctrl = fitness_array[data_sets['ctrl'],:,:]
#median_stop_ctrl = np.median(fitness_array[data_sets['ctrl'], aminotonumber['STOP'], :].flatten())
#rescale_ctrl = rescale(arbitrary_threshold(ctrl, median_stop_ctrl), -1.0*median_stop_ctrl)
#print median_stop
#plot_hmap(rescale(arbitrary_threshold(ctrl, median_stop), -1.0*median_stop)[1:,:], sequence_labels, mutation_labels[1:], max_=2.5, cmap_='Greens')

# difference (perturbed - ctrl)
#plot_hmap((rescale_perturb-rescale_ctrl)[1:,:], sequence_labels, mutation_labels[1:])


# **************************
# ddg scores for ubi monomer
# **************************

# heatmap of ddg monomer data
#ddg_array = fitness_array[data_sets['ddg'],1:,:]
#print np.amin(ddg_array), np.amax(ddg_array)
#plt.hist(ddg_array.flatten(), bins=20)
#plt.show()
#plot_hmap(ddg_array, sequence_labels, mutation_labels[1:], -20.0, 20.0, 'PuOr')



# ***************************************
# generate sequence logo pssm input files
# ***************************************

# average perturb data sets thresholded based on median stop codon value
#avg_perturb = np.mean(fitness_array[[data_sets['day1'],data_sets['day2']],:,:], axis=0)
#median_stop = np.median(fitness_array[[data_sets['day1'],data_sets['day2']], aminotonumber['STOP'], :].flatten())
#thresholded = arbitrary_threshold(avg_perturb, median_stop)
#normalized = normalize_info_content(thresholded)
#plot_hmap(normalized[1:,:], sequence_labels, mutation_labels[1:])
#pssm_out(normalized[1:,:], sequence_labels, mutation_labels[1:], 'avg_perturb_norm', 'avg_perturb_norm.pssm.txt')

# control data set thresholded based on median stop codon value
#ctrl = fitness_array[data_sets['ctrl'],:,:]
#median_stop = np.median(fitness_array[data_sets['ctrl'], aminotonumber['STOP'], :].flatten())
#thresholded = arbitrary_threshold(ctrl, median_stop)
#normalized = normalize_info_content(thresholded)
#plot_hmap(normalized[1:,:], sequence_labels, mutation_labels[1:])
#pssm_out(normalized[1:,:], sequence_labels, mutation_labels[1:], 'ctrl', 'ctrl_norm.pssm.txt')



# mask which shows where we have no data
#plot_hmap(fitness_array[data_sets['data?'],:,:],sequence_labels,mutation_labels)

# mask which shows which positions are wt
#plot_hmap(fitness_array[1,:,:],sequence_labels,mutation_labels)
