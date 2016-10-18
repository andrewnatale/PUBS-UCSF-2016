from Bio import *
from Bio import SeqIO
import scipy
from scipy import stats
import numpy as np
import numpy
import seaborn
import matplotlib.pyplot as plt
import pickle as pkl

with open('barcode_counts_retry.pkl', 'rb') as data_file:
 data = pkl.load(data_file)

with open('allele_dic_with_WT.pkl', 'r') as allele_map_file:
 allele_map = pkl.load(allele_map_file)

wt_barcodes = []
for key, value in allele_map.iteritems():
    if value[1] == 'WT':
        wt_barcodes.append(key)

wt_counts = {}
for barcode in wt_barcodes:
    wt_counts[barcode] = data[barcode]

print np.sum(np.array(wt_counts.values()), axis=0)
'''
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
   '''
