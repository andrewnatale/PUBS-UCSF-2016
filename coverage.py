#!/usr/bin/env python2
import cPickle as pic
from collections import Counter
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

allele_dict = pic.load(open('../allele_dic_with_WT.pkl','rb'))
translate = pic.load(open('../translate.pkl','rb'))
amino_to_number = pic.load(open('../aminotonumber.pkl','rb'))

#print allele_dict
#print translate
#print amino_to_number

translated_dict = {}
for barcode in allele_dict:
    position, codon = allele_dict[barcode]
    if codon != 'WT':
        codon = codon.replace('T', 'U')
        aa = translate[codon]
        code_val = amino_to_number[aa]
        #print str(position) + ' ' + str(codon) + ' ' + str(aa)
        translated_dict[barcode] = (position, code_val, aa)

frequency_data = Counter(translated_dict.values())
frequency_array = np.zeros((21,76))

for data in frequency_data:
    position = data[0]
    mutation = data[1]
    frequency_array[mutation, position - 2] = frequency_data[data]

#print frequency_array

row_labels = range(2, 78, 1)

column_labels = []
for value in sorted(amino_to_number.values()):
    for key in amino_to_number:
        if value == amino_to_number[key]:
            column_labels.append(key)

sns.set(font_scale=0.7)

grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
ax = sns.heatmap(frequency_array, ax=ax,
                  cbar_ax=cbar_ax,
                  cbar_kws={"orientation": "horizontal"},
                  xticklabels=row_labels,
                  yticklabels=column_labels)
plt.sca(ax)
plt.yticks(rotation=0)
plt.xticks(rotation=90)
plt.show()
