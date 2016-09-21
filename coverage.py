#!/usr/bin/env python2
import cPickle as pic
from collections import Counter
import numpy as np
from hmap_funct import plot_hmap

allele_dict = pic.load(open('../allele_dic_with_WT.pkl','rb'))
translate = pic.load(open('../translate.pkl','rb'))
amino_to_number = pic.load(open('../aminotonumber.pkl','rb'))

translated_dict = {}
for barcode in allele_dict:
    position, codon = allele_dict[barcode]
    if codon != 'WT':
        codon = codon.replace('T', 'U')
        aa = translate[codon]
        code_val = amino_to_number[aa]
        translated_dict[barcode] = (position, code_val, aa)

frequency_data = Counter(translated_dict.values())

frequency_array = np.zeros((21,76))

for data in frequency_data:
    position = data[0]
    mutation = data[1]
    frequency_array[mutation, position - 2] = frequency_data[data]

row_labels = range(2, 78, 1)

column_labels = []
for value in sorted(amino_to_number.values()):
    for key in amino_to_number:
        if value == amino_to_number[key]:
            column_labels.append(key)

plot_hmap(frequency_array, row_labels, column_labels)
