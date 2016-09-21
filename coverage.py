#!/usr/bin/env python2
import cPickle as pic
from collections import Counter
import numpy as np

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
#print frequency_data
#position = None
frequency_array = np.zeros((76,21))

for data in frequency_data:
    position = data[0]
    mutation = data[1]
    frequency_array[position - 2, mutation] = frequency_data[data]

print frequency_array
