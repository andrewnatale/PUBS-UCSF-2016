#!/usr/bin/env python2
import cPickle as pic

allele_dict = pic.load(open('../allele_dic_with_WT.pkl','rb'))
translate = pic.load(open('../translate.pkl','rb'))
amino_to_number = pic.load(open('../aminotonumber.pkl','rb'))

#print allele_dict
#print translate
#print amino_to_number

translated_dict = {}
for barcode in allele_dict:
    print barcode
    position, codon = allele_dict[barcode]
    if codon != 'WT':
        codon = codon.replace('T', 'U')
        aa = translate[codon]
        print str(position) + ' ' + str(codon) + ' ' + str(aa)
        translated_dict[barcode] = (position, aa)

print translated_dict
