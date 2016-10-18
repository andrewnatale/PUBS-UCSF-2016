from Bio import SeqIO
from Bio.Seq import Seq
import pickle as pic
import numpy as np
import sys

indices = {'TTTCAC':0, 'GGCCAC':1, 'CGAAAC':2, 'CGTACG':3, 'CCACTC':4, 'GCTACC':5, 'ATTCCG':6, 'AGCTAG':7, 'GTATAG':8}
map_barcodes = pic.load(open("allele_dic_with_WT.pkl", "rb"))
#print map_barcodes

with open(sys.argv[1], "rU") as data_file:
    data = SeqIO.parse(data_file, "fastq")

    barcode_counts = dict(zip(map_barcodes.keys(), np.zeros((len(map_barcodes), 9))))
    #print barcode_counts
    #print len(barcode)
    #print len(barcode_counts)
    num_loops = 0
    num_excepts = 0
    num_try = 0
    for record in data:
        num_loops += 1
        #print record.description[-6:]
        #print record.seq[:18]
        if record.description[-6:] in indices.keys():
            barcode = record.seq[:18]
            #rc_barcode = barcode.reverse_complement()
            #print barcode, rc_barcode
            try:
                barcode_counts[str(barcode.reverse_complement())][indices[record.description[-6:]]] += 1
                num_try += 1
            except KeyError:
                num_excepts += 1
            pass

print num_loops, num_excepts, num_try
print np.sum(np.array(barcode_counts.values()), axis=0)
pic.dump(barcode_counts, open("data_no_qc.pkl", 'wb'))
