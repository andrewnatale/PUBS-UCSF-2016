from Bio import SeqIO
from Bio.Seq import Seq
import pickle as pic
import numpy as np
handle = open("fastq_sample", "rU")
parse = SeqIO.parse(handle, "fastq")
indices = {'TTTCAC':0, 'GGCCAC':1, 'CGAAAC':2, 'CGTACG':3, 'CCACTC':4, 'GCTACC':5, 'ATTCCG':6, 'AGCTAG':7, 'GTATAG':8}
barcode = pic.load(open("allele_dic_with_WT.pkl", "rb"))
#print barcode
barcode_counts = dict(zip(barcode.keys(), np.zeros((len(barcode), 9))))
#print barcode_counts
#print len(barcode)
#print len(barcode_counts)
for record in parse:
	#print record.description[-6:]
	#print record.seq[:18]
	if record.description[-6:] in indices.keys():
		rc_barcode = Seq(str(record.seq[:18]))
		try:
			barcode_counts[rc_barcode.reverse_complement()][indices[record.description[-6:]]] += 1
		except KeyError:
			pass
#print np.sum(np.array(barcode_counts.values()), axis=0)
#quality is letter_annotation
handle.close()
pic.dump(barcode_counts, open("barcode_counts_retry.pkl", 'wb'))
