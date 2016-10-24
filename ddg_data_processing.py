import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
#import pickle as pkl
import sys
from collections import Counter
import pandas as pd
import os

data_path = '/Users/anatale/Documents/school/UCSF/courses/PUBS/rosetta_ubq_predictions'
data_folders = ['complexes', 'diubiquitin', 'ubiquitin_monomer']

df_list = []

for dirname in data_folders:
    folder = os.path.join(data_path, dirname)
    for filename in os.listdir(folder):
        if filename.endswith('.tsv'):
            df_list.append(pd.read_table(os.path.join(folder, filename)))

print df_list
