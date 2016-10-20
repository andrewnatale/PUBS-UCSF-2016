
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def heatmap_function(data_in):
	list_numb = range(1,79)
	list_string = []
	for i in list_numb:
		list_string.append(str(i))
	X = list_string
	Y = ['STOP', 'W', 'F', 'Y', 'L', 'I', 'M', 'V', 'C', 'A', 'G', 'P', 'S', 'T', 'N', 'Q', 'H', 'R', 'K', 'D', 'E']
	data = data_in
	df = pd.DataFrame(data, index=Y, columns=X)
	sns.heatmap(df)
