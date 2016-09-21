import matplotlib.pyplot as plt
import numpy as np
import seaborn

column_labels = [1, 2, 3]
row_labels = [6, 5, 4]

data = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
    ]
fig, ax = plt.subplots()
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

#put major ticks at middle of cell
###########ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
###########ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)
#for more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

ax.set_xticklabels(row_labels, minor = False)
ax.set_yticklabels(column_labels, minor = False)

## set up 2D grid with numpy
#
#x, y = np.meshgrid(x, y)
#
##convert intensity to numpy array
#
#intensity = np.array(intensity)
#
##plug data into colormesh
#
#plt.pcolor(column_labels, row_labels, intensity)
#plt.colorbar()
##plt determines colorbar itself but we can customize
plt.show()
