import matplotlib.pyplot as plt
import numpy as np

x = [1, 2, 3]
y = [1, 2, 3]

intensity = [
    [1, 1, 1],
    [2, 2, 2],
    ]

# set up 2D grid with numpy

x, y = np.meshgrid(x, y)

#convert intensity to numpy array

intensity = np.array(intensity)

#plug data into colormesh

plt.pcolormesh(x, y, intensity)
plt.colorbar()
#plt determines colorbar itself but we can customize
plt.show()
