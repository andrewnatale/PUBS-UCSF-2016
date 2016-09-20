import matplotlib.pyplot as plt
import numpy as np

def plot(x, y, freq)


    # set up 2D grid with numpy

    x, y = np.meshgrid(x, y)

    #convert intensity to numpy array

    freq = np.array(freq)


    #plug data into colormesh
    
    plt.pcolormesh(x, y, freq)
    plt.colorbar()
    #plt determines colorbar itself but we can customize
    plt.show()
