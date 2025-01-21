import matplotlib.pyplot as plt
import numpy as np



CSV_PATH = 'out.csv'

def plotimg():
    data = np.genfromtxt(CSV_PATH, delimiter=',')
    plt.imshow(data, cmap='gray')
    plt.savefig('out.png')

    plt.show()


plotimg()