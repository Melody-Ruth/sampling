import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt('points.txt', delimiter=',', unpack=True)

plt.plot(x, y, '.', label='Halton')
#plt.yscale('log')
#plt.xscale('log')

plt.title('First 500 points')
plt.legend()
plt.show()

plt.savefig('tempConvGraph.png')