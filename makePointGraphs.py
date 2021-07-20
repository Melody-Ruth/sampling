import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt('points.txt', delimiter=',', unpack=True)

plt.plot(x, y, '.', label='Stratified')
for i in range(17):
    plt.axvline(x=i/16, linestyle='--', color='red')
#plt.yscale('log')
#plt.xscale('log')

plt.title('16 points')
plt.legend()
plt.show()

plt.savefig('tempPoints.png')