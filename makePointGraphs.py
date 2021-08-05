import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt('points.txt', delimiter=',', unpack=True)

plt.plot(x, y, '.', label='Globally Antithetic Stratified')
for i in range(6):
    plt.axvline(x=i/5, linestyle='--', color='red')

for i in range(6):
    plt.axhline(y=i/5, linestyle='--', color='red')
    #plt.axline((0,i/12),(1,i/12), linestyle='--', color='red')
#plt.yscale('log')
#plt.xscale('log')

plt.title('50 points')
plt.legend()
plt.show()

plt.savefig('tempPoints.png')