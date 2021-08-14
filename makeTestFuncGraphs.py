import matplotlib.pyplot as plt
import numpy as np

def gaussianDerivativeWRTMean1D(x, myLambda):
    return np.e ** (-0.5 * (x - myLambda) ** 2) * (x - myLambda) / ((2 * np.pi) ** 0.5)



plt.clf()

x = np.linspace(0.0, 1.0, 100)
y = gaussianDerivativeWRTMean1D(x, 0)

samplesold = np.array([0.0731697696465,0.00343649366986,0.521980644646,0.822683116941,0.30683195669,0.69316804331,0.177316883059,0.478019355354,0.99656350633,0.926830230353])
samples = np.array([0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95])
sampleYs = gaussianDerivativeWRTMean1D(samples,0)
plt.plot(x, y)
plt.scatter(samples, sampleYs)

for i in range(len(samples)):
    plt.axhline(y=sampleYs[i], linestyle='--', color='red')

avg = 0.156971555882
guess = 0.150279908577
guess = 0.157137986409

plt.axhline(y=avg, linestyle='-', color='green')
plt.axhline(y=guess, linestyle='-', color='purple')

plt.show()

plt.savefig('temp2.png')

