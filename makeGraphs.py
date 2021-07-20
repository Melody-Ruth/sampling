import matplotlib.pyplot as plt
import numpy as np

N, mc, an, un, st, ha = np.loadtxt('conv.txt', delimiter=',', unpack=True)

logN = np.log(N, out=np.zeros_like(N), where=(N!=0))
logMc = np.log(mc, out=np.zeros_like(mc), where=(mc!=0))
logUn = np.log(un, out=np.zeros_like(un), where=(un!=0))
logAn = np.log(an, out=np.zeros_like(an), where=(an!=0))
logSt = np.log(st, out=np.zeros_like(st), where=(st!=0))
logHa = np.log(ha, out=np.zeros_like(ha), where=(ha!=0))

m, b = np.polyfit(logN, logMc, 1)
plt.plot(N, mc, linestyle='-', marker='.', label='Pure Monte Carlo (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logAn, 1)
plt.plot(N, an, linestyle='-', marker='.', label='Antithetic Monte Carlo (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logUn, 1)
plt.plot(N, un, linestyle='-', marker='.', label='Uniform (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logSt, 1)
plt.plot(N, st, linestyle='-', marker='.', label='Stratified (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logHa, 1)
plt.plot(N, ha, linestyle='-', marker='.', label='Halton (approx rate: '+"{:.4f}".format(m)+')')
plt.yscale('log')
plt.xscale('log')

plt.xlabel('Log(Number of Samples)')
plt.ylabel('Log(RMSE)')
plt.title('1D convergence graphs')
plt.legend()
plt.show()

def sampleQuad(x, myLambda):
    return 600 * (x-myLambda) * (x-myLambda) + 1

def sampleStep(x, myLambda):
    return 150 * (x > myLambda)

def gaussianDerivativeWRTMean1D(x, myLambda):
    return np.e ** (-0.5 * (x - myLambda) ** 2) * (x - myLambda) / ((2 * np.pi) ** 2)

x = np.linspace(0.0, 3.0, 100)
#for i in range(0, 30):
#    y = gaussianDerivativeWRTMean1D(x, i/10)
#
#    plt.plot(x, y)
#    plt.show()

x = np.linspace(-8.0, 8.0, 100)
y = gaussianDerivativeWRTMean1D(x, 0.5)
plt.plot(x, y)
plt.show()