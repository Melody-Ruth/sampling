import matplotlib.pyplot as plt
import numpy as np

N, mc, an, un, unji, st, stan, stan2, stan3, ha, haro, haan = np.loadtxt('conv.txt', delimiter=',', unpack=True)

logN = np.log(N, out=np.zeros_like(N), where=(N!=0))
logMc = np.log(mc, out=np.zeros_like(mc), where=(mc!=0))
logUn = np.log(un, out=np.zeros_like(un), where=(un!=0))
logUnji = np.log(unji, out=np.zeros_like(unji), where=(unji!=0))
logAn = np.log(an, out=np.zeros_like(an), where=(an!=0))
logSt = np.log(st, out=np.zeros_like(st), where=(st!=0))
logStan = np.log(stan, out=np.zeros_like(stan), where=(stan!=0))
logStan2 = np.log(stan2, out=np.zeros_like(stan2), where=(stan2!=0))
logStan3 = np.log(stan3, out=np.zeros_like(stan3), where=(stan3!=0))
logHa = np.log(ha, out=np.zeros_like(ha), where=(ha!=0))
logHaro = np.log(haro, out=np.zeros_like(haro), where=(haro!=0))
logHaan = np.log(haan, out=np.zeros_like(haan), where=(haan!=0))

m, b = np.polyfit(logN, logMc, 1)
plt.plot(N, mc, linestyle='-', marker='.', label='Pure Monte Carlo (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logAn, 1)
plt.plot(N, an, linestyle='-', marker='.', label='Antithetic Monte Carlo (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logUn, 1)
plt.plot(N, un, linestyle='-', marker='.', label='Uniform (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logUnji, 1)
plt.plot(N, unji, linestyle='-', marker='.', label='Uniform Jitter (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logSt, 1)
plt.plot(N, st, linestyle='-', marker='.', label='Stratified (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logStan, 1)
plt.plot(N, stan, linestyle='-', marker='.', label='Stratified, antithetic (local) (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logStan2, 1)
plt.plot(N, stan2, linestyle='-', marker='.', label='Stratified, antithetic (global) (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logStan3, 1)
plt.plot(N, stan3, linestyle='-', marker='.', label='Stratified, antithetic (global and local) (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logHa, 1)
#plt.plot(N, ha, linestyle='-', marker='.', label='Halton (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logHaro, 1)
#plt.plot(N, haro, linestyle='-', marker='.', label='Halton, rotated (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logHaan, 1)
#plt.plot(N, haan, linestyle='-', marker='.', label='Halton, antithetic, (approx rate: '+"{:.4f}".format(m)+')')
plt.yscale('log')
plt.xscale('log')

plt.xlabel('Log(Number of Samples)')
plt.ylabel('Log(RMSE)')
plt.title('1D convergence graphs')
plt.legend()
plt.show()

plt.savefig('tempConvGraph.png')

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