import matplotlib.pyplot as plt
import numpy as np

N, mc, un, st, ha = np.loadtxt('conv.txt', delimiter=',', unpack=True)

plt.plot(N, mc, linestyle='-', marker='.', label='Pure Monte Carlo')
plt.plot(N, un, linestyle='-', marker='.', label='Uniform')
plt.plot(N, st, linestyle='-', marker='.', label='Stratified')
plt.plot(N, ha, linestyle='-', marker='.', label='Halton')

plt.xlabel('Number of Samples')
plt.ylabel('RMSE')
plt.title('1D convergence graphs')
plt.legend()
plt.show()

logN = np.log(N, out=np.zeros_like(N), where=(N!=0))
logMc = np.log(mc, out=np.zeros_like(mc), where=(mc!=0))
logUn = np.log(un, out=np.zeros_like(un), where=(un!=0))
logSt = np.log(st, out=np.zeros_like(st), where=(st!=0))
logHa = np.log(ha, out=np.zeros_like(ha), where=(ha!=0))

m, b = np.polyfit(logN, logMc, 1)
plt.plot(logN, logMc, linestyle='-', marker='.', label='Pure Monte Carlo (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logUn, 1)
plt.plot(logN, logUn, linestyle='-', marker='.', label='Uniform (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logSt, 1)
plt.plot(logN, logSt, linestyle='-', marker='.', label='Stratified (approx rate: '+"{:.4f}".format(m)+')')
m, b = np.polyfit(logN, logHa, 1)
plt.plot(logN, logHa, linestyle='-', marker='.', label='Halton (approx rate: '+"{:.4f}".format(m)+')')

plt.xlabel('log(Number of Samples)')
plt.ylabel('log(RMSE)')
plt.title('1D convergence graphs')
plt.legend()
plt.show()