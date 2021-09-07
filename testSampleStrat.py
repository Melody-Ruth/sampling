import numpy as np
import matplotlib.pyplot as plt
import math

def w1Formula(x1, x2, x3):
    return (0.33333 - (x2 * x2)/(2 * (x2 - x3)) + (x3 * x2 * x2) / (x2 - x3) - x3 * x3 + x3 * x3 / (2 * (x2 - x3)) - x3 * x3 * x3 / (x2 - x3)) / (x1 * x1 - x3 * x3 + (x1 - x3) * (x3 * x3 - x2 * x2) / (x2 - x3))

def foo(x):
    return 0.2 * math.sin(100*x) + 0.1 * math.sin(40*x) - 2 * (x - 0.6) * (x-0.6) + 0.15 * math.cos(30*x) + 1

def pdf(x):
    #return -3 * (x - 0.5) ** 2 + 1.25
    return -1.786 * x * x + 2.099 * x + 0.546

def cdf(x):
    #return -1/8.0 - (x - 0.5) ** 3 + 1.25*x
    return -0.595 * x * x * x + 1.05 * x * x + 0.546 * x

def invCdf(x):
    start = 0
    end = 1.0
    tempY = (start + end) / 2
    tempX = cdf(tempY)
    count = 0
    while count < 250 and abs(tempX - x) > 0.0001:
        tempY = (start + end) / 2
        tempX = cdf(tempY)
        if tempX > x:
            end = tempY
        else:
            start = tempY 
        count += 1
    return tempY

def foo2(x):
    return foo(invCdf(x)) / pdf(invCdf(x))
    #return 3
    #return invCdf(x)


# plt.plot(np.arange(0,1,0.001), foo(np.arange(0,1,0.001)), linestyle='-', label='Foo')
# plt.show()

x = np.linspace(0.0, 1.0, 1000)
#for i in range(0, 1000):
foo_vec = np.vectorize(foo)
foo2_vec = np.vectorize(foo2)
pdf_vec = np.vectorize(pdf)

plt.subplot(1, 2, 1)
y = foo_vec(x)
y2 = pdf_vec(x)
plt.plot(x, y)
plt.plot(x, y2)

plt.subplot(1, 2, 2)
y = foo2_vec(x)
plt.plot(x, y)

x = np.linspace(0.0, 0.25, 100)
y = -2.08897 * x**2 - 0.560147 * x + 0.783074
plt.plot(x, y, color='orange')
plt.axvline(x=0.25, linestyle='--', color='red')

x = np.linspace(0.25, 0.5, 100)
y = -3.38676 * x**2 + 3.61422 * x + -0.169406
plt.plot(x, y, color='orange')
plt.axvline(x=0.5, linestyle='--', color='red')

x = np.linspace(0.5, 0.625, 100)
y = 8.45619 * x**2 + -6.32399 * x + 1.83896
plt.plot(x, y, color='orange')
plt.axvline(x=0.625, linestyle='--', color='red')

x = np.linspace(0.625, 0.75, 100)
y = -6.11312 * x**2 + 3.64454 * x + 1.29977
plt.plot(x, y, color='orange')
plt.axvline(x=0.75, linestyle='--', color='red')
#plt.axvline(x=0.78125, linestyle='--', color='red')
#plt.axvline(x=0.8125, linestyle='--', color='red')

x = np.linspace(0.75, 0.875, 100)
y = -50.6355 * x**2 + 83.3905 * x + -33.4658
plt.plot(x, y, color='orange')
plt.axvline(x=0.875, linestyle='--', color='red')

x = np.linspace(0.875, 0.9375, 100)
y = -317.799 * x**2 + 571.103 * x + -255.667
plt.plot(x, y, color='orange')
plt.axvline(x=0.9375, linestyle='--', color='red')

x = np.linspace(0.9375, 1, 100)
y = -360.013 * x**2 + 702.895 * x + -342.12
plt.plot(x, y, color='orange')
plt.show()

# for x1 in np.arange(0.0,1.0,0.1):
#    for x2 in np.arange(0.0,1.0,0.1):
#        for x3 in np.arange(0.0,1.0,0.1):
#            if (w1Formula(x1,x2,x3) < 0):
#                print(w1Formula(x1,x2,x3))
