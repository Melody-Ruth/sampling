import numpy as np
import matplotlib.pyplot as plt
import math

def w1Formula(x1, x2, x3):
    return (0.33333 - (x2 * x2)/(2 * (x2 - x3)) + (x3 * x2 * x2) / (x2 - x3) - x3 * x3 + x3 * x3 / (2 * (x2 - x3)) - x3 * x3 * x3 / (x2 - x3)) / (x1 * x1 - x3 * x3 + (x1 - x3) * (x3 * x3 - x2 * x2) / (x2 - x3))

def foo(x):
    return 0.2 * math.sin(100*x) + 0.1 * math.sin(40*x) - 2 * (x - 0.6) * (x-0.6) + 0.15 * math.cos(30*x) + 1

def pdf(x):
    return -3 * (x - 0.5) ** 2 + 1.25

def cdf(x):
    return -1/8.0 - (x - 0.5) ** 3 + 1.25*x

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
plt.show()

# for x1 in np.arange(0.0,1.0,0.1):
#    for x2 in np.arange(0.0,1.0,0.1):
#        for x3 in np.arange(0.0,1.0,0.1):
#            if (w1Formula(x1,x2,x3) < 0):
#                print(w1Formula(x1,x2,x3))
