import numpy as np
import matplotlib.pyplot as plt
import math

regionsBudget = 25

def w1Formula(x1, x2, x3):
    return (0.33333 - (x2 * x2)/(2 * (x2 - x3)) + (x3 * x2 * x2) / (x2 - x3) - x3 * x3 + x3 * x3 / (2 * (x2 - x3)) - x3 * x3 * x3 / (x2 - x3)) / (x1 * x1 - x3 * x3 + (x1 - x3) * (x3 * x3 - x2 * x2) / (x2 - x3))

def foo(x):
    return 0.2 * math.sin(100*x) + 0.1 * math.sin(40*x) - 2 * (x - 0.6) * (x-0.6) + 0.15 * math.cos(30*x) + 1

def pdf(x):
    #return -3 * (x - 0.5) ** 2 + 1.25
    return -1.786 * x * x + 2.099 * x + 0.546
    #return -1.4067136589943 * x * x + 1.65308956 * x + 0.43

def cdf(x):
    #return -1/8.0 - (x - 0.5) ** 3 + 1.25*x
    return -0.595 * x * x * x + 1.05 * x * x + 0.546 * x
    #return -0.468904553 * x * x * x + 0.82654478 * x * x + 0.43 * x

def invCdf(x):
    start = 0
    end = 1.0
    tempY = (start + end) / 2
    tempX = cdf(tempY)
    while abs(tempX - x) > 0.0001:
        tempY = (start + end) / 2
        tempX = cdf(tempY)
        if tempX > x:
            end = tempY
        else:
            start = tempY 
    return tempY

def primaryFoo(x):
    return foo(invCdf(x)) / pdf(invCdf(x))
    #return 3
    #return invCdf(x)

fooCoefs = [0.074100687227980882588,65.07295006044863328,0.32985005494542579685,2.3629179447093040523,0.58374381375919270898,54.776664887113206248,0.030972167439347768336,24.187096653964168524]
def foo2(x):
    return fooCoefs[0] * math.sin(fooCoefs[1] * x) + fooCoefs[2] * math.sin(fooCoefs[3] * x) + fooCoefs[4] * math.cos(fooCoefs[5] * x) + fooCoefs[6] * math.cos(fooCoefs[7] * x)

def determinant(a, b, c, d, e, f, g, h, i):
    return a * (e*i - h*f) - b * (d*i - g*f) + c * (d*h - g*e)

def polyInt(start, end, a, b, c):
    a /= 3
    b /= 2
    integral = 0
    integral += a * end * end * end
    integral += b * end * end
    integral += c * end
    integral -= a * start * start * start
    integral -= b * start * start
    integral -= c * start
    return integral

def getCoefs(x1, x3, testFunc, coefs):
    x2 = (x1 + x3) / 2.0
    y1 = testFunc(x1)
    y2 = testFunc(x2)
    y3 = testFunc(x3)

    D = determinant(x1*x1, x1, 1, x2*x2, x2, 1, x3*x3, x3, 1)
    Da = determinant(y1, x1, 1, y2, x2, 1, y3, x3, 1)
    Db = determinant(x1*x1, y1, 1, x2*x2, y2, 1, x3*x3, y3, 1)
    Dc = determinant(x1*x1, x1, y1, x2*x2, x2, y2, x3*x3, x3, y3)

    a = 0
    b = 0
    c = 0
    if (D == 0):
        a = 0
        b = 0
        c = y2
    else:
        a = Da / D
        b = Db / D
        c = Dc / D
    coefs[0] = a
    coefs[1] = b
    coefs[2] = c

def pdf2(x):
    pdfCoefs = [0,0,0]
    getCoefs(0,1,foo2,pdfCoefs)
    minValX = 0
    if pdfCoefs[0] < 0:
        minValX = 0
        if (pdfCoefs[0] * minValX * minValX + pdfCoefs[1] * minValX + pdfCoefs[2] > pdfCoefs[0] * 1 * 1 + pdfCoefs[1] * 1 + pdfCoefs[2]):
            minValX = 1
    else:
        minValX = -pdfCoefs[1] / (2 * pdfCoefs[0])
    minVal = pdfCoefs[0] * minValX * minValX + pdfCoefs[1] * minValX + pdfCoefs[2]
    if minVal < 0:
        pdfCoefs[2] -= minVal
    tempIntegral = polyInt(0,1,pdfCoefs[0],pdfCoefs[1],pdfCoefs[2])
    #print(tempIntegral)
    pdfCoefs[0] /= tempIntegral
    pdfCoefs[1] /= tempIntegral
    pdfCoefs[2] /= tempIntegral
    return pdfCoefs[0] * x * x + pdfCoefs[1] * x + pdfCoefs[2]

def cdf2(x):
    pdfCoefs = [0,0,0]
    getCoefs(0,1,foo2,pdfCoefs)
    minValX = 0
    if pdfCoefs[0] < 0:
        minValX = 0
        if (pdfCoefs[0] * minValX * minValX + pdfCoefs[1] * minValX + pdfCoefs[2] > pdfCoefs[0] * 1 * 1 + pdfCoefs[1] * 1 + pdfCoefs[2]):
            minValX = 1
    else:
        minValX = -pdfCoefs[1] / (2 * pdfCoefs[0])
    minVal = pdfCoefs[0] * minValX * minValX + pdfCoefs[1] * minValX + pdfCoefs[2]
    if minVal < 0:
        pdfCoefs[2] -= minVal
    tempIntegral = polyInt(0,1,pdfCoefs[0],pdfCoefs[1],pdfCoefs[2])
    pdfCoefs[0] /= tempIntegral
    pdfCoefs[1] /= tempIntegral
    pdfCoefs[2] /= tempIntegral
    return polyInt(0,x,pdfCoefs[0],pdfCoefs[1],pdfCoefs[2])

def invCdf2(x):
    start = 0
    end = 1.0
    tempY = (start + end) / 2
    tempX = cdf2(tempY)
    while abs(tempX - x) > 0.0001:
        tempY = (start + end) / 2
        tempX = cdf2(tempY)
        if (tempX > x):
            end = tempY
        else:
            start = tempY
    return tempY

def primaryFoo2(x):
    return foo2(invCdf2(x)) / pdf2(invCdf2(x))


def polyApprox(x, regions):
    low = 0
    high = regionsBudget
    curr =  (low + high) // 2
    #print(len(regions))

    while (regions[curr][0] > x or regions[curr][1] < x):
        if (regions[curr][0] > x):
            high = curr
            curr = (low + high) // 2
        elif (regions[curr][1] < x):
            low = curr
            curr = (low + high) // 2
        #print(curr)
    return regions[curr][2] * x * x + regions[curr][3] * x + regions[curr][4]


def residual(x, regions):
    return foo(x) - polyApprox(x,regions)

myRegions = [[0, 0.015625, -760.82259772717770829, 29.765268646910840289, 0.43000000000000004885],
[0.015625, 0.03125, -797.39287972821603034, 26.542582185135586315, 0.48928276746939758368],
[0.03125, 0.0625, 725.42187799459520647, -71.514401269124576288, 2.0664397135763459268],   
[0.0625, 0.125, 77.027572666716707772, -16.266022766368358532, 1.1462063123411083687],     
[0.125, 0.1875, -32.053875241308446675, 18.589082911561717992, -1.5062842738372586027],    
[0.1875, 0.203125, -741.23520629622180422, 307.15845205950381569, -30.680884819177094869], 
[0.203125, 0.21875, -878.09219677353394218, 358.39065656119782943, -35.440739007347019651],
[0.21875, 0.25, 751.68989904346290132, -359.17098079029409519, 43.538249343923837387],
[0.25, 0.3125, 39.462822638684201593, -23.492477234152403298, 4.1328157301871044638],
[0.3125, 0.375, -101.68767304377433902, 74.978508878005555971, -12.85513933587218105],
[0.375, 0.5, 0.67580945989917040606, -1.0481821109799369651, 1.2600050579183381672],
[0.5, 0.5625, -102.28420328698916819, 107.49012566112310196, -27.269145641411114411],
[0.5625, 0.625, -92.017014324953379401, 113.09299054874014701, -33.669359898214793247],
[0.625, 0.640625, -495.54044078520382755, 646.54640237262356095, -209.45140382710360427],
[0.640625, 0.65625, -973.71902740513905883, 1255.0993035363062518, -403.06092739776613598],
[0.65625, 0.6875, 665.66021622658900014, -904.06387040629397234, 307.86818046145629069],
[0.6875, 0.75, -65.406481926851483877, 89.140436674120337557, -29.41653660724148267],
[0.75, 0.78125, -656.09452046864635122, 1015.2726134374547655, -391.75364749998254865],
[0.78125, 0.8125, 682.08177934967534384, -1085.6714573258250311, 432.85091228924761708],
[0.8125, 0.828125, -464.01344212042749859, 780.23626513464841992, -326.59718866126149805],
[0.828125, 0.84375, -1035.7941479434375651, 1723.265651754787541, -715.42132849491554225],
[0.84375, 0.875, 660.41389718503705808, -1148.0656966409487723, 499.71013646138180775],
[0.875, 0.9375, -92.797135289142261172, 163.8194472433278861, -71.512167699316535163],
[0.9375, 0.96875, -681.63924146393856063, 1308.9388716773464694, -627.52462072601724685],
[0.96875, 1, 645.39939341336685175, -1274.8247153639649696, 630.10169785590778702]]

# plt.plot(np.arange(0,1,0.001), foo(np.arange(0,1,0.001)), linestyle='-', label='Foo')
# plt.show()

x = np.linspace(0.0, 1.0, 1000)
#for i in range(0, 1000):
foo_vec = np.vectorize(foo)
foo2_vec = np.vectorize(primaryFoo)
#foo22_vec = np.vectorize(primaryFoo2)
pdf_vec = np.vectorize(pdf)
#pdf2_vec = np.vectorize(pdf2)

plt.subplot(1, 2, 1)
y = [0] * 1000
y2 = [0] * 1000
for i in range(0,1000):
    y[i] = foo(x[i])
    y2[i] = pdf(x[i])
    #y2[i] = invCdf(x[i])
plt.plot(x, y)
plt.plot(x, y2)

plt.subplot(1, 2, 2)
y = [0] * 1000
y2 = [0] * 1000
plt.xticks(np.linspace(0.0, 1.0, 11))
for i in range(0,1000):
    #y[i] = foo(invCdf(x[i]))
    #y[i] = foo2(x[i])
    #y2[i] = pdf(invCdf(x[i]))
    y[i] = foo(x[i])
    y2[i] = polyApprox(x[i],myRegions)
plt.plot(x, y)
plt.plot(x, y2)
plt.show()

#plt.subplot(1, 2, 1)
y = foo_vec(x)
y2 = pdf_vec(x)
#plt.plot(x, y)
#plt.plot(x, y2)

est = 0
x = np.linspace(0.0, 1.0, 10000)
y = [0] * 10000
y2 = [0] * 10000
#fooCoefs = [0.34252596485104030988,40.373478800639617248,0.41481951205383099657,1.5401418457083917435,0.90731319362333173739,20.79972175384445876,0.65148028737933372234,51.340066572970592063]
#fooCoefs = [0.702400453854602147,50.504023186394000788,0.45295555709445656234,62.096501157973911234,0.22370818963690680681,24.208689887070761415,0.87206685308697196035,47.328547906014968305]
#fooCoefs = [0.60552367758166736333,10.994924155005476152,0.066656109992848516788,1.1645257770499917171,0.76824161565090720583,35.077295830035957636,0.82686686508002282814,7.6933851758133062759]
for i in range(0,10000):
    #y[i] = foo2(x[i])
    y[i] = residual(x[i],myRegions)
    est += foo2(x[i])
    y2[i] = pdf2(x[i])
est /= 10000
print(est)
plt.xticks(np.linspace(0.0, 1.0, 21))
plt.plot(x, y)
#plt.plot(x, y2)

#plt.subplot(1, 2, 2)
#y = foo22_vec(x)
#y = foo_vec(x)
y = [0] * 10000
y2 = [0] * 10000
est = 0
#for i in range(0,10000):
    #y[i] = foo2(invCdf2(x[i]))
    #y[i] = foo2(x[i])
    #y2[i] = pdf2(invCdf2(x[i]))
    #y2[i] = primaryFoo2(x[i])
    #est += primaryFoo2(x[i])
#est /= 10000
#print(est)
#plt.plot(x, y)
#plt.plot(x, y2)
plt.show()

print(foo(invCdf(0.9)))
print(pdf(invCdf(0.9)))
print(invCdf(0.9))

x = np.linspace(0.0, 0.25, 100)
y = -2.08897 * x**2 - 0.560147 * x + 0.783074
#plt.plot(x, y, color='orange')
#plt.axvline(x=0.25, linestyle='--', color='red')

x = np.linspace(0.25, 0.5, 100)
y = -3.38676 * x**2 + 3.61422 * x + -0.169406
#plt.plot(x, y, color='orange')
#plt.axvline(x=0.5, linestyle='--', color='red')

x = np.linspace(0.5, 0.625, 100)
y = 8.45619 * x**2 + -6.32399 * x + 1.83896
#plt.plot(x, y, color='orange')
#plt.axvline(x=0.625, linestyle='--', color='red')

x = np.linspace(0.625, 0.75, 100)
y = -6.11312 * x**2 + 3.64454 * x + 1.29977
#plt.plot(x, y, color='orange')
#plt.axvline(x=0.75, linestyle='--', color='red')

x = np.linspace(0.75, 0.875, 100)
y = -50.6355 * x**2 + 83.3905 * x + -33.4658
#plt.plot(x, y, color='orange')
#plt.axvline(x=0.875, linestyle='--', color='red')

x = np.linspace(0.875, 0.9375, 100)
y = -317.799 * x**2 + 571.103 * x + -255.667
#plt.plot(x, y, color='orange')
#plt.axvline(x=0.9375, linestyle='--', color='red')

x = np.linspace(0.9375, 1, 100)
y = -360.013 * x**2 + 702.895 * x + -342.12
#plt.plot(x, y, color='orange'

x2 = np.linspace(0, 1, 1000)
y2 = 1000 * [0]
#for i in range(1000):
#    y2[i] = polyApprox(x2[i], myRegions)

#for i in range(regionsBudget):
    #plt.axvline(x=myRegions[i][1], linestyle='--', color='red')

#y = polyApprox_vec(x, myRegions)
#plt.plot(x2, y2, color='red')

#plt.show()

x3 = [11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,55,57,59,63,67,69,73,77,81,85,89,93,99,103,109,115,121,127,133,141,147,155,163,171,181,189,199,209,219,231,243,255,269,283,297,313,329,345,363,381,401,421,443,465,489,515,541,569,597,629,661,695,731,767,807,849,891,937,985,1037,1089,1145,1203,1265,1331,1399,1471,1545,1625,1709,1795,1887,1985,2087,2193,2305,2423,2549,2679,2817,2961,3113,3271,3439,3617]
y3 = [0.41602118697101647671,0.33079507385298700939,0.34116705206698655362,0.28545825372269445497,0.2787970634902048328,0.22932474551458936762,0.23248299136339420357,0.18939169182180110029,0.18860963200668459661,0.15825997025812754693,0.16128216432382172507,0.13153095977781043002,0.12710719859010541355,0.11119756270141292753,0.11598242336709861655,0.095053296667747505522,0.10017444162020558429,0.098062235681958737077,0.088517370400947065212,0.075054501556389374706,0.073344034159775084447,0.066094056516837598747,0.054019821894521251326,0.056741561705405917093,0.039768606383962457063,0.022181794895447588811,0.018987269203894235892,0.0202298617983969857,0.017078327757643323864,0.017802682449896082423,0.015884257539955079497,0.015489572861679071908,0.013299746983093141059,0.013827384812109336232,0.012283652231831223428,0.010148510286803831895,0.0093038990946695768885,0.0070903520969497584037,0.0052434791631256298661,0.0022460246923352174785,0.0020420633418972550764,0.0018122182248165913301,0.0017499153605105254028,0.0017356290583200486947,0.0015090874031076868724,0.0014903186357050875535,0.0013298585758761169266,0.0012516082895643793928,0.0010524562202622479018,0.00095212024393183109904,0.00082826239600724443132,0.0006270146606538239218,0.00027682020403636299621,0.00021164978894384531803,0.00019317287388917548708,0.00018892034686190988186,0.00017441105306158032453,0.00015635501570878442825,0.00014001520058958431066,0.00013636094735690578849,0.00012069339615233080627,0.00011311811443104350858,0.00010261425597091075818,8.278675805517965306e-05,6.5686758948687391704e-05,4.663015731515285367e-05,1.8534037661981141187e-05,1.7668152640269802682e-05,1.6649031354595595652e-05,1.4338117243913883295e-05,1.3604640914906050941e-05,1.4091510696009435459e-05,1.2285179206369736533e-05,1.1262276097169131038e-05,1.0536508603029101316e-05,9.3112256118639750275e-06,8.5032779913979435605e-06,7.5530680857378206528e-06,5.5005609704458891421e-06,4.2909021690090272554e-06,1.7454202965293943049e-06,1.536187183658918145e-06,1.5017387656213667134e-06,1.3576835199268155879e-06,1.2876577266270459314e-06,1.2290758487400862149e-06,1.142076688001970804e-06,1.0296937886303714293e-06,9.9520738158853353135e-07,8.7861824972344413325e-07,7.2575070097551831145e-07,7.170168356873974697e-07,4.154350234776712254e-07,3.3756272844758381226e-07,1.3248591358415374556e-07,1.3538648451508473305e-07,1.5996395850890006362e-07,1.0218360932269399115e-07,1.2166771235501426111e-07,1.1583899292987903617e-07,9.2560173878566708582e-08,9.4524428479597291597e-08,1.0749937932210383969e-07,7.4582321688733049693e-08,9.3394184707564116745e-08,6.2826159415239562758e-08]
plt.scatter(x3, y3)
plt.yscale('log')
plt.xscale('log')
plt.show()

# for x1 in np.arange(0.0,1.0,0.1):
#    for x2 in np.arange(0.0,1.0,0.1):
#        for x3 in np.arange(0.0,1.0,0.1):
#            if (w1Formula(x1,x2,x3) < 0):
#                print(w1Formula(x1,x2,x3))
