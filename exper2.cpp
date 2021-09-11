#include <iostream>
#include <functional>
#include <random>
#include <time.h>
#include <math.h>
#include <utility>
#include <queue>
#include <bits/stdc++.h>

using namespace std;

long double determinant(long double a, long double b, long double c, long double d, long double e, long double f, long double g, long double h, long double i) {
    return a * (e*i - h*f) - b * (d*i - g*f) + c * (d*h - g*e);
}

long double foo(long double x) {
    return 0.2 * sin(100*x) + 0.1 * sin(40*x) - 2 * (x - 0.6) * (x-0.6) + 0.15 * cos(30*x) + 1;
    //return sin(15*x);
}

long double pdf(long double x) {
    return -3 * (x - 0.5) * (x - 0.5) + 1.25;
}

long double pdf2(long double x) {
    return -1.786 * x * x + 2.099 * x + 0.546;
}

long double cdf(long double x) {
    return -1/8.0 - (x - 0.5) * (x - 0.5) * (x - 0.5) + 1.25 * x;
}

long double cdf2(long double x) {
    return -0.595 * x * x * x + 1.05 * x * x + 0.546 * x;
}

long double invCdf(long double x) {
    long double start = 0;
    long double end = 1.0;
    long double tempY = (start + end) / 2;
    long double tempX = cdf(tempY);
    long double count = 0;
    while (count < 250 && abs(tempX - x) > 0.000001) {
        tempY = (start + end) / 2;
        tempX = cdf(tempY);
        if (tempX > x) {
            end = tempY;
        } else {
            start = tempY;
        }
        count ++;
    }
    return tempY;
}

long double invCdf2(long double x) {
    long double start = 0;
    long double end = 1.0;
    long double tempY = (start + end) / 2;
    long double tempX = cdf2(tempY);
    long double count = 0;
    while (abs(tempX - x) > 0.000001) {
        tempY = (start + end) / 2;
        tempX = cdf2(tempY);
        if (tempX > x) {
            end = tempY;
        } else {
            start = tempY;
        }
        count ++;
    }
    return tempY;
}

long double foo2(long double x) {
    return foo(invCdf(x)) / pdf(invCdf(x));
}

long double foo22(long double x) {
    return foo(invCdf2(x)) / pdf2(invCdf2(x));
}

int regionsBudget = 1000;
const long double epsilon = 0.00001;

long double miniSimpsons(long double intervalStart, long double intervalEnd, function<long double (long double)> testFunc) {
    return (testFunc(intervalStart) + 4 * testFunc((intervalStart + intervalEnd) / 2.0) + testFunc(intervalEnd)) * (intervalEnd - intervalStart) / 6.0;
}

long double miniTrapezoidal(long double intervalStart, long double intervalEnd, function<long double (long double)> testFunc) {
    return (intervalEnd - intervalStart) * (testFunc(intervalEnd) + testFunc(intervalStart)) / 2.0;
}

long double estimateError(long double intervalStart, long double intervalEnd, function<long double (long double)> testFunc) {
    return abs(miniSimpsons(intervalStart, intervalEnd, testFunc) - miniTrapezoidal(intervalStart, intervalEnd, testFunc)) + epsilon * (intervalEnd - intervalStart);
}

void getCoefs(long double x1, long double x3, function<long double (long double)> testFunc, long double* coefs) {
    //quadratic approximation
    long double x2 = (x1 + x3) / 2.0;
    long double y1 = testFunc(x1);
    long double y2 = testFunc(x2);
    long double y3 = testFunc(x3);

    long double D = determinant(x1*x1, x1, 1, x2*x2, x2, 1, x3*x3, x3, 1);
    long double Da = determinant(y1, x1, 1, y2, x2, 1, y3, x3, 1);
    long double Db = determinant(x1*x1, y1, 1, x2*x2, y2, 1, x3*x3, y3, 1);
    long double Dc = determinant(x1*x1, x1, y1, x2*x2, x2, y2, x3*x3, x3, y3);
    
    long double a = Da / D;
    long double b = Db / D;
    long double c = Dc / D;

    //cout << x1 << " to " << x3 << endl;
    //cout << a << ", " << b << ", " << c << endl;
    coefs[0] = a;
    coefs[1] = b;
    coefs[2] = c;
}

long double a = 0;
long double b = 1;

int compareArr(const void * arr1, const void * arr2) {
    //cout << (*((long double**) arr1))[0] << " vs " << *((long double**) arr2)[0] << endl;
    //cout << "Result: " << ((long double**) arr1)[0][0] - ((long double**) arr2)[0][0] << endl;
    return (((long double**) arr1)[0][0] > ((long double**) arr2)[0][0]) - (((long double**) arr1)[0][0] < ((long double**) arr2)[0][0]);
}

int compareArr2(const void * arr1, const void * arr2) {
    //cout << ((long double*) arr2)[0] << " vs " << ((long double*) arr1)[0] << endl;
    //cout << *((long double*) arr2) << endl;
    return (*((long double*) arr1) > *((long double*) arr2)) - (*((long double*) arr2) > *((long double*) arr1));
}

long double** splitRegions(function<long double (long double)> testFunc) {
    long double** regions = new long double*[regionsBudget];
    for (int i = 0; i < regionsBudget; i++) {
        regions[i] = new long double[5];
    }

    int numRegions = 1;
    tuple<long double, long double, long double> curr;
    priority_queue<tuple<long double, long double, long double>> myHeap;
    //Each tuple has error estimate, startX, stopX
    myHeap.push(tuple<long double, long double, long double>(estimateError(a, b, testFunc),a,b));
    long double start, end, mid;
    while (numRegions < regionsBudget) {
        curr = myHeap.top();
        start = get<1>(curr);
        end = get<2>(curr);
        mid = (start + end) / 2;
        myHeap.pop();
        myHeap.push(tuple<long double, long double, long double>(estimateError(start, mid, testFunc), start, mid));
        myHeap.push(tuple<long double, long double, long double>(estimateError(mid, end, testFunc), mid, end));
        numRegions++;
    }

    for (int i = 0; i < regionsBudget; i++) {
        curr = myHeap.top();
        regions[i][0] = get<1>(curr);
        regions[i][1] = get<2>(curr);
        myHeap.pop();
        getCoefs(regions[i][0], regions[i][1], testFunc, &regions[i][2]);
        //cout << regions[i][2] << ", " << regions[i][3] << ", " << regions[i][4] << endl;
    }

    qsort(regions, regionsBudget, sizeof(regions[0]), compareArr);

    return regions;
}

long double polyApprox(long double x, long double** regions) {
    int low = 0;
    int high = regionsBudget;
    int curr = (low + high) / 2;
    while (regions[curr][0] > x || regions[curr][1] < x) {
        if (regions[curr][0] > x) {
            high = curr;
            curr = (low + high) / 2;
        } else if (regions[curr][1] < x) {
            low = curr;
            curr = (low + high) / 2;
        }
    }
    return regions[curr][2] * x * x + regions[curr][3] * x + regions[curr][4];
}

long double polyInt(long double start, long double end, long double a, long double b, long double c) {
    a /= 3;
    b /= 2;
    long double integral = 0;
    integral += a * end * end * end;
    integral += b * end * end;
    integral += c * end;
    integral -= a * start * start * start;
    integral -= b * start * start;
    integral -= c * start;
    if (start < 0.00001) {
        //cout << "Estimating\n";
        //cout << "Start: " << start << ", end: " << end << ", a: " << a << ", b: " << b << ", c: " << c << endl;
        //cout << integral << endl;
    }
    return integral;
}

long double polyApproxEst(long double** regions) {
    long double estimate = 0;
    //cout << "Hmmm " << regions[0][2] << endl;
    for (int i = 0; i < regionsBudget; i++) {
        //cout << polyInt(regions[i][0], regions[i][1], regions[i][2], regions[i][3], regions[i][4]) << endl;
        estimate += polyInt(regions[i][0], regions[i][1], regions[i][2], regions[i][3], regions[i][4]);
        if (i > regionsBudget - 50) {
            //cout << estimate << endl;
        }
        //estimate += (regions[i][2] * regions[i][1] * regions[i][1] + regions[i][3] * regions[i][1] + regions[i][4]) - (regions[i][2] * regions[i][0] * regions[i][0] + regions[i][3] * regions[i][0] + regions[i][4]);
    }
    return estimate;
}

long double groundTruth = 0.812835882622424311126150;//From wolfram alpha

int main() {
    cout.precision(12);
    //getCoefs(0, 1, foo);
    //cout << determinant(2, 3, 7, -5, 23, 10, 4, 5, -1);
    regionsBudget = 10000;
    long double** myRegions = splitRegions(foo22);

    /*for (int i = 0; i < regionsBudget; i++) {
        cout << "[";
        for (int j = 0; j < 4; j++) {
            cout << myRegions[i][j] << ", ";
        }
        cout << myRegions[i][4];
        cout << "]," << endl;
    }*/

    //cout << polyApproxEst(myRegions) << endl;
    for (double j = 1.5; j < 15; j += 0.5) {
        if ((int) exp(j) == (int) exp(j-0.05)) {
            continue;
        }
        regionsBudget = (int) exp(j);
        //cout << regionsBudget << endl;
        myRegions = splitRegions(foo);
        cout << ((int) exp(j)) << ": " << polyApproxEst(myRegions) << endl;
        //cout << abs(polyApproxEst(myRegions)-groundTruth) << ",";
    }
    cout << endl;

    delete myRegions;
    myRegions = nullptr;
}