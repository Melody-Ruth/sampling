#include <iostream>
#include <functional>
#include <random>
#include <time.h>
#include <math.h>
#include <utility>
#include <queue>
#include <bits/stdc++.h>

using namespace std;

double determinant(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    return a * (e*i - h*f) - b * (d*i - g*f) + c * (d*h - g*e);
}

//Returns the integral from x=start to x=end of f(x) = ax^2 + bx + c
double polyInt(double start, double end, double a, double b, double c) {
    a /= 3;
    b /= 2;
    double integral = 0;
    integral += a * end * end * end;
    integral += b * end * end;
    integral += c * end;
    integral -= a * start * start * start;
    integral -= b * start * start;
    integral -= c * start;
    //if (start < 0.25) {
        //cout << "Estimating\n";
        //cout << "Start: " << start << ", end: " << end << ", a: " << a << ", b: " << b << ", c: " << c << endl;
        //cout << integral << endl;
    //}
    return integral;
}


double fooCoefs[8] = {0.2,100,0.1,40,0.15,30,0.5,20};

double foo(double x) {
    return 0.2 * sin(100*x) + 0.1 * sin(40*x) - 2 * (x - 0.6) * (x-0.6) + 0.15 * cos(30*x) + 1;
    //return sin(15*x);
}

//customizable version of the test function
double foo2(double x) {
    return fooCoefs[0] * sin(fooCoefs[1] * x) + fooCoefs[2] * sin(fooCoefs[3] * x) + fooCoefs[4] * cos(fooCoefs[5] * x) + fooCoefs[6] * cos(fooCoefs[7] * x);
}

long double fooIntegralGroundTruth2() {
    long double result = -fooCoefs[0] * cos(fooCoefs[1] * 1) / fooCoefs[1] - fooCoefs[2] * cos(fooCoefs[3] * 1) / fooCoefs[3] + fooCoefs[4] * sin(fooCoefs[5] * 1) / fooCoefs[5] + fooCoefs[6] * sin(fooCoefs[6] * 1) / fooCoefs[7];
    result -= -fooCoefs[0] * cos(fooCoefs[1] * 0) / fooCoefs[1] - fooCoefs[2] * cos(fooCoefs[3] * 0) / fooCoefs[3] + fooCoefs[4] * sin(fooCoefs[5] * 0) / fooCoefs[5] + fooCoefs[6] * sin(fooCoefs[6] * 0) / fooCoefs[7];
    //cout << "Part 1: " << (-fooCoefs[0] * cos(fooCoefs[1] * 1) / fooCoefs[1] + fooCoefs[0] * cos(fooCoefs[1] * 0) / fooCoefs[1]) << endl;
    //cout << "Part 2: " << (-fooCoefs[2] * cos(fooCoefs[3] * 1) / fooCoefs[3] + fooCoefs[2] * cos(fooCoefs[3] * 0) / fooCoefs[3]) << endl;
    //cout << "Part 3: " << (fooCoefs[4] * sin(fooCoefs[5] * 1) / fooCoefs[5] - fooCoefs[4] * sin(fooCoefs[5] * 0) / fooCoefs[5]) << endl;
    //cout << "Part 4: " << (fooCoefs[6] * sin(fooCoefs[7] * 1) / fooCoefs[7] - fooCoefs[6] * sin(fooCoefs[7] * 0) / fooCoefs[7]) << endl;
    result = 0;
    result += (-fooCoefs[0] * cos(fooCoefs[1] * 1) / fooCoefs[1] + fooCoefs[0] * cos(fooCoefs[1] * 0) / fooCoefs[1]);
    result += (-fooCoefs[2] * cos(fooCoefs[3] * 1) / fooCoefs[3] + fooCoefs[2] * cos(fooCoefs[3] * 0) / fooCoefs[3]);
    result += (fooCoefs[4] * sin(fooCoefs[5] * 1) / fooCoefs[5] - fooCoefs[4] * sin(fooCoefs[5] * 0) / fooCoefs[5]);
    result += (fooCoefs[6] * sin(fooCoefs[7] * 1) / fooCoefs[7] - fooCoefs[6] * sin(fooCoefs[7] * 0) / fooCoefs[7]);
    return result;
}


void getCoefs(double x1, double x3, function<double (double)> testFunc, double* coefs) {
    //quadratic approximation
    double x2 = (x1 + x3) / 2.0;
    double y1 = testFunc(x1);
    double y2 = testFunc(x2);
    double y3 = testFunc(x3);

    double D = determinant(x1*x1, x1, 1, x2*x2, x2, 1, x3*x3, x3, 1);
    double Da = determinant(y1, x1, 1, y2, x2, 1, y3, x3, 1);
    double Db = determinant(x1*x1, y1, 1, x2*x2, y2, 1, x3*x3, y3, 1);
    double Dc = determinant(x1*x1, x1, y1, x2*x2, x2, y2, x3*x3, x3, y3);

    double a;
    double b;
    double c;
    if (D == 0) {
        /*cout << endl << "Oh no!! x1: " << x1 << ", x2: " << x2 << ", x3: " << x3 << ", y1: " << y1 << ", y2: " << y3 << ", y3: " << y3 << endl << endl;
        cout << "Matrix:\n";
        cout << x1*x1 << " " << x1 << " " << 1 << endl << x2*x2 << " " << x2 << " " << 1 << endl << x3*x3 << " " << x3 << " " << 1 << endl;*/
        //Just do a constant line
        //cout << endl << "Oh no!! x1: " << x1 << ", x2: " << x2 << ", x3: " << x3 << ", y1: " << y1 << ", y2: " << y3 << ", y3: " << y3 << endl << endl;
        a = 0;
        b = 0;
        c = y2;
    } else {
        a = Da / D;
        b = Db / D;
        c = Dc / D;
    }

    //cout << x1 << " to " << x3 << endl;
    //cout << a << ", " << b << ", " << c << endl;
    if (isinf(a)) {
        cout << "What????" << endl;
    }
    coefs[0] = a;
    coefs[1] = b;
    coefs[2] = c;
}

double pdf(double x) {
    return -1.786 * x * x + 2.099 * x + 0.546;
}

//customizable version of the pdf
double pdf2(double x) {
    double pdfCoefs[3];
    getCoefs(0,1,foo2,pdfCoefs);
    double tempIntegral = polyInt(0,1,pdfCoefs[0],pdfCoefs[1],pdfCoefs[2]);
    //We want the integral to be 1 (so it's a pdf)
    pdfCoefs[0] /= tempIntegral;
    pdfCoefs[1] /= tempIntegral;
    pdfCoefs[2] /= tempIntegral;
    return pdfCoefs[0] * x * x + pdfCoefs[1] * x + pdfCoefs[2];
}

double cdf(double x) {
    return -0.595 * x * x * x + 1.05 * x * x + 0.546 * x;
}

//customizable version of the cdf
double cdf2(double x) {
    double pdfCoefs[3];
    getCoefs(0,1,foo2,pdfCoefs);
    double tempIntegral = polyInt(0,1,pdfCoefs[0],pdfCoefs[1],pdfCoefs[2]);
    //We want the integral to be 1 (so it's a pdf)
    pdfCoefs[0] /= tempIntegral;
    pdfCoefs[1] /= tempIntegral;
    pdfCoefs[2] /= tempIntegral;
    return polyInt(0,x,pdfCoefs[0],pdfCoefs[1],pdfCoefs[2]);
}

double invCdf(double x) {
    double start = 0;
    double end = 1.0;
    double tempY = (start + end) / 2;
    double tempX = cdf(tempY);
    double count = 0;
    while (count < 250 && abs(tempX - x) > 0.0001) {
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

//customizable version of the inverse cdf
double invCdf2(double x) {
    double start = 0;
    double end = 1.0;
    double tempY = (start + end) / 2;
    double tempX = cdf2(tempY);
    double count = 0;
    while (count < 250 && abs(tempX - x) > 0.0001) {
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

//Foo in primary space
double primaryFoo(double x) {
    return foo(invCdf(x)) / pdf(invCdf(x));
}

//(customizable) Foo in primary space
double primaryFoo2(double x) {
    return foo2(invCdf2(x)) / pdf2(invCdf2(x));
}

int regionsBudget = 30000;
//const double epsilon = 0.00001;
const double epsilon = 0.001;

double miniSimpsons(double intervalStart, double intervalEnd, function<double (double)> testFunc) {
    return (testFunc(intervalStart) + 4 * testFunc((intervalStart + intervalEnd) / 2.0) + testFunc(intervalEnd)) * (intervalEnd - intervalStart) / 6.0;
}

double miniTrapezoidal(double intervalStart, double intervalEnd, function<double (double)> testFunc) {
    return (intervalEnd - intervalStart) * (testFunc(intervalEnd) + testFunc(intervalStart)) / 2.0;
}

double estimateError(double intervalStart, double intervalEnd, function<double (double)> testFunc) {
    //return abs(miniSimpsons(intervalStart, intervalEnd, testFunc) - miniTrapezoidal(intervalStart, intervalEnd, testFunc)) + epsilon * (intervalEnd - intervalStart);
    return (intervalEnd - intervalStart);
}

double a = 0;
double b = 1;

int compareArr(const void * arr1, const void * arr2) {
    //cout << (*((double**) arr1))[0] << " vs " << *((double**) arr2)[0] << endl;
    //cout << "Result: " << ((double**) arr1)[0][0] - ((double**) arr2)[0][0] << endl;
    return (((double**) arr1)[0][0] > ((double**) arr2)[0][0]) - (((double**) arr1)[0][0] < ((double**) arr2)[0][0]);
}

int compareArr2(const void * arr1, const void * arr2) {
    //cout << ((double*) arr2)[0] << " vs " << ((double*) arr1)[0] << endl;
    //cout << *((double*) arr2) << endl;
    return (*((double*) arr1) > *((double*) arr2)) - (*((double*) arr2) > *((double*) arr1));
}

double** splitRegions(function<double (double)> testFunc) {
    double** regions = new double*[regionsBudget];
    for (int i = 0; i < regionsBudget; i++) {
        regions[i] = new double[5];
    }

    int numRegions = 1;
    tuple<double, double, double> curr;
    priority_queue<tuple<double, double, double>> myHeap;
    //Each tuple has error estimate, startX, stopX
    myHeap.push(tuple<double, double, double>(estimateError(a, b, testFunc),a,b));
    double start, end, mid;
    while (numRegions < regionsBudget) {
        curr = myHeap.top();
        start = get<1>(curr);
        end = get<2>(curr);
        mid = (start + end) / 2;
        myHeap.pop();
        myHeap.push(tuple<double, double, double>(estimateError(start, mid, testFunc), start, mid));
        myHeap.push(tuple<double, double, double>(estimateError(mid, end, testFunc), mid, end));
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

double polyApprox(double x, double** regions) {
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

double polyApproxEst(double** regions) {
    double estimate = 0;
    for (int i = 0; i < regionsBudget; i++) {
        estimate += polyInt(regions[i][0], regions[i][1], regions[i][2], regions[i][3], regions[i][4]);
        if (isnan(estimate) && regionsBudget == 134) {
            cout << "Estimate is no longer a number!\n";
            cout << regions[i][0] << ", " << regions[i][1] << ", " << regions[i][2] << ", " << regions[i][3] << ", " << regions[i][4] << endl;
        }
        //estimate += (regions[i][2] * regions[i][1] * regions[i][1] + regions[i][3] * regions[i][1] + regions[i][4]) - (regions[i][2] * regions[i][0] * regions[i][0] + regions[i][3] * regions[i][0] + regions[i][4]);
    }
    return estimate;
}

double residual(double x, double** regions, function<double (double)> testFunc) {
    return testFunc(x) - polyApprox(x,regions);
}

//pure monte carlo estimate of the original integral using the polynomial as a control variate
double residualMCEst(double** regions, function<double (double)> testFunc, double numSamples, mt19937 gen) {
    double est = 0;
    uniform_real_distribution<> dist(0, 1);
    for (int i = 0; i < numSamples; i++) {
        est += residual(dist(gen),regions,testFunc);
    }
    est /= numSamples;
    est += polyApproxEst(regions);
    return est;
}

long double groundTruth = 0.812835882622424311126150;//From wolfram alpha

int main() {
    random_device r;
    mt19937 gen(r());
    uniform_real_distribution<> dist1(0, 1);
    uniform_real_distribution<> dist2(1, 70);

    cout.precision(20);
    double temp [3];
    getCoefs(0, 1, foo,temp);
    cout << temp[0] << ", " << temp[1] << ", " << temp[2] << endl;
    //cout << determinant(2, 3, 7, -5, 23, 10, 4, 5, -1);
    regionsBudget = 25;
    fooCoefs[0] = dist1(gen);
    fooCoefs[2] = dist1(gen);
    fooCoefs[4] = dist1(gen);
    fooCoefs[6] = dist1(gen);
    fooCoefs[1] = dist2(gen);
    fooCoefs[3] = dist2(gen);
    fooCoefs[5] = dist2(gen);
    fooCoefs[7] = dist2(gen);
    double** myRegions = splitRegions(foo);

    //cout << fooCoefs[0] << "*sin(" << fooCoefs[1] << "*x) + " << fooCoefs[2] << "*sin(" << fooCoefs[3] << "*x) + " << fooCoefs[4] << "*cos(" << fooCoefs[5] << "*x) + " << fooCoefs[6] << "*cos(" << fooCoefs[7] << "*x) + " << endl;
    for (int i = 0; i < regionsBudget; i++) {
        cout << "[";
        for (int j = 0; j < 4; j++) {
            cout << myRegions[i][j] << ", ";
        }
        cout << myRegions[i][4];
        cout << "]," << endl;
    }

    for (double j = 1.5; j < 7.5; j += 0.05) {
        if ((int) exp(j) == (int) exp(j-0.05)) {
            continue;
        }
        cout << (((int) exp(j)) * 2 + 1) << ",";
    }
    cout << endl;

    //cout << fooCoefs[0] << "*sin(" << fooCoefs[1] << "*x) + " << fooCoefs[2] << "*sin(" << fooCoefs[3] << "*x) + " << fooCoefs[4] << "*cos(" << fooCoefs[5] << "*x) + " << fooCoefs[6] << "*cos(" << fooCoefs[7] << "*x) + " << endl;
    //cout << fooIntegralGroundTruth2() << endl;

    int numRandomCoefs = 1000;
    
    for (double j = 1.5; j < 7.5; j += 0.05) {
        if ((int) exp(j) == (int) exp(j-0.05)) {
            continue;
        }
        regionsBudget = (int) exp(j);
        //Remove later:
        regionsBudget /= 2;
        //cout << regionsBudget << endl;
        double RMSE = 0;
        //cout << "j: " << j << endl;
        for (int i = 0; i < numRandomCoefs; i++) {
            fooCoefs[0] = dist1(gen);
            fooCoefs[2] = dist1(gen);
            fooCoefs[4] = dist1(gen);
            fooCoefs[6] = dist1(gen);
            fooCoefs[1] = dist2(gen);
            fooCoefs[3] = dist2(gen);
            fooCoefs[5] = dist2(gen);
            fooCoefs[7] = dist2(gen);
            
            myRegions = splitRegions(foo2);
            //RMSE += (polyApproxEst(myRegions)-fooIntegralGroundTruth2()) * (polyApproxEst(myRegions)-fooIntegralGroundTruth2());
            RMSE += pow(residualMCEst(myRegions,foo2,2*regionsBudget+1,gen)-fooIntegralGroundTruth2(),2);
            
            if (j < 1.7 && i < 20) {
                //cout << "hi\n";
                //cout << endl << fooCoefs[0] << "*sin(" << fooCoefs[1] << "*x) + " << fooCoefs[2] << "*sin(" << fooCoefs[3] << "*x) + " << fooCoefs[4] << "*cos(" << fooCoefs[5] << "*x) + " << fooCoefs[6] << "*cos(" << fooCoefs[7] << "*x) + " << endl;
                //cout << fooCoefs[0] << "," << fooCoefs[1] << "," << fooCoefs[2] << "," << fooCoefs[3] << "," << fooCoefs[4] << "," << fooCoefs[5] << "," << fooCoefs[6] << "," << fooCoefs[7] << endl;
                // << ((int) exp(j)) << ": " << polyApproxEst(myRegions) << " versus " << fooIntegralGroundTruth2() << endl;
            }
        }
        RMSE /= numRandomCoefs;
        if (isnan(RMSE)) {
            cout << endl <<  regionsBudget << endl;
        }
        RMSE = sqrt(RMSE);
        //cout << ((int) exp(j)) << ": " << polyApproxEst(myRegions) << endl;
        //cout << abs(polyApproxEst(myRegions)-groundTruth) << ",";
        cout << RMSE << ",";
    }
    cout << endl;
    //cout << polyApproxEst(myRegions) << endl;

    delete myRegions;
    myRegions = nullptr;
}