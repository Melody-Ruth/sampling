#include <iostream>
#include <functional>
#include <random>
#include <time.h>
#include <math.h>
#include <utility>
#include <queue>

using namespace std;

double determinant(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    return a * (e*i - h*f) - b * (d*i - g*f) + c * (d*h - g*e);
}

double foo(double x) {
    return 0.2 * sin(100*x) + 0.1 * sin(40*x) - 2 * (x - 0.6) * (x-0.6) + 0.15 * cos(30*x) + 1;
}

double pdf(double x) {
    return -3 * (x - 0.5) * (x - 0.5) + 1.25;
}

double pdf2(double x) {
    return -1.786 * x * x + 2.099 * x + 0.546;
}

double cdf(double x) {
    return -1/8.0 - (x - 0.5) * (x - 0.5) * (x - 0.5) + 1.25 * x;
}

double cdf2(double x) {
    return -0.595 * x * x * x + 1.05 * x * x + 0.546 * x;
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

double foo2(double x) {
    return foo(invCdf(x)) / pdf(invCdf(x));
}

double foo22(double x) {
    return foo(invCdf2(x)) / pdf2(invCdf2(x));
}

const int regionsBudget = 7;
const double epsilon = 0.00001;

double miniSimpsons(double intervalStart, double intervalEnd, function<double (double)> testFunc) {
    return (testFunc(intervalStart) + 4 * testFunc((intervalStart + intervalEnd) / 2.0) + testFunc(intervalEnd)) * (intervalEnd - intervalStart) / 6.0;
}

double miniTrapezoidal(double intervalStart, double intervalEnd, function<double (double)> testFunc) {
    return (intervalEnd - intervalStart) * (testFunc(intervalEnd) + testFunc(intervalStart)) / 2.0;
}

double estimateError(double intervalStart, double intervalEnd, function<double (double)> testFunc) {
    return abs(miniSimpsons(intervalStart, intervalEnd, testFunc) - miniTrapezoidal(intervalStart, intervalEnd, testFunc)) + epsilon * (intervalEnd - intervalStart);
}

void getCoefs(double x1, double x3, function<double (double)> testFunc) {
    //quadratic approximation
    double x2 = (x1 + x3) / 2.0;
    double y1 = testFunc(x1);
    double y2 = testFunc(x2);
    double y3 = testFunc(x3);

    double D = determinant(x1*x1, x1, 1, x2*x2, x2, 1, x3*x3, x3, 1);
    double Da = determinant(y1, x1, 1, y2, x2, 1, y3, x3, 1);
    double Db = determinant(x1*x1, y1, 1, x2*x2, y2, 1, x3*x3, y3, 1);
    double Dc = determinant(x1*x1, x1, y1, x2*x2, x2, y2, x3*x3, x3, y3);
    
    double a = Da / D;
    double b = Db / D;
    double c = Dc / D;

    cout << x1 << " to " << x3 << endl;
    cout << a << ", " << b << ", " << c << endl;
}

double a = 0;
double b = 1;

void splitRegions(function<double (double)> testFunc) {
    double regions[regionsBudget][2];

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
        getCoefs(regions[i][0], regions[i][1], testFunc);
        //cout << regions[i][0] << ", " << regions[i][1] << endl;
    }
}

int main() {
    //getCoefs(0, 1, foo);
    //cout << determinant(2, 3, 7, -5, 23, 10, 4, 5, -1);
    splitRegions(foo22);
}