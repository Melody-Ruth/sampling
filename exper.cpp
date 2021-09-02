#include <iostream>
#include <functional>
#include <random>
#include <time.h>
#include <math.h>

double foo(double x) {
    return 0.2 * sin(100*x) + 0.1 * sin(40*x) - 2 * (x - 0.6) * (x-0.6) + 0.15 * cos(30*x) + 1;
}

double pdf(double x) {
    return -3 * (x - 0.5) * (x - 0.5) + 1.25;
}

double cdf(double x) {
    -1/8.0 - (x - 0.5) * (x - 0.5) * (x - 0.5) + 1.25 * x;
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

double foo2(double x) {
    return foo(invCdf(x)) / pdf(invCdf(x));
}