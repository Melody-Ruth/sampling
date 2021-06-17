#include <iostream>
#include <functional>
#include <random>
#include <time.h>
#include <math.h>
#include <fstream> 
#include <bits/stdc++.h>
#include <complex.h>
using namespace std;

double sampleQuad(double x, double lambda) {
    return 600 * (x-lambda) * (x-lambda) + 1;
}

double groundTruthQuad(double lambda) {
    //just for (x-lambda)^2 + 2 from 0 to 1
    return 201 - 600 * lambda + 600 * lambda * lambda;
}

double sampleStep(double x, double lambda) {
    return (double) 150 * (x > lambda);
}

double groundTruthStep(double lambda) {
    return 150 * (1-lambda);
}

double disk(double x, double y, double lambda) {
    //lambda isn't used
    return (double) (x*x + y*y < 2/M_PI);
}

double groundTruthDisk(double lambda) {
    //lambda isn't used
    return 0.5;
}

double triangle(double x, double y, double lambda) {
    //lambda isn't used
    return (double) (y > x);
}

double groundTruthTriangle(double lambda) {
    //lambda isn't used
    return 0.5;
}

double step2D(double x, double y, double lambda) {
    //lambda isn't used
    return (double) (x < 1/M_PI);
}

double groundTruthStep2D(double lambda) {
    //lambda isn't used
    return 1/M_PI;
}

double gaussian(double x, double y, double lambda) {
    //lambda isn't used
    return exp(-x*x - y*y);
}

double groundTruthGaussian(double lambda) {
    //lambda isn't used
    return 0.55774628535;
}

double bilinear(double x, double y, double lambda) {
    //lambda isn't used
    return x*y;
}

double groundTruthBilinear(double lambda) {
    //lambda isn't used
    return 0.25;
}

double pureMonteCarlo(double a, double b, int N, double lambda, function<double (double,double)> F, mt19937 & gen) {
    double estimate = 0;
    uniform_real_distribution<> dist(a, b);
    double temp;
    for (int i = 0; i < N; i++) {
        temp = F(dist(gen),lambda);
        estimate += temp;
    }
    estimate /= N;
    estimate *= (b-a);
    return estimate;
}

/**
 * Estimates the integral of F(x,y,lambda) over the interval a <= x <= b, c <= y <= d
 * Does this by placing N random points in the interval
 */
double pureMonteCarlo2D(double a, double b, double c, double d, int N, double lambda, function<double (double,double,double)> F, mt19937 & gen) {
    //cout << N << endl;
    double estimate = 0;
    uniform_real_distribution<> distX(a, b);
    uniform_real_distribution<> distY(c,d);
    double temp;
    double tempX;
    double tempY;
    for (int i = 0; i < N; i++) {
        tempX = distX(gen);
        tempY = distY(gen);
        temp = F(distX(gen),distY(gen),lambda);
        estimate += temp;
    }
    estimate /= N;
    estimate *= (b-a);
    estimate *= (d-c);
    return estimate;
}

double* pureMonteCarloFourierCoefs(double a, double b, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    double samples[N];
    uniform_real_distribution<> dist(a,b);
    uniform_real_distribution<> dist2(0,2*M_PI);

    double* coefs = new double[numW];
    for (int i = 0; i < numW; i++) {
        coefs[i] = 0;
    }
    complex<double> temp(0,0);
    double tempAngle;
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < N; j++) {
            samples[j] = dist(gen);
        }
        /*for (int j = 0; j < N; j++) {
            cout << samples[j] << " ";
        }
        cout << endl;*/
        for (int w = 0; w < numW; w++) {
            temp = (0,0);
            for (int j = 0; j < N; j++) {
                temp += exp(-2 * M_PI * wStep * w * samples[j] * complex<double>(0,1));
                //tempAngle = dist2(gen);
                //temp += (cos(tempAngle),sin(tempAngle));
                if (w == 2) {
                    //cout << samples[j] << " " << norm(exp(-2 * M_PI * wStep * w * samples[j] * complex<double>(0,1))) << endl;
                }
                //cout << temp << endl;
            }
            if (w == 1) {
                //cout << (norm(temp) / (N * N)) << endl;
            }
            coefs[w] += norm(temp) / (N * N * numTrials);
            //cout << "This coefficient: " << coefs[w] << endl;
        }
    }
    return coefs;
}

//at 0.5 in each strata
double uniform(double a, double b, int N, double lambda, function<double (double,double)> F) {
    double estimate = 0;
    double strataSize = (b-a)/N;
    for (int i = 0; i < N; i++) {
        estimate += F(a + strataSize/2 + i * strataSize,lambda);
    }
    estimate /= N;
    estimate *= (b-a);
    return estimate;
}

double* uniformFourierCoefs(double a, double b, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    double samples[N];
    uniform_real_distribution<> dist(0, 1);
    double strataSize = (b-a)/N;

    double* coefs = new double[numW];
    for (int i = 0; i < numW; i++) {
        coefs[i] = 0;
    }
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < N; j++) {
            samples[j] = a + (j + 0.5) * strataSize;
        }
        /*for (int j = 0; j < N; j++) {
            cout << samples[j] << " ";
        }
        cout << endl;*/
        for (int w = 0; w < numW; w++) {
            temp = (0,0);
            for (int j = 0; j < N; j++) {
                temp += exp(-2 * M_PI * wStep * w * samples[j] * complex<double>(0,1));
                //cout << exp(-2 * M_PI * wStep * w * samples[j] * complex<double>(0,1)) << endl;
            }
            //cout << temp << endl;
            coefs[w] += norm(temp) / (N * N * numTrials);
            //cout << "This coefficient: " << coefs[w] << endl;
        }
    }
    return coefs;
}

double uniform2D(double a, double b, double c, double d, int N, double lambda, function<double (double,double,double)> F) {
    //will have less than N points if N isn't a perfect square
    double estimate = 0;
    int numRows = (int) sqrt(N);
    double strataSizeX = (b-a)/numRows;
    double strataSizeY = (d-c)/numRows;
    double tempX;
    double tempY;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numRows; j++) {
            tempX = a + strataSizeX/2 + i * strataSizeX;
            tempY = c + strataSizeY/2 + j * strataSizeY;
            //cout << tempX << ", " << tempY << endl;
            estimate += F(tempX,tempY,lambda);
        }
    }
    estimate /= numRows * numRows;
    estimate *= (b-a);
    estimate *= (d-c);
    return estimate;
}

//random position in each strata
double stratified(double a, double b, int N, double lambda, function<double (double,double)> F, mt19937 & gen) {
    double estimate = 0;
    uniform_real_distribution<> dist(0,1);
    double strataSize = (b-a)/N;
    double temp;
    for (int i = 0; i < N; i++) {
        estimate += F(a + (i + dist(gen)) * strataSize, lambda);
    }
    estimate /= N;
    estimate *= (b-a);
    return estimate;
}

double* stratifiedFourierCoefs(double a, double b, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    double samples[N];
    uniform_real_distribution<> dist(0, 1);
    double strataSize = (b-a)/N;

    double* coefs = new double[numW];
    for (int i = 0; i < numW; i++) {
        coefs[i] = 0;
    }
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < N; j++) {
            samples[j] = a + (j + dist(gen)) * strataSize;
        }
        /*for (int j = 0; j < N; j++) {
            cout << samples[j] << " ";
        }
        cout << endl;*/
        for (int w = 0; w < numW; w++) {
            temp = (0,0);
            for (int j = 0; j < N; j++) {
                temp += exp(-2 * M_PI * wStep * w * samples[j] * complex<double>(0,1));
                //cout << temp << endl;
            }
            //cout << coefs[w] << endl;
            coefs[w] += norm(temp) / (N * N * numTrials);
            //cout << "This coefficient: " << coefs[w] << endl;
        }
    }
    return coefs;
}

double stratified2D(double a, double b, double c, double d, int N, double lambda, function<double (double,double,double)> F, mt19937 & gen) {
    //will have less than N points if N isn't a perfect square
    uniform_real_distribution<> dist(0,1);
    double estimate = 0;
    int numRows = (int) sqrt(N);
    double strataSizeX = (b-a)/numRows;
    double strataSizeY = (d-c)/numRows;
    double tempX;
    double tempY;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numRows; j++) {
            tempX = a + (i + dist(gen)) * strataSizeX;
            tempY = c + (j + dist(gen)) * strataSizeY;
            //cout << tempX << ", " << tempY << endl;
            estimate += F(tempX,tempY,lambda);
        }
    }
    estimate /= numRows * numRows;
    estimate *= (b-a);
    estimate *= (d-c);
    return estimate;
}

double* stratified2DFourierCoefs(double a, double b, double c, double d, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    double sampleXs[N];
    double sampleYs[N];
    uniform_real_distribution<> dist(0, 1);
    double strataSizeX = (b-a)/N;
    double strataSizeY = (d-c)/N;

    double* coefs = new double[numW];
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < N; j++) {
            sampleXs[j] = a + (j + dist(gen)) * strataSizeX;
            sampleYs[j] = c + (j + dist(gen)) * strataSizeY;
        }
        for (int w = 0; w < numW; w++) {
            temp = (0,0);
            for (int j = 0; j < N; j++) {
                /*for (double k = 0; k < 2 * M_PI; k += M_PI/6) {
                    temp += exp(-2 * M_PI * wStep * w * (cos(k) * sampleXs[j] + sin(k) * sampleYs[j]) * (1/12) * complex<double>(0,1));
                }*/
                temp += exp(-2 * M_PI * wStep * w * (sampleXs[j]/sqrt(2) + sampleYs[j]/sqrt(2)) * complex<double>(0,1));
            }
            //cout << temp << endl;
            coefs[w] += norm(temp) * norm(temp) / (N * N * numTrials);
        }
    }
    return coefs;
}

int main(int argc, char** argv) {
    string fileName = argv[1];
    random_device r;
    mt19937 gen(r());
    int imgWidth = 500;
    int imgHeight = 500;
    int numSamples = 15;
    int numTrials = 100;
    int numLambdas = 250;
    unsigned char image[imgHeight*imgWidth];
    for (int i = 0; i < imgHeight*imgWidth; i++) {
        image[i] = 255;
    }
    double avgError = 0;
    double avgError2 = 0;
    double temp;
    /*
    cout << "1D:" << endl;
    cout << "Quadratic:\n";
    //cout << "Ground truth:\n";
    for (int i = 0; i < 100; i++) {
        for (double j = 0; j < numLambdas; j++) {
            //image[(unsigned int) (i*imgWidth + j + 20)] = groundTruthQuad(j/numLambdas);
            //cout << groundTruth(j) << " ";
        }
        //cout << endl;
    }
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = pureMonteCarlo(0,1,numSamples,j/numLambdas,sampleQuad,gen);
            //image[(unsigned int) ((i + 110)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas);
            avgError2 += abs(temp-groundTruthQuad(j/numLambdas));
        }
    }
    cout << "Average error: " << (avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (avgError/(numTrials*numLambdas)) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = uniform(0,1,numSamples,j/numLambdas,sampleQuad);
            //image[(unsigned int) ((i + 220)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas);
            avgError2 += abs(temp-groundTruthQuad(j/numLambdas));
        }
    }
    cout << "Average error: " << (avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (avgError/(numTrials*numLambdas)) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = stratified(0,1,numSamples,j/numLambdas,sampleQuad,gen);
            //image[(unsigned int) ((i + 330)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas);
            avgError2 += abs(temp-groundTruthQuad(j/numLambdas));
        }
    }
    cout << "Average error: " << (avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (avgError/(numTrials*numLambdas)) << "%\n";
    cout << endl;

    cout << "Step:\n";
    for (int i = 0; i < 100; i++) {
        for (double j = 0; j < numLambdas; j++) {
            image[(unsigned int) (i*imgWidth + j + 20)] = groundTruthStep(j/numLambdas);
        }
    }
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = pureMonteCarlo(0,1,numSamples,j/numLambdas,sampleStep,gen);
            image[(unsigned int) ((i + 110)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthStep(j/numLambdas))/groundTruthStep(j/numLambdas);
            avgError2 += abs(temp-groundTruthStep(j/numLambdas));
        }
    }
    cout << "Average error: " << (avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (avgError/(numTrials*numLambdas)) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = uniform(0,1,numSamples,j/numLambdas,sampleStep);
            //cout << temp << " vs " << groundTruthStep(j/numLambdas) << endl;
            image[(unsigned int) ((i + 220)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthStep(j/numLambdas))/groundTruthStep(j/numLambdas);
            avgError2 += abs(temp-groundTruthStep(j/numLambdas));
        }
    }
    cout << "Average error: " << (avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (avgError/(numTrials*numLambdas)) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = stratified(0,1,numSamples,j/numLambdas,sampleStep,gen);
            image[(unsigned int) ((i + 330)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthStep(j/numLambdas))/groundTruthStep(j/numLambdas);
            avgError2 += abs(temp-groundTruthStep(j/numLambdas));
        }
    }
    cout << "Average error: " << (avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (avgError/(numTrials*numLambdas)) << "%\n";
    cout << endl;


    cout << "--------------------------" << endl;
    cout << endl;
    cout << "2D:" << endl;
    cout << "Disk:\n";
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,disk,gen);
        avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
        avgError2 += abs(temp-groundTruthDisk(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,disk);
        avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,disk,gen);
        avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
        avgError2 += abs(temp-groundTruthDisk(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << endl;

    cout << "Triangle:\n";
    avgError = 0;
    avgError2 = 0;
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,triangle,gen);
        avgError += abs(temp-groundTruthTriangle(0))/groundTruthTriangle(0);
        avgError2 += abs(temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,triangle);
        avgError += abs(temp-groundTruthTriangle(0))/groundTruthTriangle(0);
        avgError2 += abs(temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,triangle,gen);
        avgError += abs(temp-groundTruthTriangle(0))/groundTruthTriangle(0);
        avgError2 += abs(temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << endl;

    cout << "Step:\n";
    avgError = 0;
    avgError2 = 0;
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,step2D,gen);
        avgError += abs(temp-groundTruthStep2D(0))/groundTruthStep2D(0);
        avgError2 += abs(temp-groundTruthStep2D(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,step2D);
        avgError += abs(temp-groundTruthStep2D(0))/groundTruthStep2D(0);
        avgError2 += abs(temp-groundTruthStep2D(0));
    }
   cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,step2D,gen);
        avgError += abs(temp-groundTruthStep2D(0))/groundTruthStep2D(0);
        avgError2 += abs(temp-groundTruthStep2D(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << endl;

    cout << "Gaussian:\n";
    cout << "Monte Carlo:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,gaussian,gen);
        avgError += abs(temp-groundTruthGaussian(0))/groundTruthGaussian(0);
        avgError2 += abs(temp-groundTruthGaussian(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,gaussian);
        avgError += abs(temp-groundTruthGaussian(0))/groundTruthGaussian(0);
        avgError2 += abs(temp-groundTruthGaussian(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,gaussian,gen);
        avgError += abs(temp-groundTruthGaussian(0))/groundTruthGaussian(0);
        avgError2 += abs(temp-groundTruthGaussian(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << endl;

    cout << "Bilinear:\n";
    cout << "Monte Carlo:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,bilinear,gen);
        avgError += abs(temp-groundTruthBilinear(0))/groundTruthBilinear(0);
        avgError2 += abs(temp-groundTruthBilinear(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,bilinear);
        avgError += abs(temp-groundTruthBilinear(0))/groundTruthBilinear(0);
        avgError2 += abs(temp-groundTruthBilinear(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,bilinear,gen);
        avgError += abs(temp-groundTruthBilinear(0))/groundTruthBilinear(0);
        avgError2 += abs(temp-groundTruthBilinear(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (avgError/numTrials) << "%\n";
    cout << endl;*/

    //Fourier analysis (Not working yet)
    double* testCoefs = stratifiedFourierCoefs(0,1,numSamples,100000,1,40,gen);
    for (int i = 0; i < 40; i++) {
        cout << (i) <<  "," << testCoefs[i]*numSamples << "\n";
    }
    cout << endl;
    delete testCoefs;

    //Make image file
    ofstream ofs(fileName, ios::out | ios::binary);
    ofs << "P6\n" << imgWidth << " " << imgHeight << "\n255\n"; 
    for (int i = 0; i < imgWidth * imgHeight; ++i) { 
        ofs << image[i] << image[i] << image[i];
    } 
    ofs.close();
}