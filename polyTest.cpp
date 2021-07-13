#include <iostream>
#include <functional>
#include <random>
#include <time.h>
#include <math.h>
#include <fstream> 
#include <bits/stdc++.h>
#include <complex.h>
#include <algorithm>    //For random_shuffle
using namespace std;

/*
struct sampleSet1D {
    string samplingStrategy;//"Pure Monte Carlo", "Stratified", etc
    int numTrials;//Number of runs for each lambda
    int numLambdas;//Number of values for lambda tested
    double** xValues;//X values generated by the sampling strategy
    //xValues is a numTrials x numLambdas array
};*/

const int primes[10] = {2,3,5,7,11,13,17,19,23,29};

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

//Standard deviation (sigma) is 1
//Mean (mu) is lambda
//Returns the derivative of the gaussian w.r.t. its mean
double gaussianDerivativeWRTMean1D(double x, double lambda) {
    return exp(-0.5 * (x - lambda) * (x - lambda)) * (x - lambda) / sqrt(2 * M_PI);
}

//Standard deviation (sigma) is 1
//Mean (mu) is lambda
//Returns the integral from 0 to 1 of the derivative of the gaussian w.r.t. its mean
double groundTruthGaussianDerivativeWRTMean1D(double lambda) {
    return -(exp(-0.5 * (1 - lambda) * (1 - lambda)) - exp(-0.5 * lambda * lambda))/sqrt(2 * M_PI);
}

double disk(double x, double y) {
    return (double) (x*x + y*y < 2/M_PI);
}

double triangle(double x, double y) {
    return (double) (y > x);
}

double step2D(double x, double y) {
    return (double) (x < 1/M_PI);
}

double gaussian(double x, double y) {
    return exp(-x*x - y*y);
}

double bilinear(double x, double y) {
    return x*y;
}

double radicalInverse(int base, int a) {
    //Could be optimized based on https://www.pbr-book.org/3ed-2018/Sampling_and_Reconstruction/The_Halton_Sampler, but currently isn't
    int temp;
    double result = 0;
    double multBy = 1.0/base;
    while (a != 0) {
        temp = a % base;
        result += temp * multBy;
        //cout << (temp * multBy) << endl;
        multBy /= base;
        a /= base;
    }
    return result;
}

double* genHaltonSeq1D(int N, double a, double b, mt19937 & gen) {
    double* seqArr = new double[N];
    for (int i = 0; i < N; i++) {
        seqArr[i] = (b-a) * radicalInverse(2, i+1) + a;
    }
    return seqArr;
}

double** genHaltonSeq2D(int N, double a, double b, double c, double d, mt19937 & gen) {
    double** seqArr = new double*[N];
    for (int i = 0; i < N; i++) {
        seqArr[i] = new double[2];
        seqArr[i][0] = (b-a) * radicalInverse(primes[0], i+1) + a;
        seqArr[i][1] = (d-c) * radicalInverse(primes[1], i+1) + c;
    }
    return seqArr;
}

double** genHaltonSeqRot2D(int N, double a, double b, double c, double d, mt19937 & gen) {
    //uniform_real_distribution<> dist1(0, b-a);
    //uniform_real_distribution<> dist2(0, d-c);
    uniform_real_distribution<> dist3(0, 1);
    double** seqArr = new double*[N];
    double offset = dist3(gen);
    //double rot1 = dist1(gen);
    //double rot2 = dist2(gen);
    for (int i = 0; i < N; i++) {
        seqArr[i] = new double[2];
        seqArr[i][0] = (b-a) * (radicalInverse(primes[0], i+1) + offset) + a;
        seqArr[i][1] = (d-c) * (radicalInverse(primes[1], i+1) + offset) + c;
        if (seqArr[i][0] > b) {
            seqArr[i][0] -= b-a;
        }
        if (seqArr[i][1] > d) {
            seqArr[i][1] -= d-c;
        }
    }
    return seqArr;
}

double* genPureMonteCarlo1D(int N, double a, double b, mt19937 & gen) {
    uniform_real_distribution<> dist(a, b);
    double* seqArr = new double[N];
    for (int i = 0; i < N; i++) {
        seqArr[i] = dist(gen);
    }
    return seqArr;
}

double** genPureMonteCarlo2D(int N, double a, double b, double c, double d, mt19937 & gen) {
    uniform_real_distribution<> distX(a, b);
    uniform_real_distribution<> distY(c, d);
    double** seqArr = new double*[N];
    for (int i = 0; i < N; i++) {
        seqArr[i] = new double[2];
        seqArr[i][0] = distX(gen);
        seqArr[i][1] = distY(gen);
    }
    return seqArr;
}

double* genAntitheticMonteCarlo1D(int N, double a, double b, mt19937 & gen) {
    uniform_real_distribution<> dist(a, b);
    double* seqArr = new double[N];
    for (int i = 0; i < N/2.0; i++) {
        seqArr[i] = dist(gen);
        seqArr[N-i-1] = b-(seqArr[i]-a);
    }
    return seqArr;
}

double** genAntitheticMonteCarlo2D(int N, double a, double b, double c, double d, mt19937 & gen) {
    uniform_real_distribution<> distX(a, b);
    uniform_real_distribution<> distY(c, d);
    double** seqArr = new double*[N];
    for (int i = 0; i < N; i++)
    	seqArr[i] = new double[2];
    for (int i = 0; i < N/2.0; i++) {
        seqArr[i] = new double[2];
        seqArr[i][0] = distX(gen);
        seqArr[i][1] = distY(gen);
        seqArr[N-i-1][0] = b-(seqArr[i][0]-a);
        seqArr[N-i-1][1] = d-(seqArr[i][1]-c);
    }
    return seqArr;
}

double* genUniform1D(int N, double a, double b, mt19937 & gen) {
    double strataSize = (b-a)/N;
    double* seqArr = new double[N];
    for (int i = 0; i < N; i++) {
        seqArr[i] = a + strataSize/2 + i * strataSize;
    }
    return seqArr;
}

double** genUniform2D(int N, double a, double b, double c, double d, mt19937 & gen) {
    int dimN = sqrt(N);
    double strataSizeX = (b-a)/dimN;
    double strataSizeY = (d-c)/dimN;
    double** seqArr = new double*[dimN * dimN];
    for (int i = 0; i < dimN; i++) {
        for (int j = 0; j < dimN; j++) {
            seqArr[i * dimN + j] = new double[dimN * dimN];
            seqArr[i * dimN + j][0] = a + (i + 0.5) * strataSizeX;
            seqArr[i * dimN + j][1] = c + (j + 0.5) * strataSizeY;
        }
    }
    return seqArr;
}

double* genStratified1D(int N, double a, double b, mt19937 & gen) {
    uniform_real_distribution<> dist(0, 1);
    double strataSize = (b-a)/N;
    double* seqArr = new double[N];
    for (int i = 0; i < N; i++) {
        seqArr[i] = a + strataSize * (i + dist(gen));
    }
    return seqArr;
}

double** genStratified2D(int N, double a, double b, double c, double d, mt19937 & gen) {
    uniform_real_distribution<> dist(0, 1);
    int dimN = sqrt(N);
    double strataSizeX = (b-a)/dimN;
    double strataSizeY = (d-c)/dimN;
    double** seqArr = new double*[dimN * dimN];
    for (int i = 0; i < dimN; i++) {
        for (int j = 0; j < dimN; j++) {
            seqArr[i * dimN + j] = new double[2];
            seqArr[i * dimN + j][0] = a + strataSizeX * (i + dist(gen));
            seqArr[i * dimN + j][1] = c + strataSizeY * (j + dist(gen));
        }
    }
    return seqArr;
}

double** genNRooks2D(int N, double a, double b, double c, double d, mt19937 & gen) {
    //somewhat based on pseudo code from https://cs.dartmouth.edu/~wjarosz/publications/subr16fourier.html
    uniform_real_distribution<> dist(0,1);
    double strataSizeX = (b-a)/N;
    double strataSizeY = (d-c)/N;

    vector<double> sampleXs(N,0);
    vector<double> sampleYs(N,0);
    for (int i = 0; i < N; i++) {
        sampleXs[i] = a + (i + dist(gen)) * strataSizeX;
        sampleYs[i] = c + (i + dist(gen)) * strataSizeY;
    }
    random_shuffle(sampleXs.begin(),sampleXs.end());//Shuffle the sample x coordinates
    random_shuffle(sampleYs.begin(),sampleYs.end());//Shuffle the sample y coordinates

    double** seqArr = new double*[N];
    for (int i = 0; i < N; i++) {
        seqArr[i] = new double[2];
        seqArr[i][0] = sampleXs[i];
        seqArr[i][1] = sampleYs[i];
    }
    return seqArr;
}

double** genMultiJitter2D(int N, double a, double b, double c, double d, mt19937 & gen) {
    //somewhat based on pseudo code from https://cs.dartmouth.edu/~wjarosz/publications/subr16fourier.html
    N = sqrt(N);//Code written for N*N samples, so take square root to adjust to that
    uniform_real_distribution<> dist(0,1);
    double estimate = 0;
    double strataSizeX = (b-a)/(N*N);
    double strataSizeY = (d-c)/(N*N);
    //The outer vector contains the columns
    //Each inner vector contains the x values of each of the points in the columns
    vector<vector<double>> sampleXs(N,vector<double>(N,0));
    //The outer vector contains the rows
    //Each inner vector contains the y values of each of the points in the rows
    vector<vector<double>> sampleYs(N,vector<double>(N,0));

    for (int i = 0; i < N; i++) {
        //column i
        //Also simultaneously row i
        for (int j = 0; j < N; j++) {
            //point j in column i
            //Also simultaneously point j in row i
            sampleXs[i][j] = a + (i*N + j + dist(gen)) * strataSizeX;
            sampleYs[i][j] = c + (i*N + j + dist(gen)) * strataSizeY;
        }
    }

    for (int i = 0; i < N; i++) {
        random_shuffle(sampleXs[i].begin(),sampleXs[i].end());//Shuffle the sample x coordinates
        random_shuffle(sampleYs[i].begin(),sampleYs[i].end());//Shuffle the sample y coordinates
    }

    double** seqArr = new double*[N*N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            seqArr[i * N + j] = new double[2];
            seqArr[i * N + j][0] = sampleXs[i][j];
            seqArr[i * N + j][1] = sampleYs[j][i];
        }
    }

    return seqArr;
}

double estimateIntegral1D(double a, double b, int N, double lambda, function<double (double,double)> F, function<double* (int, double, double, mt19937 &)> sampleGen, mt19937 & gen) {
    double* samplePoints = sampleGen(N, a, b, gen);
    double estimate = 0;
    double temp;
    for (int i = 0; i < N; i++) {
        temp = F(samplePoints[i],lambda);
        estimate += temp;
    }
    estimate /= N;
    estimate *= (b-a);
    delete samplePoints;
    return estimate;
}

double estimateIntegral2D(double a, double b, double c, double d, int N, function<double (double,double)> F, function<double** (int, double, double, double, double, mt19937 &)> sampleGen, mt19937 & gen) {
    double** samplePoints = sampleGen(N, a, b, c, d, gen);
    double estimate = 0;
    double temp;
    for (int i = 0; i < N; i++) {
        temp = F(samplePoints[i][0],samplePoints[i][1]);
        estimate += temp;
    }
    estimate /= N;
    estimate *= (b-a);
    for (int i = 0; i < N; i++) {
        delete samplePoints[i];
    }
    delete samplePoints;
    return estimate;
}

double* fourierCoefs1D(double a, double b, int N, int numTrials, double wStep, int numW, function<double* (int, double, double, mt19937 &)> sampleGen, mt19937 & gen) {
    double* samples;

    double* coefs = new double[numW];
    for (int i = 0; i < numW; i++) {
        coefs[i] = 0;
    }
    complex<double> temp(0,0);
    double tempAngle;
    for (int i = 0; i < numTrials; i++) {
        samples = sampleGen(N, a, b, gen);
        for (int w = 0; w < numW; w++) {
            temp = (0,0);
            for (int j = 0; j < N; j++) {
                temp += exp(-2 * M_PI * wStep * w * samples[j] * complex<double>(0,1));
            }
            coefs[w] += norm(temp) / (N * N * numTrials);
        }
    }
    return coefs;
}

double** fourierCoefs2D(double a, double b, double c, double d, int N, int numTrials, double wStep, int numW, function<double** (int, double, double, double, double, mt19937 &)> sampleGen, mt19937 & gen) {
    //Round N to close perfect square
    //N = sqrt(N);
    //N = N * N;
    
    double** samples;
    double** spectra = new double*[numW*2+1];
    for (int i = 0; i < numW*2+1; i++) {
        spectra[i] = new double[numW*2+1];
        for (int j = 0; j < numW*2+1; j++) {
            spectra[i][j] = 0;
        }
    }
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
        samples = sampleGen(N, a, b, c, d, gen);
        for (int wX = -numW; wX <= numW; wX++) {
            for (int wY = -numW; wY <= numW; wY++) {
                temp = (0,0);
                for (int j = 0; j < N; j++) {
                    temp += exp(-2 * M_PI * wStep * (wX * samples[j][0] + wY * samples[j][1]) * complex<double>(0,1));
                }
                spectra[-wY + numW][wX + numW] += norm(temp) / (N * numTrials);//Each coefficient multiplied by sqrt(N) compared to regular (?) formula
            }
        }
        for (int i = 0; i < N; i++) {
            delete samples[i];
        }
        delete samples;
    }
    return spectra;
}

/**
 * Creates intensity strip images based on estimating the integral of a particular 1D function over the interval [0, 1]
 * testFunc is the function to approximate the integral of, and groundTruthFunc returns the correct value of the integral for a particular lambda
 * Lambda will range from 0 to 1, with numLambdas intermediate values tested, each corresponding to a particular column of the strips
 * numTrials trials will be done for each value of lambda, corresponding to the rows of the strips
 * A strip will be generated for ground truth, pure Monte Carlo (uniformly random sampling), uniform (sample points at the center of each strata),
 * and stratified/jittered (sample points at a random location within each strata)
 */
void makeIntensityStrips(int numLambdas, int numSamples, int numTrials, int imgWidth, int imgHeight, function<double (double,double)> testFunc, function<double (double)> groundTruthFunc, mt19937 & gen, string fileName) {
    double* samples;
    unsigned char image[imgHeight*imgWidth];
    for (int i = 0; i < imgHeight*imgWidth; i++) {
        image[i] = 255;
    }
    double temp;
    double minVal = 10000000000;
    double maxVal = -10000000000;
    for (double j = 0; j < numLambdas; j++) {
        temp = groundTruthFunc(j/numLambdas);
        if (temp < minVal) {
            minVal = temp;
        }
        if (temp > maxVal) {
            maxVal = temp;
        }
    }
    //Ground Truth
    for (int i = 0; i < 100; i++) {
        for (double j = 0; j < numLambdas; j++) {
            image[(unsigned int) (i*imgWidth + j + 20)] = (unsigned char) max(0.0, min(255.0, 30 + (180 / (maxVal - minVal)) * (groundTruthFunc(j/numLambdas) - minVal)));
            //cout << groundTruthFunc(j/numLambdas) << endl;
        }
    }
    cout << "ground truth done" << endl;
    //Monte Carlo
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genPureMonteCarlo1D, gen);
            image[(unsigned int) ((i + 110)*imgWidth + j + 20)] = (unsigned char) max(0.0, min(255.0, 30 + (180 / (maxVal - minVal)) * (temp - minVal)));
        }
    }
    cout << "pure monte done" << endl;
    //Uniform
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genUniform1D, gen);
            image[(unsigned int) ((i + 220)*imgWidth + j + 20)] = (unsigned char) max(0.0, min(255.0, 30 + (180 / (maxVal - minVal)) * (temp - minVal)));
        }
    }
    //Stratified
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genStratified1D, gen);
            image[(unsigned int) ((i + 330)*imgWidth + j + 20)] = (unsigned char) max(0.0, min(255.0, 30 + (180 / (maxVal - minVal)) * (temp - minVal)));
        }
    }
    //Halton (Since it's 1D this is a.k.a. van der Corput)
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genHaltonSeq1D, gen);
            image[(unsigned int) ((i + 440)*imgWidth + j + 20)] = (unsigned char) max(0.0, min(255.0, 30 + (180 / (maxVal - minVal)) * (temp - minVal)));
        }
    }
    //Antithetic Monte Carlo
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genAntitheticMonteCarlo1D, gen);
            image[(unsigned int) ((i + 550)*imgWidth + j + 20)] = (unsigned char) max(0.0, min(255.0, 30 + (180 / (maxVal - minVal)) * (temp - minVal)));
        }
    }
    //Make image file
    ofstream ofs(fileName, ios::out | ios::binary);
    ofs << "P6\n" << imgWidth << " " << imgHeight << "\n255\n"; 
    for (int i = 0; i < imgWidth * imgHeight; ++i) { 
        ofs << image[i] << image[i] << image[i];
    } 
    ofs.close();
}

/**
 * Prints the RMSE from estimating the integral of a particular 1D function over the interval [0, 1]
 * testFunc is the function to approximate the integral of, and groundTruthFunc returns the correct value of the integral for a particular lambda
 * Lambda will range from 0 to 1, with numLambdas intermediate values tested
 * numTrials trials will be done for each value of lambda
 * RMSE will be calculated for pure Monte Carlo (uniformly random sampling), uniform (sample points at the center of each strata),
 * and stratified/jittered (sample points at a random location within each strata)
 */
void printRMSE1D(int numLambdas, int numSamples, int numTrials, function<double (double,double)> testFunc, function<double (double)> groundTruthFunc, mt19937 & gen) {
    double avgError = 0;
    double temp;
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genPureMonteCarlo1D, gen);
            avgError += (temp-groundTruthFunc(j/numLambdas)) * (temp-groundTruthFunc(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials*numLambdas)) << endl << endl;
    avgError = 0;
    cout << "Antithetic Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genAntitheticMonteCarlo1D, gen);
            avgError += (temp-groundTruthFunc(j/numLambdas)) * (temp-groundTruthFunc(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials*numLambdas)) << endl << endl;
    avgError = 0;
    cout << "Uniform:\n";
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genUniform1D, gen);
            avgError += (temp-groundTruthFunc(j/numLambdas)) * (temp-groundTruthFunc(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials*numLambdas)) << endl << endl;
    avgError = 0;
    cout << "Stratified:\n";
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genStratified1D, gen);
            avgError += (temp-groundTruthFunc(j/numLambdas)) * (temp-groundTruthFunc(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials*numLambdas)) << endl << endl;
    avgError = 0;
    cout << "Halton:\n";
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genHaltonSeq1D, gen);
            avgError += (temp-groundTruthFunc(j/numLambdas)) * (temp-groundTruthFunc(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials*numLambdas)) << endl << endl;
}

/**
 * Prints the variance from estimating the integral of a particular 1D function over the interval [0, 1]
 * testFunc is the function to approximate the integral of, and groundTruthFunc returns the correct value of the integral for a particular lambda
 * Lambda will range from 0 to 1, with numLambdas intermediate values tested
 * numTrials trials will be done for each value of lambda
 * Variance will be calculated for pure Monte Carlo (uniformly random sampling),
 * antithetic Monte Carlo (points generated in pairs as x and 1-x, where x is generated uniformly randomly)
 * uniform (sample points at the center of each strata),
 * and stratified/jittered (sample points at a random location within each strata)
 */
void printVariance1D(int numLambdas, int numSamples, int numTrials, function<double (double,double)> testFunc, function<double (double)> groundTruthFunc, mt19937 & gen) {
    double avgVar = 0;
    double tempVar = 0;
    double avgEst = 0;
    double ests[numTrials];
    double temp;
    cout << "Pure Monte Carlo:\n";
    for (double j = 0; j < numLambdas; j++) {
        tempVar = 0;
        avgEst = 0;
        for (int i = 0; i < numTrials; i++) {
            ests[i] = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genPureMonteCarlo1D, gen);
            avgEst += ests[i];
        }
        avgEst /= numTrials;
        for (int i = 0; i < numTrials; i++) {
            tempVar += (ests[i]-avgEst) * (ests[i]-avgEst);
        }
        avgVar += tempVar/(numTrials - 1);
    }
    cout << "Variance: " << avgVar/numLambdas << endl << endl;
    avgVar = 0;
    cout << "Antithetic Monte Carlo:\n";
    for (double j = 0; j < numLambdas; j++) {
        tempVar = 0;
        avgEst = 0;
        for (int i = 0; i < numTrials; i++) {
            ests[i] = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genAntitheticMonteCarlo1D, gen);
            avgEst += ests[i];
        }
        avgEst /= numTrials;
        for (int i = 0; i < numTrials; i++) {
            tempVar += (ests[i]-avgEst) * (ests[i]-avgEst);
        }
        avgVar += tempVar/(numTrials - 1);
    }
    cout << "Variance: " << avgVar/numLambdas << endl << endl;
    avgVar = 0;
    cout << "Uniform:\n";
    for (double j = 0; j < numLambdas; j++) {
        tempVar = 0;
        avgEst = 0;
        for (int i = 0; i < numTrials; i++) {
            ests[i] = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genUniform1D, gen);
            avgEst += ests[i];
        }
        avgEst /= numTrials;
        for (int i = 0; i < numTrials; i++) {
            tempVar += (ests[i]-avgEst) * (ests[i]-avgEst);
        }
        avgVar += tempVar/(numTrials - 1);
    }
    cout << "Variance: " << avgVar/numLambdas << endl << endl;
    avgVar = 0;
    cout << "Stratified:\n";
    for (double j = 0; j < numLambdas; j++) {
        tempVar = 0;
        avgEst = 0;
        for (int i = 0; i < numTrials; i++) {
            ests[i] = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genStratified1D, gen);
            avgEst += ests[i];
        }
        avgEst /= numTrials;
        for (int i = 0; i < numTrials; i++) {
            tempVar += (ests[i]-avgEst) * (ests[i]-avgEst);
        }
        avgVar += tempVar/(numTrials - 1);
    }
    cout << "Variance: " << avgVar/numLambdas << endl << endl;
    avgVar = 0;
    cout << "Halton:\n";
    for (double j = 0; j < numLambdas; j++) {
        tempVar = 0;
        avgEst = 0;
        for (int i = 0; i < numTrials; i++) {
            ests[i] = estimateIntegral1D(0, 1, numSamples, j/numLambdas, testFunc, genHaltonSeq1D, gen);
            avgEst += ests[i];
        }
        avgEst /= numTrials;
        for (int i = 0; i < numTrials; i++) {
            tempVar += (ests[i]-avgEst) * (ests[i]-avgEst);
        }
        avgVar += tempVar/(numTrials - 1);
    }
    cout << "Variance: " << avgVar/numLambdas << endl << endl;
}

/**
 * Prints the RMSE from estimating the integral of a particular 2D function over the interval x in [0, 1], y in [0, 1]
 * testFunc is the function to approximate the integral of, and groundTruth is the correct value for the integral
 * numTrials trials will be done for each value of lambda
 * RMSE will be calculated for pure Monte Carlo (uniformly random sampling), uniform (sample points at the center of each strata),
 * and stratified/jittered (sample points at a random location within each strata)
 */
void printRMSE2D(int numSamples, int numTrials, function<double (double,double)> testFunc, double groundTruth, mt19937 & gen) {
    double avgError = 0;
    double temp;
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genPureMonteCarlo2D, gen);
        avgError += (temp-groundTruth) * (temp-groundTruth);
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials)) << endl << endl;
    avgError = 0;
    cout << "Uniform:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genUniform2D, gen);
        avgError += (temp-groundTruth) * (temp-groundTruth);
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials)) << endl << endl;
    avgError = 0;
    cout << "Stratified:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genStratified2D, gen);
        avgError += (temp-groundTruth) * (temp-groundTruth);
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials)) << endl << endl;
    avgError = 0;
    cout << "N-Rooks:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genNRooks2D, gen);
        avgError += (temp-groundTruth) * (temp-groundTruth);
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials)) << endl << endl;
    avgError = 0;
    cout << "Multi-Jitter:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genMultiJitter2D, gen);
        avgError += (temp-groundTruth) * (temp-groundTruth);
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials)) << endl << endl;
    avgError = 0;
    cout << "Halton:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genHaltonSeq2D, gen);
        avgError += (temp-groundTruth) * (temp-groundTruth);
    }
    cout << "RMSE: " << sqrt(avgError/(numTrials)) << endl << endl;
}

/**
 * Prints the average error (not RMSE) from estimating the integral of a particular 2D function over the interval x in [0, 1], y in [0, 1]
 * testFunc is the function to approximate the integral of, and groundTruth is the correct value for the integral
 * numTrials trials will be done for each value of lambda
 * Error will be calculated for pure Monte Carlo (uniformly random sampling), uniform (sample points at the center of each strata),
 * and stratified/jittered (sample points at a random location within each strata)
 */
void printError2D(int numSamples, int numTrials, function<double (double,double)> testFunc, double groundTruth, mt19937 & gen) {
    //Mostly to recreate https://www.semanticscholar.org/paper/Progressive-Multi-Jittered-Sample-Sequences-%3A-Christensen-Kensler/c19a94e814c8bd56a12ee108daacd9c243c2b7ef/figure/1
    double avgError = 0;
    double temp;
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genPureMonteCarlo2D, gen);
        avgError += abs(temp-groundTruth);
    }
    cout << "Average Error: " << avgError/(numTrials) << endl << endl;
    avgError = 0;
    cout << "Uniform:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genUniform2D, gen);
        avgError += abs(temp-groundTruth);
    }
    cout << "RMSE: " << avgError/(numTrials) << endl << endl;
    avgError = 0;
    cout << "Stratified:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genStratified2D, gen);
        avgError += abs(temp-groundTruth);
    }
    cout << "RMSE: " << avgError/(numTrials) << endl << endl;
    avgError = 0;
    cout << "N-Rooks:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genNRooks2D, gen);
        avgError += abs(temp-groundTruth);
    }
    cout << "RMSE: " << avgError/(numTrials) << endl << endl;
    avgError = 0;
    cout << "Multi-Jitter:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genMultiJitter2D, gen);
        avgError += abs(temp-groundTruth);
    }
    cout << "RMSE: " << avgError/(numTrials) << endl << endl;
    avgError = 0;
    cout << "Halton:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = estimateIntegral2D(0, 1, 0, 1, numSamples, testFunc, genHaltonSeqRot2D, gen);
        avgError += abs(temp-groundTruth);
    }
    cout << "RMSE: " << avgError/(numTrials) << endl << endl;
}

void makePowerSpectra(int numSamples, int numTrials, int imgWidth, int imgHeight, int maxW, function<double** (int, double, double, double, double, mt19937 &)> sampleGen, mt19937 & gen, string fileName) {
    unsigned char image[imgHeight*imgWidth];
    for (int i = 0; i < imgHeight*imgWidth; i++) {
        image[i] = 255;
    }

    double** testSpectra = fourierCoefs2D(0, 1, 0, 1, numSamples, numTrials, 1, maxW, sampleGen, gen);
    for (int i = 0; i < maxW * 2 + 1; i++) {
        for (int j = 0; j < maxW * 2 + 1; j++) {
            image[i * imgWidth + j] = (unsigned char) min(255, (int) (200 * testSpectra[i][j]));
        }
    }
    for (int i = 0; i < maxW * 2 + 1; i++) {
        delete testSpectra[i];
    }
    delete testSpectra;

    //Make image file
    ofstream ofs(fileName, ios::out | ios::binary);
    ofs << "P6\n" << imgWidth << " " << imgHeight << "\n255\n"; 
    for (int i = 0; i < imgWidth * imgHeight; ++i) { 
        ofs << image[i] << image[i] << image[i];
    } 
    ofs.close();
}

void printConvergenceRates1D(int startN, int endN, int numLambdas, int numTrials, function<double (double,double)> testFunc, function<double (double)> groundTruthFunc, mt19937 & gen) {
    double avgError = 0;
    double temp;
    ofstream ofs("conv.txt", ios::out);
    //ofs << "Number of samples,Pure Monte Carlo,Uniform,Stratified,Halton\n";
    for (int n = startN; n <= endN; n++) {
        ofs << n << ",";
        avgError = 0;
        //Monte Carlo
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < 1; j += 1.0/numLambdas) {
                temp = estimateIntegral1D(0, 1, n, j, testFunc, genPureMonteCarlo1D, gen);
                avgError += (temp-groundTruthFunc(j)) * (temp-groundTruthFunc(j));
            }
            //ofs << "trial done" << endl;
        }
        ofs << sqrt(avgError/(numTrials * numLambdas)) << ",";
        avgError = 0;
        //Antithetic Monte Carlo
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < 1; j += 1.0/numLambdas) {
                temp = estimateIntegral1D(0, 1, n, j, testFunc, genAntitheticMonteCarlo1D, gen);
                avgError += (temp-groundTruthFunc(j)) * (temp-groundTruthFunc(j));
            }
            //ofs << "trial done" << endl;
        }
        ofs << sqrt(avgError/(numTrials * numLambdas)) << ",";
        avgError = 0;
        //Uniform
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < 1; j += 1.0/numLambdas) {
                temp = estimateIntegral1D(0, 1, n, j, testFunc, genUniform1D, gen);
                avgError += (temp-groundTruthFunc(j)) * (temp-groundTruthFunc(j));
            }
        }
        ofs << sqrt(avgError/(numTrials * numLambdas)) << ",";
        avgError = 0;
        //Stratified
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < 1; j += 1.0/numLambdas) {
                temp = estimateIntegral1D(0, 1, n, j, testFunc, genStratified1D, gen);
                avgError += (temp-groundTruthFunc(j)) * (temp-groundTruthFunc(j));
            }
        }
        ofs << sqrt(avgError/(numTrials * numLambdas)) << ",";
        avgError = 0;
        //Halton
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < 1; j += 1.0/numLambdas) {
                temp = estimateIntegral1D(0, 1, n, j, testFunc, genHaltonSeq1D, gen);
                avgError += (temp-groundTruthFunc(j)) * (temp-groundTruthFunc(j));
            }
        }
        ofs << sqrt(avgError/(numTrials * numLambdas)) << endl;
        cout << "n = " << n << " done" << endl;
    }
}

//Tests perfect square between startN^2 and endN^2, inclusive
void printConvergenceRates2D(int startN, int endN, int numTrials, function<double (double,double)> testFunc, double groundTruth, mt19937 & gen) {
    double avgError = 0;
    double temp;
    int n;
    ofstream ofs("conv.txt", ios::out);
    //cout << "Number of samples,Pure Monte Carlo,Uniform,Stratified,Halton\n";
    for (int j = startN; j <= endN; j++) {
        n = j * j;
        ofs << n << ",";
        avgError = 0;
        //Monte Carlo
        for (int i = 0; i < numTrials; i++) {
            temp = estimateIntegral2D(0, 1, 0, 1, n, testFunc, genPureMonteCarlo2D, gen);
            avgError += abs(temp-groundTruth);
        }
        ofs << avgError/(numTrials) << ",";
        avgError = 0;
        //Monte Carlo w/ antithetic sampling
        for (int i = 0; i < numTrials; i++) {
            temp = estimateIntegral2D(0, 1, 0, 1, n, testFunc, genAntitheticMonteCarlo2D, gen);
            avgError += abs(temp-groundTruth);
        }
        ofs << avgError/(numTrials) << ",";
        avgError = 0;
        //Uniform
        for (int i = 0; i < numTrials; i++) {
            temp = estimateIntegral2D(0, 1, 0, 1, n, testFunc, genUniform2D, gen);
            avgError += abs(temp-groundTruth);
        }
        ofs << avgError/(numTrials) << ",";
        avgError = 0;
        //Stratified
        for (int i = 0; i < numTrials; i++) {
            temp = estimateIntegral2D(0, 1, 0, 1, n, testFunc, genStratified2D, gen);
            avgError += abs(temp-groundTruth);
        }
        ofs << avgError/(numTrials) << ",";
        avgError = 0;
        //N-Rooks
        for (int i = 0; i < numTrials; i++) {
            temp = estimateIntegral2D(0, 1, 0, 1, n, testFunc, genNRooks2D, gen);
            avgError += abs(temp-groundTruth);
        }
        ofs << avgError/(numTrials) << ",";
        avgError = 0;
        //Multi-jitter
        for (int i = 0; i < numTrials; i++) {
            temp = estimateIntegral2D(0, 1, 0, 1, n, testFunc, genMultiJitter2D, gen);
            avgError += abs(temp-groundTruth);
        }
        ofs << avgError/(numTrials) << ",";
        avgError = 0;
        //Halton
        for (int i = 0; i < 1; i++) {
            temp = estimateIntegral2D(0, 1, 0, 1, n, testFunc, genHaltonSeqRot2D, gen);
            cout << "error for n = " << n << " is " << abs(temp-groundTruth) << endl;
            avgError += abs(temp-groundTruth);
        }
        ofs << avgError << endl;
    }
}

void printPoints2D(int N, function<double** (int, double, double, double, double, mt19937 &)> sampleGen, mt19937 & gen) {
    ofstream ofs("points.txt", ios::out);
    //cout << "Number of samples,Pure Monte Carlo,Uniform,Stratified,Halton\n";
    double** points = sampleGen(N, 0, 1, 0, 1, gen);
    for (int i = 0; i < N; i++) {
        ofs << points[i][0] << "," << points[i][1] << endl;
    }
    delete points;
}

int main(int argc, char** argv) {
    string fileName = argv[1];
    random_device r;
    mt19937 gen(r());
    int imgWidth = 500;
    int imgHeight = 700;
    int numSamples = 16;
    int numTrials = 100;
    int numLambdas = 250;
    double avgError = 0;
    double avgError2 = 0;
    double temp;

    const double groundTruthDisk = 0.5;
    const double groundTruthTriangle = 0.5;
    const double groundTruthStep2D = 1/M_PI;
    const double groundTruthGaussian = 0.55774628535;
    const double groundTruthBilinear = 0.25;

    makeIntensityStrips(numLambdas,numSamples,numTrials,imgWidth,imgHeight,gaussianDerivativeWRTMean1D,groundTruthGaussianDerivativeWRTMean1D,gen,fileName);
    //printRMSE1D(numLambdas,numSamples,numTrials,gaussianDerivativeWRTMean1D,groundTruthGaussianDerivativeWRTMean1D,gen);
    //printVariance1D(numLambdas,numSamples,numTrials,gaussianDerivativeWRTMean1D,groundTruthGaussianDerivativeWRTMean1D,gen);
    //printError2D(256,1000,gaussian,groundTruthGaussian,gen);
    //makePowerSpectra(numSamples,numTrials,imgWidth,imgHeight,60,genHaltonSeq2D,gen,"test2.ppm");
    //radicalInverse(3,7);
    //printConvergenceRates1D(6,150,numLambdas,numTrials,gaussianDerivativeWRTMean1D,groundTruthGaussianDerivativeWRTMean1D,gen);
    //printConvergenceRates2D(2,40,numTrials,gaussian,groundTruthGaussian,gen);
    //printPoints2D(500,genPureMonteCarlo2D,gen);
}
