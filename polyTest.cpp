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

double** pureMonteCarlo2DFourierCoefs(double a, double b, double c, double d, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    double sampleXs[N];
    double sampleYs[N];
    uniform_real_distribution<> distX(a,b);
    uniform_real_distribution<> distY(c,d);
    double strataSizeX = (b-a)/N;
    double strataSizeY = (d-c)/N;

    double** spectra = new double*[numW*2+1];
    for (int i = 0; i < numW*2+1; i++) {
        spectra[i] = new double[numW*2+1];
        for (int j = 0; j < numW*2+1; j++) {
            spectra[i][j] = 0;
        }
    }
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < N; j++) {
            sampleXs[j] = distX(gen);
            sampleYs[j] = distY(gen);
        }
        for (int wX = -numW; wX <= numW; wX++) {
            for (int wY = -numW; wY <= numW; wY++) {
                temp = (0,0);
                for (int j = 0; j < N; j++) {
                    temp += exp(-2 * M_PI * wStep * (wX * sampleXs[j] + wY * sampleYs[j]) * complex<double>(0,1));
                }
               // cout << temp << endl;
                //cout << (-wY + numW) << " " << (wX + numW) << endl;
                spectra[-wY + numW][wX + numW] += norm(temp) / (N * N * numTrials);
                //cout << "point done" << endl;
            }
            //cout << "strip done" << endl;
        }
    }
    return spectra;
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

double** uniform2DFourierCoefs(double a, double b, double c, double d, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    int numRows = (int) sqrt(N);
    double sampleXs[numRows*numRows] = {0};
    double sampleYs[numRows*numRows] = {0};
    uniform_real_distribution<> dist(0, 1);
    double strataSizeX = (b-a)/numRows;
    double strataSizeY = (d-c)/numRows;

    double** spectra = new double*[numW*2+1];
    for (int i = 0; i < numW*2+1; i++) {
        spectra[i] = new double[numW*2+1];
        for (int j = 0; j < numW*2+1; j++) {
            spectra[i][j] = 0;
        }
    }
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
        for (int k = 0; k < numRows; k++) {
            for (int j = 0; j < numRows; j++) {
                sampleXs[k * numRows + j] = a + (j + 0.5) * strataSizeX;
                sampleYs[k * numRows + j] = c + (k + 0.5) * strataSizeY;
            }
        }
        for (int wX = -numW; wX <= numW; wX++) {
            for (int wY = -numW; wY <= numW; wY++) {
                temp = (0,0);
                for (int j = 0; j < numRows*numRows; j++) {
                    temp += exp(-2 * M_PI * wStep * (wX * sampleXs[j] + wY * sampleYs[j]) * complex<double>(0,1));
                }
                spectra[-wY + numW][wX + numW] += norm(temp) / (numRows * numTrials);//extra gone as normalization
            }
        }
    }
    return spectra;
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

double** stratified2DFourierCoefs(double a, double b, double c, double d, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    int numRows = (int) sqrt(N);
    double sampleXs[numRows*numRows] = {0};
    double sampleYs[numRows*numRows] = {0};
    uniform_real_distribution<> dist(0, 1);
    double strataSizeX = (b-a)/numRows;
    double strataSizeY = (d-c)/numRows;

    double** spectra = new double*[numW*2+1];
    for (int i = 0; i < numW*2+1; i++) {
        spectra[i] = new double[numW*2+1];
        for (int j = 0; j < numW*2+1; j++) {
            spectra[i][j] = 0;
        }
    }
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
        for (int k = 0; k < numRows; k++) {
            for (int j = 0; j < numRows; j++) {
                sampleXs[k * numRows + j] = a + (j + dist(gen)) * strataSizeX;
                sampleYs[k * numRows + j] = c + (k + dist(gen)) * strataSizeY;
            }
        }
        /*if (i == 0) {
            for (int k = 0; k < numRows*numRows; k++) {
                cout << sampleXs[k] << ", " << sampleYs[k] << endl;
            }
        }*/
        for (int wX = -numW; wX <= numW; wX++) {
            for (int wY = -numW; wY <= numW; wY++) {
                temp = (0,0);
                for (int j = 0; j < numRows * numRows; j++) {
                    temp += exp(-2 * M_PI * wStep * (wX * sampleXs[j] + wY * sampleYs[j]) * complex<double>(0,1));
                }
               // cout << temp << endl;
                //cout << (-wY + numW) << " " << (wX + numW) << endl;
                spectra[-wY + numW][wX + numW] += norm(temp) / (numRows * numRows * numTrials);//extra gone as normalization
                //cout << "point done" << endl;
            }
            //cout << "strip done" << endl;
        }
    }
    return spectra;
}

double NRooks2D(double a, double b, double c, double d, int N, double lambda, function<double (double,double,double)> F, mt19937 & gen) {
    //somewhat based on pseudo code from https://cs.dartmouth.edu/~wjarosz/publications/subr16fourier.html
    uniform_real_distribution<> dist(0,1);
    double estimate = 0;
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

    for (int i = 0; i < N; i++) {
        estimate += F(sampleXs[i],sampleYs[i],lambda);
    }
    estimate /= N;
    estimate *= (b-a);
    estimate *= (d-c);
    return estimate;
}

double** NRooks2DFourierCoefs(double a, double b, double c, double d, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    uniform_real_distribution<> dist(0,1);
    double estimate = 0;
    double strataSizeX = (b-a)/N;
    double strataSizeY = (d-c)/N;
    vector<double> sampleXs(N,0);
    vector<double> sampleYs(N,0);

    double** spectra = new double*[numW*2+1];
    for (int i = 0; i < numW*2+1; i++) {
        spectra[i] = new double[numW*2+1];
        for (int j = 0; j < numW*2+1; j++) {
            spectra[i][j] = 0;
        }
    }
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < N; j++) {
            sampleXs[j] = a + (j + dist(gen)) * strataSizeX;
            sampleYs[j] = c + (j + dist(gen)) * strataSizeY;
        }

        random_shuffle(sampleXs.begin(),sampleXs.end());//Shuffle the sample x coordinates
        random_shuffle(sampleYs.begin(),sampleYs.end());//Shuffle the sample y coordinates

        /*if (i == 0) {
            for (int k = 0; k < numRows*numRows; k++) {
                cout << sampleXs[k] << ", " << sampleYs[k] << endl;
            }
        }*/
        for (int wX = -numW; wX <= numW; wX++) {
            for (int wY = -numW; wY <= numW; wY++) {
                temp = (0,0);
                for (int j = 0; j < N; j++) {
                    temp += exp(-2 * M_PI * wStep * (wX * sampleXs[j] + wY * sampleYs[j]) * complex<double>(0,1));
                }
               // cout << temp << endl;
                //cout << (-wY + numW) << " " << (wX + numW) << endl;
                spectra[-wY + numW][wX + numW] += norm(temp) / (N * numTrials);//extra gone as normalization
                //cout << "point done" << endl;
            }
            //cout << "strip done" << endl;
        }
    }
    return spectra;
}

double multiJitter2D(double a, double b, double c, double d, int N, double lambda, function<double (double,double,double)> F, mt19937 & gen) {
    //somewhat based on pseudo code from https://cs.dartmouth.edu/~wjarosz/publications/subr16fourier.html
    //Has N*N samples, not N samples
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

    /*cout << "Start: " << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << sampleXs[i][j] << ", " << sampleYs[j][i] << endl;
        }
    }*/

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            estimate += F(sampleXs[i][j],sampleYs[j][i],lambda);
        }
    }
    estimate /= N*N;
    estimate *= (b-a);
    estimate *= (d-c);
    return estimate;
}

double** multiJitter2DFourierCoefs(double a, double b, double c, double d, int N, int numTrials, double wStep, int numW, mt19937 & gen) {
    //Has N*N samples, not N samples
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

    double** spectra = new double*[numW*2+1];
    for (int i = 0; i < numW*2+1; i++) {
        spectra[i] = new double[numW*2+1];
        for (int j = 0; j < numW*2+1; j++) {
            spectra[i][j] = 0;
        }
    }
    complex<double> temp(0,0);
    for (int i = 0; i < numTrials; i++) {
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

        for (int wX = -numW; wX <= numW; wX++) {
            for (int wY = -numW; wY <= numW; wY++) {
                temp = (0,0);
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < N; k++) {
                        temp += exp(-2 * M_PI * wStep * (wX * sampleXs[j][k] + wY * sampleYs[k][j]) * complex<double>(0,1));
                    }
                }
               // cout << temp << endl;
                //cout << (-wY + numW) << " " << (wX + numW) << endl;
                spectra[-wY + numW][wX + numW] += norm(temp) / (N * N * numTrials);//extra gone as normalization
                //cout << "point done" << endl;
            }
            //cout << "strip done" << endl;
        }
    }
    return spectra;
}

int main(int argc, char** argv) {
    string fileName = argv[1];
    random_device r;
    mt19937 gen(r());
    int imgWidth = 500;
    int imgHeight = 500;
    int numSamples = 20;
    int numTrials = 100;
    int numLambdas = 250;
    unsigned char image[imgHeight*imgWidth];
    for (int i = 0; i < imgHeight*imgWidth; i++) {
        image[i] = 0;
    }
    double avgError = 0;
    double avgError2 = 0;
    double temp;
    
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
            //cout << abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas) << endl;
            avgError2 += (temp-groundTruthQuad(j/numLambdas)) * (temp-groundTruthQuad(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    /*for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < numLambdas; j++) {
                temp = pureMonteCarlo(0,1,numSamples,j/numLambdas,sampleQuad,gen);
                //image[(unsigned int) ((i + 110)*imgWidth + j + 20)] = temp;
                avgError += abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas);
                //cout << abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas) << endl;
                avgError2 += (temp-groundTruthQuad(j/numLambdas)) * (temp-groundTruthQuad(j/numLambdas));
            }
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/(numTrials*numLambdas)) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }*/
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < numLambdas; j++) {
                temp = uniform(0,1,numSamples,j/numLambdas,sampleQuad);
                //image[(unsigned int) ((i + 110)*imgWidth + j + 20)] = temp;
                avgError += abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas);
                //cout << abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas) << endl;
                avgError2 += (temp-groundTruthQuad(j/numLambdas)) * (temp-groundTruthQuad(j/numLambdas));
            }
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/(numTrials*numLambdas)) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }*/
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = uniform(0,1,numSamples,j/numLambdas,sampleQuad);
            //image[(unsigned int) ((i + 220)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas);
            avgError2 += (temp-groundTruthQuad(j/numLambdas)) * (temp-groundTruthQuad(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < numLambdas; j++) {
                temp = stratified(0,1,numSamples,j/numLambdas,sampleQuad,gen);
                //image[(unsigned int) ((i + 110)*imgWidth + j + 20)] = temp;
                avgError += abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas);
                //cout << abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas) << endl;
                avgError2 += (temp-groundTruthQuad(j/numLambdas)) * (temp-groundTruthQuad(j/numLambdas));
            }
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/(numTrials*numLambdas)) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }*/
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = stratified(0,1,numSamples,j/numLambdas,sampleQuad,gen);
            //image[(unsigned int) ((i + 330)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthQuad(j/numLambdas))/groundTruthQuad(j/numLambdas);
            avgError2 += (temp-groundTruthQuad(j/numLambdas)) * (temp-groundTruthQuad(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    cout << endl;

    cout << "Step:\n";
    for (int i = 0; i < 100; i++) {
        for (double j = 0; j < numLambdas; j++) {
            image[(unsigned int) (i*imgWidth + j + 20)] = groundTruthStep(j/numLambdas);
        }
    }
    cout << "Monte Carlo:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = pureMonteCarlo(0,1,numSamples,j/numLambdas,sampleStep,gen);
            image[(unsigned int) ((i + 110)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthStep(j/numLambdas))/groundTruthStep(j/numLambdas);
            avgError2 += (temp-groundTruthStep(j/numLambdas)) * (temp-groundTruthStep(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    /*for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            for (double j = 0; j < numLambdas; j++) {
                temp = pureMonteCarlo(0,1,numSamples,j/numLambdas,sampleStep,gen);
                //image[(unsigned int) ((i + 110)*imgWidth + j + 20)] = temp;
                avgError += abs(temp-groundTruthStep(j/numLambdas))/groundTruthStep(j/numLambdas);
                avgError2 += (temp-groundTruthStep(j/numLambdas)) * (temp-groundTruthStep(j/numLambdas));
            }
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/(numTrials*numLambdas)) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }*/
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = uniform(0,1,numSamples,j/numLambdas,sampleStep);
            //cout << temp << " vs " << groundTruthStep(j/numLambdas) << endl;
            image[(unsigned int) ((i + 220)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthStep(j/numLambdas))/groundTruthStep(j/numLambdas);
            avgError2 += (temp-groundTruthStep(j/numLambdas)) * (temp-groundTruthStep(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        for (double j = 0; j < numLambdas; j++) {
            temp = stratified(0,1,numSamples,j/numLambdas,sampleStep,gen);
            image[(unsigned int) ((i + 330)*imgWidth + j + 20)] = temp;
            avgError += abs(temp-groundTruthStep(j/numLambdas))/groundTruthStep(j/numLambdas);
            avgError2 += (temp-groundTruthStep(j/numLambdas)) * (temp-groundTruthStep(j/numLambdas));
        }
    }
    cout << "RMSE: " << sqrt(avgError2/(numTrials*numLambdas)) << endl;
    cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    cout << endl;

    
    cout << "--------------------------" << endl;
    cout << endl;
    cout << "2D:" << endl;
    cout << "Disk:\n";
    cout << "Monte Carlo:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,disk,gen);
        avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
        avgError2 += (temp-groundTruthDisk(0)) * (temp-groundTruthDisk(0));
    }
    cout << "RMSE: " << sqrt(avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,disk,gen);
            avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
            avgError2 += (temp-groundTruthDisk(0)) * (temp-groundTruthDisk(0));
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,disk);
        avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
        avgError2 += (temp-groundTruthDisk(0)) * (temp-groundTruthDisk(0));
    }
    cout << "RMSE: " << sqrt(avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = uniform2D(0,1,0,1,numSamples,0,disk);
            avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
            avgError2 += (temp-groundTruthDisk(0)) * (temp-groundTruthDisk(0));
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,disk,gen);
        avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
        avgError2 += abs(temp-groundTruthDisk(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = stratified2D(0,1,0,1,numSamples,0,disk,gen);
            avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
            avgError2 += abs(temp-groundTruthDisk(0)) * abs(temp-groundTruthDisk(0));
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }
    cout << "N-Rooks:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = NRooks2D(0,1,0,1,numSamples,0,disk,gen);
        avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
        avgError2 += abs(temp-groundTruthDisk(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = NRooks2D(0,1,0,1,numSamples,0,disk,gen);
            avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
            avgError2 += abs(temp-groundTruthDisk(0)) * abs(temp-groundTruthDisk(0));
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }
    cout << "Multi-Jitter:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = multiJitter2D(0,1,0,1,(int) sqrt(numSamples),0,disk,gen);
        avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
        avgError2 += abs(temp-groundTruthDisk(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k ++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = multiJitter2D(0,1,0,1,(int) sqrt(numSamples),0,disk,gen);
            avgError += abs(temp-groundTruthDisk(0))/groundTruthDisk(0);
            avgError2 += abs(temp-groundTruthDisk(0)) * abs(temp-groundTruthDisk(0));
        }
        //RMSE
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
        //cout << "Average percent error: " << (100*avgError/(numTrials*numLambdas)) << "%\n";
    }
    cout << endl;

    cout << "Triangle:\n";
    avgError = 0;
    avgError2 = 0;
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,triangle,gen);
        avgError += abs(temp-groundTruthTriangle(0))/groundTruthTriangle(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,triangle);
        avgError += abs(temp-groundTruthTriangle(0))/groundTruthTriangle(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,triangle,gen);
        avgError += abs(temp-groundTruthTriangle(0))/groundTruthTriangle(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "N-Rooks:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = NRooks2D(0,1,0,1,numSamples,0,triangle,gen);
        avgError += abs(temp-groundTruthTriangle(0))/groundTruthTriangle(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Multi-Jitter:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = multiJitter2D(0,1,0,1,(int) sqrt(numSamples),0,triangle,gen);
        avgError += abs(temp-groundTruthTriangle(0))/groundTruthTriangle(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << endl;

    cout << "Step:\n";
    avgError = 0;
    avgError2 = 0;
    cout << "Monte Carlo:\n";
    for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,step2D,gen);
        avgError += abs(temp-groundTruthStep2D(0))/groundTruthStep2D(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,step2D);
        avgError += abs(temp-groundTruthStep2D(0))/groundTruthStep2D(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
   cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,step2D,gen);
        avgError += abs(temp-groundTruthStep2D(0))/groundTruthStep2D(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "NRooks:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = NRooks2D(0,1,0,1,numSamples,0,step2D,gen);
        avgError += abs(temp-groundTruthStep2D(0))/groundTruthStep2D(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Multi-Jitter:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = multiJitter2D(0,1,0,1,(int) sqrt(numSamples),0,step2D,gen);
        avgError += abs(temp-groundTruthStep2D(0))/groundTruthStep2D(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << endl;

    cout << "Gaussian:\n";
    cout << "Monte Carlo:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,gaussian,gen);
        avgError += abs(temp-groundTruthGaussian(0))/groundTruthGaussian(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,gaussian,gen);
            avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
        }
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
    }
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,gaussian);
        avgError += abs(temp-groundTruthGaussian(0))/groundTruthGaussian(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = uniform2D(0,1,0,1,numSamples,0,gaussian);
            avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
        }
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
    }
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,gaussian,gen);
        avgError += abs(temp-groundTruthGaussian(0))/groundTruthGaussian(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = stratified2D(0,1,0,1,numSamples,0,gaussian,gen);
            avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
        }
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
    }
    cout << "N-Rooks:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = NRooks2D(0,1,0,1,numSamples,0,gaussian,gen);
        avgError += abs(temp-groundTruthGaussian(0))/groundTruthGaussian(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = NRooks2D(0,1,0,1,numSamples,0,gaussian,gen);
            avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
        }
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
    }
    cout << "Multi-Jitter:\n";
    avgError = 0;
    avgError2 = 0;
    /*for (int i = 0; i < numTrials; i++) {
        temp = multiJitter2D(0,1,0,1,(int) sqrt(numSamples),0,gaussian,gen);
        avgError += abs(temp-groundTruthGaussian(0))/groundTruthGaussian(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";*/
    for (int k = 5; k <= 150; k++) {
        avgError = 0;
        avgError2 = 0;
        numSamples = k;
        for (int i = 0; i < numTrials; i++) {
            temp = multiJitter2D(0,1,0,1,(int) sqrt(numSamples),0,gaussian,gen);
            avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
        }
        cout << k << "," << sqrt(avgError2/numTrials) << endl;
    }
    cout << endl;

    cout << "Bilinear:\n";
    cout << "Monte Carlo:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = pureMonteCarlo2D(0,1,0,1,numSamples,0,bilinear,gen);
        avgError += abs(temp-groundTruthBilinear(0))/groundTruthBilinear(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Uniform:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = uniform2D(0,1,0,1,numSamples,0,bilinear);
        avgError += abs(temp-groundTruthBilinear(0))/groundTruthBilinear(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Stratified:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = stratified2D(0,1,0,1,numSamples,0,bilinear,gen);
        avgError += abs(temp-groundTruthBilinear(0))/groundTruthBilinear(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "N-Rooks:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = NRooks2D(0,1,0,1,numSamples,0,bilinear,gen);
        avgError += abs(temp-groundTruthBilinear(0))/groundTruthBilinear(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << "Multi-jitter:\n";
    avgError = 0;
    avgError2 = 0;
    for (int i = 0; i < numTrials; i++) {
        temp = multiJitter2D(0,1,0,1,(int) sqrt(numSamples),0,bilinear,gen);
        avgError += abs(temp-groundTruthBilinear(0))/groundTruthBilinear(0);
        avgError2 += (temp-groundTruthTriangle(0)) * (temp-groundTruthTriangle(0));
    }
    cout << "Average error: " << (avgError2/numTrials) << endl;
    cout << "Average percent error: " << (100*avgError/numTrials) << "%\n";
    cout << endl;

    //Fourier analysis
    /*double* testCoefs = stratifiedFourierCoefs(0,1,numSamples,100000,1,40,gen);
    for (int i = 0; i < 40; i++) {
        cout << (i) <<  "," << testCoefs[i]*numSamples << "\n";
    }
    cout << endl;
    delete testCoefs;*/

    //Power spectra
    /*double** testSpectra = multiJitter2DFourierCoefs(0,1,0,1,15,1000,1,60,gen);
    //cout << "Hi" << endl;
    for (int i = 0; i < 121; i++) {
        //cout << "Hello???" << endl;
        for (int j = 0; j < 121; j++) {
            if (i == 5) {
                //cout << numSamples * testSpectra[i][j] << " ";
            }
            //cout << testSpectra[i][j] << " ";
            image[i * imgWidth + j] = (unsigned char) min(255, (int) (200 * testSpectra[i][j]));
        }
        //cout << endl;
    }
    for (int i = 0; i < 33; i++) {
        delete testSpectra[i];
    }
    delete testSpectra;*/

    //Make image file
    ofstream ofs(fileName, ios::out | ios::binary);
    ofs << "P6\n" << imgWidth << " " << imgHeight << "\n255\n"; 
    for (int i = 0; i < imgWidth * imgHeight; ++i) { 
        ofs << image[i] << image[i] << image[i];
    } 
    ofs.close();
}