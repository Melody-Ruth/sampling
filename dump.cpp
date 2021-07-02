//Old code from polyTest.cpp that I don't think I want anymore but don't want to delete yet.

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
double pureMonteCarlo2D(double a, double b, double c, double d, int N, function<double (double,double)> F, mt19937 & gen) {
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
        temp = F(distX(gen),distY(gen));
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
                spectra[-wY + numW][wX + numW] += norm(temp) / (N * N * numTrials);
            }
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

double uniform2D(double a, double b, double c, double d, int N, function<double (double,double)> F) {
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
            estimate += F(tempX,tempY);
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

double stratified2D(double a, double b, double c, double d, int N, function<double (double,double)> F, mt19937 & gen) {
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
            estimate += F(tempX,tempY);
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
        for (int wX = -numW; wX <= numW; wX++) {
            for (int wY = -numW; wY <= numW; wY++) {
                temp = (0,0);
                for (int j = 0; j < numRows * numRows; j++) {
                    temp += exp(-2 * M_PI * wStep * (wX * sampleXs[j] + wY * sampleYs[j]) * complex<double>(0,1));
                }
                spectra[-wY + numW][wX + numW] += norm(temp) / (numRows * numRows * numTrials);//extra gone as normalization
            }
        }
    }
    return spectra;
}

double NRooks2D(double a, double b, double c, double d, int N, function<double (double,double)> F, mt19937 & gen) {
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
        estimate += F(sampleXs[i],sampleYs[i]);
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

        for (int wX = -numW; wX <= numW; wX++) {
            for (int wY = -numW; wY <= numW; wY++) {
                temp = (0,0);
                for (int j = 0; j < N; j++) {
                    temp += exp(-2 * M_PI * wStep * (wX * sampleXs[j] + wY * sampleYs[j]) * complex<double>(0,1));
                }
                spectra[-wY + numW][wX + numW] += norm(temp) / (N * numTrials);//extra gone as normalization
            }
        }
    }
    return spectra;
}

double multiJitter2D(double a, double b, double c, double d, int N, function<double (double,double)> F, mt19937 & gen) {
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

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            estimate += F(sampleXs[i][j],sampleYs[j][i]);
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
                spectra[-wY + numW][wX + numW] += norm(temp) / (N * N * numTrials);//extra gone as normalization
            }
        }
    }
    return spectra;
}