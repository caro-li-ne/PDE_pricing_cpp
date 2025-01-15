#ifndef PDE_H
#define PDE_H

#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "Matrix.h"

class PDEPricer {
private:
    vector<double> timeGrid;
    vector<double> spaceGrid;
    vector<double> terminalCondition;
    vector<double> boundaryLow, boundaryHigh;
    vector<double> A, B, C, D;
    Matrix P, Q, V;
    double S0, T, K, r, sigma, dx, dt, theta, lambda;
    int spaceSteps, timeSteps;

public:

    PDEPricer(double S0, double maturity, double K, int timeSteps, int spaceSteps, double volatility, double rate, double lambda, double theta);

    void computeA();
    void computeB();
    void computeC();
    void computeD();
    void setBoundaryConditions(const vector<double>& low, const vector<double>& high);
    void setTerminalCondition(const vector<double>& terminal);
    void computeMatrices();
    Matrix solve();
};

double blackScholesCall(double S, double K, double r, double sigma, double T);

#endif // PDE_H

