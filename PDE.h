#ifndef PDE_H
#define PDE_H

#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "Matrix.h"

class PDEPricer {
private:
    std::vector<double> timeGrid;
    std::vector<double> spaceGrid;
    std::vector<double> terminalCondition;
    std::vector<double> boundaryLow, boundaryHigh;
    std::vector<double> A, B, C, D;
    Matrix P, Q, V;
    double S0, T, K, r, sigma, dx, dt, theta, lambda;
    int spaceSteps, timeSteps;

public:
//PDEPricer(double maturity, int timeSteps, int spaceSteps, double volatility, double rate, double lambda);

    PDEPricer(double S0, double maturity, double K, int timeSteps, int spaceSteps, double volatility, double rate, double lambda, double theta);

    void computeA();
    void computeB();
    void computeC();
    void computeD();
    void setBoundaryConditions(const std::vector<double>& low, const std::vector<double>& high);
    void setTerminalCondition(const std::vector<double>& terminal);
    void computeMatrices();
   //void compute();
    Matrix solve();
};

double blackScholesCall(double S, double K, double r, double sigma, double T);

#endif // PDE_H

