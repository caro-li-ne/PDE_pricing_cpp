#include "PDE2.h"
//#include "Matrix2.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <iostream>

PDEPricer::PDEPricer(double S0, double maturity, double K, int timeSteps, int spaceSteps, double volatility, double rate, double lambda, double theta)
    : S0(S0), T(maturity), K(K), r(rate), sigma(volatility), timeSteps(timeSteps), spaceSteps(spaceSteps), theta(theta), lambda(lambda),
    P(spaceSteps - 1, spaceSteps - 1), Q(spaceSteps - 1, spaceSteps - 1), V(spaceSteps - 1, 1)  {
    dx = 2.0 * lambda * sigma / (spaceSteps - 1);
    dt = T / (timeSteps - 1);

    spaceGrid.resize(spaceSteps);
    for (int i = 0; i < spaceSteps; ++i) {
        spaceGrid[i] = - lambda * sigma + i * dx;
    }

    timeGrid.resize(timeSteps);
    for (int i = 0; i < timeSteps; ++i) {
        timeGrid[i] = i * dt;
    }
    // Initialize matrices P and Q
    P = Matrix(spaceSteps - 1, spaceSteps - 1);
    Q = Matrix(spaceSteps - 1, spaceSteps - 1);
    V = Matrix(spaceSteps - 1, 1);

    computeA();
    computeB();
    computeC();
    computeD();

    computeMatrices();

}


void PDEPricer::computeA() {
    A.assign(spaceSteps - 1, 0.0);
}

void PDEPricer::computeB() {
    B.assign(spaceSteps - 1, - 0.5 * sigma * sigma);
}

void PDEPricer::computeC() {
    C.assign(spaceSteps - 1, 0.5 * sigma * sigma);
}

void PDEPricer::computeD() {
    D.assign(spaceSteps - 1, 0.0);
}

void PDEPricer::computeMatrices() {
    for (int j = 0; j < spaceSteps - 1; ++j) {       
        if (j > 0) {
            P.at(j, j - 1) = - 0.5 * B[j] / dx + theta * C[j] / (dx * dx);
            Q.at(j, j - 1) = (1.0 - theta) * C[j] / (dx * dx);
        }

        P.at(j, j) = A[j] - (1.0 / dt + 2 * theta * C[j] / (dx * dx));
        Q.at(j, j) = 1.0 / dt - 2 * (1.0 - theta) * C[j] / (dx * dx);

        if (j < spaceSteps - 2) {
            P.at(j, j + 1) = 0.5 * B[j] / dx + theta * C[j] / (dx * dx);
            Q.at(j, j + 1) = (1.0 - theta) * C[j] / (dx * dx);        
    }
}}

Matrix PDEPricer::solve() {
    Matrix U(spaceSteps - 1, 1); // Initial condition
    for (int i=0;i < spaceSteps - 1;++i){
        U.at(i,0)=terminalCondition[i];
        //terminal values initialized
    }
    double U_first_0=boundaryLow[timeSteps-1];
    double U_first_m=boundaryHigh[timeSteps-1];
    double U_prev_0, U_prev_m;
    Matrix P_inv = P.inverse();

    for (int n = timeSteps - 2; n > 0; --n) {
        U_prev_0=U_first_0;
        U_prev_m=U_first_m;
        U_first_0=boundaryLow[n];
        U_first_m=boundaryHigh[n];

        // Compute V
        for (int j = 0; j < spaceSteps - 1; ++j) {
            V.at(j,0) = D[j] ; 
            if (j==0){  // Compute V including boundary conditions
            V.at(j,0) += (- 0.5 * B[0] / dx + theta * C[0] / (dx * dx)) * U_first_0 + (1 - theta) * C[0] / (dx * dx) * U_prev_0;
            }
            if (j==spaceSteps-2){
            V.at(j, 0) += (0.5 * B[spaceSteps - 2] / dx + theta * C[spaceSteps - 2] / (dx * dx)) * U_first_m
            + (1 - theta) * C[spaceSteps - 2] / (dx * dx) *  U_prev_m ;
            }}

        //Compute the new U
        Matrix product(spaceSteps - 1,spaceSteps - 1);
        product = P_inv *(Q * U + V); // Assuming Q * U is valid

        for (size_t i = 0; i < spaceSteps - 1; ++i) {
            U.at(i,0) = -1.0 * product.at(i,0);
            }
    }
    return U;
}

// Note: `theta` is the parameter determining the numerical scheme:
// theta = 0 corresponds to explicit scheme
// theta = 1 corresponds to implicit scheme
// 0 < theta < 1 corresponds to Crank-Nicolson scheme.

void PDEPricer::setBoundaryConditions(const std::vector<double>& low, const std::vector<double>& high) {
    boundaryLow = low;
    boundaryHigh = high;
}

void PDEPricer::setTerminalCondition(const std::vector<double>& terminal) {
    terminalCondition = terminal;
}

double blackScholesCall(double S, double K, double r, double sigma, double T) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    return S * 0.5 * erfc(-d1 / sqrt(2)) - K * exp(-r * T) * 0.5 * erfc(- d2/sqrt(2));
}
