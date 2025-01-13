#include "PDE.h"
#include "Matrix.h"
#include <iostream>
#include <vector>
#include <cmath>

int main() {
    std::cout << "Program started..." << std::endl;

    // Parameters
    double maturity = 1; // Time to maturity (T, in years)
    int timeSteps = 100; // Number of time steps (n)
    int spaceSteps = 400; // Number of space steps (m)
    double sigma = 0.15; // Volatility
    double r = 0.03; // Risk-free interest rate
    double lambda = 3.0; // Multiplier for space grid boundaries
    double strike = 40.0; // Strike price
    double S_0 = 70; // Current stock price 
    double q = 0; // Dividend yield
    double theta = 0.5;

double dt = maturity / timeSteps;

// Compute forward price
    double F_0_T = S_0 * exp((r - q) * maturity);

// Print forward price
std::cout << "Forward Price F_0_T: " << F_0_T << std::endl;

// Uniform space grid with boundaries -lambda*sigma to +lambda*sigma
std::vector<double> spaceGrid(spaceSteps);
double dx = (2 * lambda * sigma) / (spaceSteps - 1); // Step size
for (int j = 0; j < spaceSteps; ++j) {
    spaceGrid[j] = - lambda * sigma + j * dx;
}

double dx2 = spaceGrid[1] - spaceGrid[0];
double stability_criterion = dx2 * dx2 / (sigma * sigma);
if (dt > stability_criterion) {
    std::cerr << "Warning: Time step size exceeds stability limit!" <<dt<<stability_criterion<< std::endl;
}

std::cout << "Initializing PDEPricer..." << std::endl;

PDEPricer pricer(S_0, maturity, strike, timeSteps, spaceSteps, sigma, r, lambda, theta);

// Terminal condition: max(S - K, 0) -> aligned with spaceGrid
std::vector<double> terminal(spaceSteps,0.0);
for (int j = 0; j < spaceSteps; ++j) {
    double S = exp(spaceGrid[j]) * S_0;// Convert back to S
    terminal[j] = std::max(S - strike * exp(- r * maturity), 0.0);
}

// Boundary conditions
std::vector<double> boundaryLow(timeSteps); 
std::vector<double> boundaryHigh(timeSteps);

for (int i = 0; i < timeSteps; ++i) {
    double t = i * maturity / (timeSteps - 1);
    
    double S_min = S_0 * exp(spaceGrid.front());
    S_min = std::max(S_min - strike * exp(- r * (maturity - t)), 0.0);
    boundaryLow[i] = S_min; // At S = -lambda*sigma, u(t, S_min) = 0

    double S_max = S_0 * exp(spaceGrid.back());
    boundaryHigh[i] = S_max - strike * exp(- r * (maturity - t));
    }

// Set boundary conditions in the PDE pricer
std::cout << "Setting boundary and terminal conditions..." << std::endl;
pricer.setBoundaryConditions(boundaryLow, boundaryHigh);
pricer.setTerminalCondition(terminal);

// Compute Matrix P and Q 
std::cout << "Computing matrices..." << std::endl;
pricer.computeMatrices();

// Solving PDE backwards
std::cout << "Solving PDE..." << std::endl;
Matrix prices(spaceSteps,1);
prices = pricer.solve();

// Compute with Black-Scholes price at S_0 = 50.0
double bsPrice = blackScholesCall(S_0, strike, r, sigma, maturity);

// Interpolate PDE price at S0, calculated grid index may not align perfectly with S0=50 due to discretisation, 
//this mismatch could lead to the PDE price being evaluated at a point far from S0.
int idx = static_cast<int> ((log(S_0 / strike) + lambda * sigma)* spaceSteps / (2 * lambda * sigma)) ;
idx = std::max(1, std::min(idx, spaceSteps - 2));

double S_left = exp(spaceGrid[idx - 1]) * S_0;
double S_right = exp(spaceGrid[idx]) * S_0;
double price = prices.at(idx - 1,0) + (prices.at(idx,0) - prices.at(idx - 1,0)) * (S_0 - S_left) / (S_right - S_left);

std::cout << "Interpolated PDE Price at S0: " << price <<std::endl;
std::cout << "Black-Scholes Price: " << bsPrice << std::endl;

return 0;
}
