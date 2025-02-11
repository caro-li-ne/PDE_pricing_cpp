# Pricing vanilla options through PDEs

The aim of this project is to build a pricer of vanilla option from the PDE of the derivatives payoff. Key implementations include matrix inversion algorithm and Crank-Nicholson scheme resolution. Comments are written along the way. *(last update: Jan 13th)*


*Matrix inversion*: uses the **Gauss** pivot algorithm to find the inverse of a matrix. Checks for legality of inversion *(matrix is squared and determinant is not zero)*. In case of non invertibility, epsilon added to make inversion possible.

```cpp
matrix::inverse()
```

*Matrix determinant calculator*: computes determinant. Used to check if a matrix can be inverted.

```cpp
matrix::computeDeterminantGaussian()
```

*Theoretical call option pricing*: uses **Black-Scholes** formula to compute the price of a vanilla call option. Used to study the convergence of a price obtained through finite difference.

```cpp
blackScholesCall()
```

*Boundaries and terminal conditions*: Prepares the boundary and terminal conditions.

```cpp
PDEpricer::setBoundaryConditions()
```

```cpp
PDEpricer::setTerminalCondition()
```

*Matrix initialization*: prepares the Matrixes P and Q for our finite difference resolution. (Those are constant through the iteration so we only compute them once.)

```cpp
PDEpricer::computeMatrices()
```

*Finite difference call option pricing*: uses **Crank-Nicholson** scheme to build the whole final prices grid backwards in time.

```cpp
PDEpricer::solve()
```

# Results

All results and conclusion are documented in the pdf **CPP_PDE_Pricer.pdf**. Have a nice day and enjoy your reading :relaxed:!

