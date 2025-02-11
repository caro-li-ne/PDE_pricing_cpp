#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

class Matrix {
private:
    vector<vector<double>> mat;
    int rows, cols;
    double cachedDeterminant = 0.0;
    bool determinantComputed = false;

public:
    // Constructor to initialize the matrix
    Matrix(int rows, int cols) ;

    // Getter for matrix dimensions
    int getRows() const;
    int getCols() const;

    // Access matrix element with bounds checking
    double& at(int row, int col) ;

    // Const access to matrix elements
    const double& at(int row, int col)  const;

    // Function to set the matrix elements
    void input();

    // Function to check if the matrix is invertible
    bool isInvertible() ;

    //Compute inverse
    Matrix inverse() ;

    // Function to display the matrix
    void display() const;

    double determinant();
    double computeDeterminantGaussian();        
    vector<vector<double>> adjoint();

    // Overloaded * operator for matrix multiplication
    Matrix operator*(const Matrix& other) const;
    
    // Overloaded + operator for matrix addition
    Matrix operator+(const Matrix& other) const;
};

#endif // MATRIX_H
