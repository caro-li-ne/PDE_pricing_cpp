#ifndef MATRIX2_H
#define MATRIX2_H
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
    Matrix(int rows, int cols) : mat(rows, vector<double>(cols, 0.0)), rows(rows), cols(cols) {}

    // Getter for matrix dimensions
    int getRows() const { return rows; }
    int getCols() const { return cols; }

    // Access matrix element with bounds checking
    double& at(int row, int col) {
        if (row < 0 || row >= rows || col < 0 || col >= cols) {
            throw std::out_of_range("Matrix indices out of bounds");
        }
        return mat[row][col];
    }

    // Const access to matrix elements
    const double& at(int row, int col) const {
        if (row < 0 || row >= rows || col < 0 || col >= cols) {
            throw std::out_of_range("Matrix indices out of bounds");
        }
        return mat[row][col];
    }

    // Function to set the matrix elements
    void input() {
        cout << "Enter the elements of the matrix row by row:\n";
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                cin >> mat[i][j];
            }
        }
    }

    // Function to check if the matrix is invertible
    bool isInvertible() {
        double epsilon =1e-8;
        cout<<"Determinant : "<<determinant()<<endl;
        return abs(determinant()) > epsilon;
    }

    Matrix inverse() {
        if (rows != cols) {
            throw std::runtime_error("Matrix must be square to compute its inverse.");
        }

        int n = rows;
        Matrix result(n, n); // Identity matrix to store the inverse
        std::vector<std::vector<double>> augmented(n, std::vector<double>(2 * n));

        // Initialize the augmented matrix [A | I]
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                augmented[i][j] = mat[i][j];
                augmented[i][j + n] = (i == j) ? 1.0 : 0.0; // Identity matrix
            }
        }

        // Perform Gaussian elimination with partial pivoting and diagonal correction
        for (int i = 0; i < n; ++i) {
            // Find the pivot
            int pivot = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(augmented[k][i]) > std::abs(augmented[pivot][i])) {
                    pivot = k;
                }
            }

            // Swap rows if necessary
            if (pivot != i) {
                std::swap(augmented[i], augmented[pivot]);
            }

            // Check if the pivot is near zero; add epsilon if necessary
            if (std::abs(augmented[i][i]) < 1e-16) {
                std::cout << "Adding epsilon to diagonal element at (" << i << ", " << i << ")\n";
                augmented[i][i] += 1e-16;
            }

            // Normalize the pivot row
            double pivotValue = augmented[i][i];
            for (int j = 0; j < 2 * n; ++j) {
                augmented[i][j] /= pivotValue;
            }

            // Eliminate below and above
            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    double factor = augmented[k][i];
                    for (int j = 0; j < 2 * n; ++j) {
                        augmented[k][j] -= factor * augmented[i][j];
                    }
                }
            }
        }

        // Extract the inverse matrix from the augmented matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result.at(i, j) = augmented[i][j + n];
            }
        }

        return result;
}



    // Function to display the matrix
    void display() const {
        for (const auto& row : mat) {
            for (double elem : row) {
                cout << setw(10) << elem << " ";
            }
            cout << endl;
        }
    }


double determinant() {
    if (!determinantComputed) {
        cachedDeterminant = computeDeterminantGaussian(); // Switch to Gaussian elimination
        determinantComputed = true;
    }
    return cachedDeterminant;
}

double computeDeterminantGaussian() {
    if (rows != cols) throw std::runtime_error("Matrix must be square to compute determinant.");

    double det = 1.0;
    std::vector<std::vector<double>> tempMat = mat;

    for (int i = 0; i < rows; i++) {
        int pivot = i;
        for (int j = i + 1; j < rows; j++) {
            if (abs(tempMat[j][i]) > abs(tempMat[pivot][i])) pivot = j;
        }

        if (pivot != i) {
            std::swap(tempMat[i], tempMat[pivot]);
            det *= -1;
        }

        if (tempMat[i][i] == 0) return 0;

        det *= tempMat[i][i];
        for (int j = i + 1; j < rows; j++) {
            double factor = tempMat[j][i] / tempMat[i][i];
            for (int k = i; k < cols; k++) {
                tempMat[j][k] -= factor * tempMat[i][k];
            }
        }
    }

    return det;
}


vector<vector<double>> adjoint() {
    vector<vector<double>> adj(rows, vector<double>(cols));

    if (rows == 1) {
        adj[0][0] = 1;
        return adj;
    }

    int n = rows;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Create the sub-matrix by excluding the current row and column
            vector<vector<double>> subMatrix(n - 1, vector<double>(n - 1));
            int sub_i = 0;
            for (int row = 0; row < n; row++) {
                if (row == i) continue; // Skip the current row
                int sub_j = 0;
                for (int col = 0; col < n; col++) {
                    if (col == j) continue; // Skip the current column
                    subMatrix[sub_i][sub_j] = mat[row][col];
                    sub_j++;
                }
                sub_i++;
            }

            // Calculate determinant of the sub-matrix
            Matrix tempMatrix(n - 1, n - 1);
            for (int subRow = 0; subRow < n - 1; ++subRow) {
                for (int subCol = 0; subCol < n - 1; ++subCol) {
                    tempMatrix.at(subRow, subCol) = subMatrix[subRow][subCol];
                }
            }
            double cofactorDeterminant = tempMatrix.computeDeterminantGaussian();

            // Calculate cofactor and assign to adjoint matrix
            adj[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * cofactorDeterminant;
        }
    }

    return adj;
}
/*
    // Overloaded * operator for matrix-vector multiplication
    std::vector<double> operator*(const std::vector<double>& vec) const {
        if (vec.size() != cols) {
            throw std::runtime_error("Matrix and vector size mismatch for multiplication.");
        }

        std::vector<double> result(rows, 0.0);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[i] += mat[i][j] * vec[j];
            }
        }
        return result;
    }*/
   Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::runtime_error("Matrix size mismatch for addition.");
        }
        Matrix result(rows, other.cols);
        #pragma omp parallel for
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
            double temp=0.0;
            for (int k = 0; k < cols; ++k) {
                temp += this->at(i,k) * other.at(k,j);

            }
            result.at(i,j)=temp;
        }}
        return result;
    };
    
    // Overloaded + operator for matrix addition
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::runtime_error("Matrix size mismatch for addition.");
        }

        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.at(i, j) = mat[i][j] + other.at(i, j);
            }
        }
        return result;
    }
};

#endif // MATRIX2_H
