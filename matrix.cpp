#include "Matrix.h"

// Constructor
Matrix::Matrix(int rows, int cols) : mat(rows, vector<double>(cols, 0.0)), rows(rows), cols(cols) {}

// Getters
int Matrix::getRows() const {
    return rows;
}
int Matrix::getCols() const {
    return cols;
}

// Access elements
double& Matrix::at(int row, int col) {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw out_of_range("Matrix indices out of bounds");
    }
    return mat[row][col];
}

const double& Matrix::at(int row, int col) const {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw out_of_range("Matrix indices out of bounds");
    }
    return mat[row][col];
}

// Input elements
void Matrix::input() {
    cout << "Enter the elements of the matrix row by row:\n";
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cin >> mat[i][j];
        }
    }
}

// Check if invertible
bool Matrix::isInvertible() {
    double epsilon = 1e-8;
    cout << "Determinant : " << determinant() << endl;
    return abs(determinant()) > epsilon;
}

// Compute inverse
Matrix Matrix::inverse() {
    if (rows != cols) {
        throw runtime_error("Matrix must be square to compute its inverse.");
    }

    int n = rows;
    Matrix result(n, n);
    vector<vector<double>> augmented(n, vector<double>(2 * n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented[i][j] = mat[i][j];
            augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(augmented[k][i]) > abs(augmented[pivot][i])) {
                pivot = k;
            }
        }

        if (pivot != i) {
            swap(augmented[i], augmented[pivot]);
        }

        if (abs(augmented[i][i]) < 1e-16) {
            cout << "Adding epsilon to diagonal element at (" << i << ", " << i << ")\n";
            augmented[i][i] += 1e-16;
        }

        double pivotValue = augmented[i][i];
        for (int j = 0; j < 2 * n; ++j) {
            augmented[i][j] /= pivotValue;
        }

        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = augmented[k][i];
                for (int j = 0; j < 2 * n; ++j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result.at(i, j) = augmented[i][j + n];
        }
    }

    return result;
}

// Display the matrix
void Matrix::display() const {
    for (const auto& row : mat) {
        for (double elem : row) {
            cout << setw(10) << elem << " ";
        }
        cout << endl;
    }
}

// Compute determinant
double Matrix::determinant() {
    if (!determinantComputed) {
        cachedDeterminant = computeDeterminantGaussian();
        determinantComputed = true;
    }
    return cachedDeterminant;
}

double Matrix::computeDeterminantGaussian() {
    if (rows != cols) throw runtime_error("Matrix must be square to compute determinant.");

    double det = 1.0;
    vector<vector<double>> tempMat = mat;

    for (int i = 0; i < rows; i++) {
        int pivot = i;
        for (int j = i + 1; j < rows; j++) {
            if (abs(tempMat[j][i]) > abs(tempMat[pivot][i])) pivot = j;
        }

        if (pivot != i) {
            swap(tempMat[i], tempMat[pivot]);
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

// Compute adjoint
vector<vector<double>> Matrix::adjoint() {
    vector<vector<double>> adj(rows, vector<double>(cols));

    if (rows == 1) {
        adj[0][0] = 1;
        return adj;
    }

    int n = rows;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vector<vector<double>> subMatrix(n - 1, vector<double>(n - 1));
            int sub_i = 0;
            for (int row = 0; row < n; row++) {
                if (row == i) continue;
                int sub_j = 0;
                for (int col = 0; col < n; col++) {
                    if (col == j) continue;
                    subMatrix[sub_i][sub_j] = mat[row][col];
                    sub_j++;
                }
                sub_i++;
            }

            Matrix tempMatrix(n - 1, n - 1);
            for (int subRow = 0; subRow < n - 1; ++subRow) {
                for (int subCol = 0; subCol < n - 1; ++subCol) {
                    tempMatrix.at(subRow, subCol) = subMatrix[subRow][subCol];
                }
            }
            double cofactorDeterminant = tempMatrix.computeDeterminantGaussian();
            adj[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * cofactorDeterminant;
        }
    }

    return adj;
}

// Operator overloads
Matrix Matrix::operator*(const Matrix& other) const {
    if (cols != other.rows) {
        throw runtime_error("Matrix size mismatch for multiplication.");
    }
    Matrix result(rows, other.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            double temp = 0.0;
            for (int k = 0; k < cols; ++k) {
                temp += this->at(i, k) * other.at(k, j);
            }
            result.at(i, j) = temp;
        }
    }
    return result;
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw runtime_error("Matrix size mismatch for addition.");
    }

    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.at(i, j) = mat[i][j] + other.at(i, j);
        }
    }
    return result;
}
