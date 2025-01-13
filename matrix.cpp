#if !defined(MATRIX_H)
#define MATRIX_H
#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <math.h>
#include <stdlib.h>
#define ZERO 1e-8

#include <vector>
using namespace std;
class Matrix
{
    private:
        int n;
        char m_name[128];

        std::vector<std::vector<double>> mat;
        Matrix();
            double computeDeterminant(const vector<vector<double>>& matrix, int size) {
        if (size == 1) {
            return matrix[0][0];
        }

        if (size == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }

        double det = 0;
        std::vector<std::vector<double>> subMatrix(size - 1, std::vector<double>(size - 1));

        for (int x = 0; x < size; x++) {
            // Create the sub-matrix
            for (int i = 1; i < size; i++) {
                int sub_j = 0;
                for (int j = 0; j < size; j++) {
                    if (j == x) continue;
                    subMatrix[i - 1][sub_j] = matrix[i][j];
                    sub_j++;
                }
            }

            // Recursive call
            det += (x % 2 == 0 ? 1 : -1) * matrix[0][x] * computeDeterminant(subMatrix, size - 1);
        }

        return det;
    };
    public:
        double **m_pData;

    double determinant() {
        return computeDeterminant(mat, n);
    }

    Matrix(int size) : n(size), mat(size, vector<double>(size)) 
        /*Matrix(const char *name, int rows, int cols) :
            n(rows), n(cols)*/
        {
            strcpy(m_name, name);
            m_pData = new double*[n];
            for (int i = 0; i < n; i++)
                m_pData[i] = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    m_pData[i][j] = 0.0;
                }
            }
        }
        Matrix(const Matrix &other)
        {
            strcpy(m_name, other.m_name);
            n = other.n;
            m_pData = new double*[n];
            for (int i = 0; i < n; i++)
                m_pData[i] = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    m_pData[i][j] = other.m_pData[i][j];
                }
            }
        }
        ~Matrix()
        {
            for (int i = 0; i < n; i++)
                delete[] m_pData[i];
            delete[] m_pData;
            n = 0;
        }
        void SetName(const char *name)
        {
            strcpy(m_name, name);
        }
        const char* GetName() const
        {
            return m_name;
        }
        void GetInput()
        {
            std::cin >> *this;
        }
        void FillSimulatedInput()
        {
            //static int factor1 = 1, factor2 = 2;
            std::cout << "\n\nEnter Input For Matrix : " << std::endl;
            std::cin>>m_name;
            std::cout<<"Matrix Dimension : "<<std::endl;
            int nb;
            std::cin>> nb;
            n=nb;
            n=nb;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    std::cout << "Input For Row: " << i + 1 << " Col: " << j
                            + 1 << " = ";
                    std::cin>>m_pData[i][j];
                    //int data = ((i + 1) * factor1) + (j + 1) * factor2;
                    //m_pData[i][j] = data / 10.2;
                    //std::cout << m_pData[i][j] << "\n";
                    //factor1 += (rand() % 4);
                    //factor2 += (rand() % 3);
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }

/*
        double Determinant(int a[MAX][MAX],int n) {
        
        double det = 0;
        double **a = m_pData;

        int p, h, k, i, j, temp[n][n];
  if(n==1) {
    return a[0][0];
  } else if(n==2) {
    det=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
    return det;
  } else {
    for(p=0;p<n;p++) {
      h = 0;
      k = 0;
      for(i=1;i<n;i++) {
        for( j=0;j<n;j++) {
          if(j==p) {
            continue;
          }
          temp[h][k] = a[i][j];
          k++;
          if(k==n-1) {
            h++;
            k = 0;
          }
        }
      }
      det=det+a[0][p]*pow(-1,p)*Determinant(temp,n-1);
    }
    return det;
  }
} */

double determinant() {
    
    double **mat = m_pData;
    if (n == 1) {
        return mat[0][0];
    }

    if (n == 2) {
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }

    double det = 0;
    Matrix subMatrix(n - 1, vector<double>(n - 1));

    for (int x = 0; x < n; x++) {
        // Create the sub-matrix
        for (int i = 1; i < n; i++) {
            int sub_j = 0;
            for (int j = 0; j < n; j++) {
                if (j == x) continue;
                subMatrix[i - 1][sub_j] = mat[i][j];
                sub_j++;
            }
        }

        // Recursive call
        det += (x % 2 == 0 ? 1 : -1) * mat[0][x] * determinant(subMatrix, n - 1);
    }

    return det;
}



double Determinant() {
    int n=n;
    double **a = m_pData;
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(a[j][i]) > abs(a[pivot][i])) {
                pivot = j;
            }
        }
        if (pivot != i) {
            std::swap(a[i], a[pivot]);
            det *= -1;
        }
        if (a[i][i] == 0) {
            det *= 0;
        }
        det *= a[i][i];
        for (int j = i + 1; j < n; j++) {
            double factor = a[j][i] / a[i][i];
            for (int k = i + 1; k < n; k++) {
                a[j][k] -= factor * a[i][k];
            }
        }
    }
    return det;
}



        /* double Determinant()
        {
            double det = 0;
            double **pd = m_pData;
            switch (n)
            {
                /*
                case 2:
                {
                    det = pd[0][0] * pd[1][1] - pd[0][1] * pd[1][0];
                    return det;
                }
                    break;
                case 3:
                {
                    /***
                     a b c
                     d e f
                     g h i
 
                     a b c a b c
                     d e f d e f
                     g h i g h i
 
                     // det (A) = aei + bfg + cdh - afh - bdi - ceg.
                     ***/
                  /*  double a = pd[0][0];
                    double b = pd[0][1];
                    double c = pd[0][2];
                    double d = pd[1][0];
                    double e = pd[1][1];
                    double f = pd[1][2];
                    double g = pd[2][0];
                    double h = pd[2][1];
                    double i = pd[2][2];
                    double det = (a * e * i + b * f * g + c * d * h);
                    det = det - a * f * h;
                    det = det - b * d * i;
                    det = det - c * e * g;
                    return det;
                }
                    break;
                case 4:
                {
                    Matrix *temp[4];
                    for (int i = 0; i < 4; i++)
                        temp[i] = new Matrix("", 3, 3);
                    for (int k = 0; k < 4; k++)
                    {
                        for (int i = 1; i < 4; i++)
                        {
                            int j1 = 0;
                            for (int j = 0; j < 4; j++)
                            {
                                if (k == j)
                                    continue;
                                temp[k]->m_pData[i - 1][j1++]
                                        = this->m_pData[i][j];
                            }
                        }
                    }
                    double det = this->m_pData[0][0] * temp[0]->Determinant()
                            - this->m_pData[0][1] * temp[1]->Determinant()
                            + this->m_pData[0][2] * temp[2]->Determinant()
                            - this->m_pData[0][3] * temp[3]->Determinant();
                    return det;
                }
                    break;
                case 5:
                {
                    Matrix *temp[5];
                    for (int i = 0; i < 5; i++)
                        temp[i] = new Matrix("", 4, 4);
                    for (int k = 0; k < 5; k++)
                    {
                        for (int i = 1; i < 5; i++)
                        {
                            int j1 = 0;
                            for (int j = 0; j < 5; j++)
                            {
                                if (k == j)
                                    continue;
                                temp[k]->m_pData[i - 1][j1++]
                                        = this->m_pData[i][j];
                            }
                        }
                    }
                    double det = this->m_pData[0][0] * temp[0]->Determinant()
                            - this->m_pData[0][1] * temp[1]->Determinant()
                            + this->m_pData[0][2] * temp[2]->Determinant()
                            - this->m_pData[0][3] * temp[3]->Determinant()
                            + this->m_pData[0][4] * temp[4]->Determinant();
                    return det;
                }
                case 1:
                case 2:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                default:
                {
                    int DIM = n;
                    Matrix **temp = new Matrix*[DIM];
                    for (int i = 0; i < DIM; i++)
                        temp[i] = new Matrix("", DIM - 1, DIM - 1);
                    for (int k = 0; k < DIM; k++)
                    {
                        for (int i = 1; i < DIM; i++)
                        {
                            int j1 = 0;
                            for (int j = 0; j < DIM; j++)
                            {
                                if (k == j)
                                    continue;
                                temp[k]->m_pData[i - 1][j1++]
                                        = this->m_pData[i][j];
                            }
                        }
                    }
                    double det = 0;
                    for (int k = 0; k < DIM; k++)
                    {
                        if ((k % 2) == 0)
                            det = det + (this->m_pData[0][k]
                                    * temp[k]->Determinant());
                        else
                            det = det - (this->m_pData[0][k]
                                    * temp[k]->Determinant());
                    }
                    for (int i = 0; i < DIM; i++)
                        delete temp[i];
                    delete[] temp;
                    return det;
                }
                    break;
            }
        }*/
        Matrix& operator =(const Matrix &other)
        {
            if (this->n != other.n || this->n != other.n)
            {
                std::cout
                        << "WARNING: Assignment is taking place with by changing the number of rows and columns of the matrix";
            }
            for (int i = 0; i < n; i++)
                delete[] m_pData[i];
            delete[] m_pData;
            n = n = 0;
            strcpy(m_name, other.m_name);
            n = other.n;
            n = other.n;
            m_pData = new double*[n];
            for (int i = 0; i < n; i++)
                m_pData[i] = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    m_pData[i][j] = other.m_pData[i][j];
                }
            }
            return *this;
        }
        Matrix CoFactor()
        {
            Matrix cofactor("COF", n, n);
            if (n != n)
                return cofactor;
            if (n < 2)
                return cofactor;
            else if (n == 2)
            {
                cofactor.m_pData[0][0] = m_pData[1][1];
                cofactor.m_pData[0][1] = -m_pData[1][0];
                cofactor.m_pData[1][0] = -m_pData[0][1];
                cofactor.m_pData[1][1] = m_pData[0][0];
                return cofactor;
            }
            else if (n >= 3)
            {
                int DIM = n;
                Matrix ***temp = new Matrix**[DIM];
                for (int i = 0; i < DIM; i++)
                    temp[i] = new Matrix*[DIM];
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        temp[i][j] = new Matrix("", DIM - 1, DIM - 1);
                for (int k1 = 0; k1 < DIM; k1++)
                {
                    for (int k2 = 0; k2 < DIM; k2++)
                    {
                        int i1 = 0;
                        for (int i = 0; i < DIM; i++)
                        {
                            int j1 = 0;
                            for (int j = 0; j < DIM; j++)
                            {
                                if (k1 == i || k2 == j)
                                    continue;
                                temp[k1][k2]->m_pData[i1][j1++]
                                        = this->m_pData[i][j];
                            }
                            if (k1 != i)
                                i1++;
                        }
                    }
                }
                bool flagPositive = true;
                for (int k1 = 0; k1 < DIM; k1++)
                {
                    flagPositive = ((k1 % 2) == 0);
                    for (int k2 = 0; k2 < DIM; k2++)
                    {
                        if (flagPositive == true)
                        {
                            cofactor.m_pData[k1][k2]
                                    = temp[k1][k2]->Determinant();
                            flagPositive = false;
                        }
                        else
                        {
                            cofactor.m_pData[k1][k2]
                                    = -temp[k1][k2]->Determinant();
                            flagPositive = true;
                        }
                    }
                }
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        delete temp[i][j];
                for (int i = 0; i < DIM; i++)
                    delete[] temp[i];
                delete[] temp;
            }
            return cofactor;
        }
        Matrix Adjoint()
        {
            Matrix cofactor("COF", n, n);
            Matrix adj("ADJ", n, n);
            if (n != n)
                return adj;
            cofactor = this->CoFactor();
            // adjoint is transpose of a cofactor of a matrix
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    adj.m_pData[j][i] = cofactor.m_pData[i][j];
                }
            }
            return adj;
        }
        Matrix Transpose()
        {
            Matrix trans("TR", n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    trans.m_pData[j][i] = m_pData[i][j];
                }
            }
            return trans;
        }
        Matrix Inverse()
        {
            Matrix cofactor("COF", n, n);
            Matrix inv("INV", n, n);
            if (n != n)
                return inv;
            // to find out Determinant
            double det = determinant();
            std::cout<<det<<std::endl;
            if (fabs(det)<ZERO){
                std::cout<<"Non invertible Matrix"<<std::endl;
            }
            else{
            cofactor = this->CoFactor();
            // inv = transpose of cofactor / Determinant
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    inv.m_pData[j][i] = cofactor.m_pData[i][j]/ det;
                }
            }
            return inv;}
        }
        Matrix operator +(const Matrix &other)
        {
            if (this->n != other.n || this->n != other.n)
            {
                std::cout
                        << "Addition could not take place because number of rows and columns are different between the two matrices";
                return *this;
            }
            Matrix result("", n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    result.m_pData[i][j] = this->m_pData[i][j]
                            + other.m_pData[i][j];
                }
            }
            return result;
        }
        Matrix operator -(const Matrix &other)
        {
            if (this->n != other.n || this->n != other.n)
            {
                std::cout
                        << "Subtraction could not take place because number of rows and columns are different between the two matrices";
                return *this;
            }
            Matrix result("", n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    result.m_pData[i][j] = this->m_pData[i][j]
                            - other.m_pData[i][j];
                }
            }
            return result;
        }
        Matrix operator *(const Matrix &other)
        {
            if (this->n != other.n)
            {
                std::cout
                        << "Multiplication could not take place because number of columns of 1st Matrix and number of rows in 2nd Matrix are different";
                return *this;
            }
            Matrix result("", this->n, other.n);
            for (int i = 0; i < this->n; i++)
            {
                for (int j = 0; j < other.n; j++)
                {
                    for (int k = 0; k < this->n; k++)
                    {
                        result.m_pData[i][j] += this->m_pData[i][k]
                                * other.m_pData[k][j];
                    }
                }
            }
            return result;
        }
        bool operator ==(const Matrix &other)
        {
            if (this->n != other.n || this->n != other.n)
            {
                std::cout
                        << "Comparision could not take place because number of rows and columns are different between the two matrices";
                return false;
            }
            Matrix result("", n, n);
            bool bEqual = true;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (this->m_pData[i][j] != other.m_pData[i][j])
                        bEqual = false;
                }
            }
            return bEqual;
        }
        friend std::istream& operator >>(std::istream &is, Matrix &m);
        friend std::ostream& operator <<(std::ostream &os, const Matrix &m);
};
std::istream& operator >>(std::istream &is, Matrix &m)
{
    std::cout << "\n\nEnter Input For Matrix : " << m.m_name << " Rows: "
            << m.n << " Cols: " << m.n << "\n";
    for (int i = 0; i < m.n; i++)
    {
        for (int j = 0; j < m.n; j++)
        {
            std::cout << "Input For Row: " << i + 1 << " Col: " << j + 1
                    << " = ";
            is >> m.m_pData[i][j];
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    return is;
}
std::ostream& operator <<(std::ostream &os, const Matrix &m)
{
 
    os << "\n\nMatrix : " << m.m_name << " Rows: " << m.n << " Cols: "
            << m.n << "\n\n";
    for (int i = 0; i < m.n; i++)
    {
        os << " | ";
        for (int j = 0; j < m.n; j++)
        {
            char buf[32];
            double data = m.m_pData[i][j];
            if (m.m_pData[i][j] > -0.00001 && m.m_pData[i][j] < 0.00001)
                data = 0;
            sprintf(buf, "%10.2lf ", data);
            os << buf;
        }
        os << "|\n";
    }
    os << "\n\n";
    return os;

#endif

private:


    double computeDeterminant(const vector<vector<double>>& matrix, int size) {
        if (size == 1) {
            return matrix[0][0];
        }

        if (size == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }

        double det = 0;
        vector<vector<double>> subMatrix(size - 1, vector<double>(size - 1));

        for (int x = 0; x < size; x++) {
            // Create the sub-matrix
            for (int i = 1; i < size; i++) {
                int sub_j = 0;
                for (int j = 0; j < size; j++) {
                    if (j == x) continue;
                    subMatrix[i - 1][sub_j] = matrix[i][j];
                    sub_j++;
                }
            }

            // Recursive call
            det += (x % 2 == 0 ? 1 : -1) * matrix[0][x] * computeDeterminant(subMatrix, size - 1);
        }

        return det;
    }
};



int main()
{
    Matrix a("A", 5, 5);
    //std::cin >> a;
    a.FillSimulatedInput();
    Matrix aadj = a.Inverse();
    std::cout << a;
    std::cout << aadj;
    Matrix unit = (a * aadj);
    unit.SetName("A * A-Inv");
    std::cout << unit;
}