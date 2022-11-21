#include "mh.h"
#include <iostream>
#import <cmath>

double** dyn_array(size_t rows, size_t columns){
    double** C=new double*[rows];
    for(size_t i = 0; i < rows; i++){
        C[i] = new double[columns];
    }
    return C;
}

void dyn_array_destroy(double** M,size_t rows){
    for (size_t i = 0; i < rows; ++i)
        delete[] M[i];
    delete[] M;
}

void print_matrix(double **M, size_t rows, size_t columns){
    for(size_t i = 0; i < rows; ++i){
        for(size_t j = 0; j < columns; ++j)
            std::cout << M[i][j] << "\t";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

double** multMatrix(double lambda, double** M, size_t rows, size_t columns){
    auto A = dyn_array(rows,columns);
    for(size_t i = 0; i < rows; ++i)
        for(size_t j = 0; j < columns; ++j)
            A[i][j] = lambda * M[i][j];
    return A;
}

double** plusMatrix(double** A, double** B, size_t rows, size_t columns){
    auto C= dyn_array(rows,columns);
    for(size_t i = 0; i < rows; ++i)
        for(size_t j = 0; j < columns; ++j)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

double** minusMatrix(double** A, double** B, size_t rows, size_t columns){
    auto C= dyn_array(rows,columns);
    for(size_t i = 0; i < rows; ++i)
        for(size_t j = 0; j < columns; ++j)
            C[i][j] = A[i][j] - B[i][j];
    return C;

}

double** multMatrix(double** A, double** B, size_t rowsA, size_t columnsA, size_t columnsB){
    auto C = dyn_array(rowsA,columnsB);
    for(size_t i = 0; i < rowsA; ++i)
        for(size_t j = 0; j < columnsB; ++j)
            C[i][j] = 0;

    for(size_t i = 0; i < rowsA; ++i)
        for(size_t j = 0; j < columnsB; ++j)
            for(size_t k = 0; k < columnsA; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

double& maximumMatrix(double** M, size_t rows, size_t columns){
    size_t max_i=0,max_j=0;
    double max = M[0][0];
    for(size_t i = 0; i < rows; ++i)
        for(size_t j = 0; j < columns; ++j)
            if (M[i][j] > max){
                max=M[i][j];
                max_i=i;
                max_j=j;
            }

    return M[max_i][max_j];
}

double& minmaxMatrix(double** M, size_t rows, size_t columns){
    double min[rows];
    size_t min_j[rows];
    for (size_t i = 0; i < rows; ++i){
        min[i] = M[i][0];
        min_j[i] = 0;
        for(size_t j = 1; j < columns; ++j)
            if (M[i][j] < min[i]){
                min[i] = M[i][j];
                min_j[rows] = j;
            }
    }

    double max = min[0];
    size_t max_j;
    size_t max_i;

    for (size_t i = 1; i < rows; ++i)
        if (max < min[i]) {
            max = min[i];
            max_i = i;
            max_j = min_j[i];
        }

    return M[max_i][max_j];
}

void mixMatrix(double** M, size_t rows, size_t columns, size_t K_1, size_t K_2)
{
    double temp;
    for(size_t i = 0; i < columns; ++i){
        temp = M[K_1][i];
        M[K_1][i] = M[K_2][i];
        M[K_2][i] = temp;
    }
}

void localMinimum(double** M, size_t rows, size_t columns){
    bool A[rows][columns];

    for(size_t i = 0; i < rows; ++i)
        for(size_t j = 0; j < columns; ++j){
            A[i][j] = true;
            for (size_t a = -1; a < 2 && A[i][j]; ++a)
                for(size_t b = -1; b < 2 && A[i][j];++b)
                    if ( i + a >= 0 && i + a < rows && j + b >= 0 && j + b < columns && !(a == 0 && b == 0))
                        if(M[i][j] >= M[i+a][j+b])
                            A[i][j]=false;
        }

    for(size_t i = 0; i < rows; ++i)
        for(size_t j = 0; j < columns; ++j)
            if (A[i][j])
                M[i][j]=0;

}

double det(double** M, size_t size){
    if (size == 2)
        return M[0][0] * M[1][1] - M[0][1] * M[1][0];

    double d = 0;
    //поищем нулевые строки и столбцы
    bool zero = true;
    for(size_t i = 0;i < size; ++i){
        zero = true;
        for (size_t j = 0; j < size && zero; ++j)
            if (M[i][j] != 0)
                zero = false;
        if (zero) return 0;
    }

    for(size_t i = 0;i < size; ++i){
        zero = true;
        for (size_t j = 0; j < size && zero; ++j)
            if (M[j][i] != 0)
                zero = false;
        if (zero) return 0;
    }

    for(size_t n = 0; n < size; ++n){
        if (M[0][n] != 0){
            //построение минора
            auto minor = dyn_array(size - 1, size - 1);
            for (size_t i = 0; i < size - 1; ++i)
                for (size_t j = 0; j < size - 1; ++j){
                    if (n == 0)
                        minor[i][j]=M[i+1][j+1];
                    else if(n == size-1)
                        minor[i][j]=M[i+1][j];
                    else if (j < n)
                        minor[i][j] = M[i + 1][j];
                    else if (j >= n)
                        minor[i][j] = M[i + 1][j + 1];
                }
            d += pow(-1, 1 + n) * M[0][n] * det(minor, size - 1);
            dyn_array_destroy(minor,size-1);
        }
    }

    return d;
}

double det_gauss(double ** M, size_t size){
    for(size_t n = 0; n < size; ++n) {
        if (M[n][n] == 0){
            bool zero = true;
            for (int i = n + 1; M[i][n] != 0 && i < size; ++i)
                if (M[i][n] != 0) {
                    zero = false;
                    mixMatrix(M, size, size, n, i);
                }
            if (zero) return 0;
        }

        //ведущая строка
        double row[size - n];
        for (size_t i = 0; i < size - n; ++i)
            row[i] = M[n][n + i] / M[n][n];

        for (size_t i = n + 1; i < size; ++i){
            double multiple = M[i][n];
            for (size_t j = n; j < size; ++j)
                M[i][j] += (-1) * multiple * row[j - n];
        }
    }

    double det = 1;
    for(size_t i = 0; i < size; ++i)
        det = det * M[i][i];

    return det;
}

double** inv(double** M, size_t size){
    if (det(M,size) == 0){
        std::cout << "Non invertible" << std::endl;
        return M;
    }

    //заполняем единичную матрицу для метода гаусса
    auto inv= dyn_array(size,size);
    for (size_t i = 0; i < size; ++i)
        for(size_t j = 0; j < size; ++j)
            if (i == j)
                inv[i][j]=1;
            else
                inv[i][j]=0;

    //инициализация копии матрицы
    auto M_copy = dyn_array(size,size);
    for(size_t i = 0; i < size; ++i)
        for(size_t j = 0; j < size; ++j)
            M_copy[i][j]=M[i][j];

    //выводит
    for (size_t i = 0; i < size; ++i){
        for(size_t j = 0; j < size; ++j)
            std::cout << M_copy[i][j] <<"\t";
        std::cout << "|";
        for(size_t j = 0; j < size; ++j)
            std::cout << inv[i][j] <<"\t";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    for(size_t n = 0; n < size; ++n){
        double multiple = M_copy[n][n];
        for(size_t j = 0; j < size; ++j){
            if (j>=n)
                M_copy[n][j] = M_copy[n][j] / multiple;
            inv[n][j] = inv[n][j] / multiple;
        }

        //инициализация копии ведущей строки
        double row[size - n];
        double row_inv[size];
        for (size_t i = 0; i < size - n; ++i)
            row[i] = M_copy[n][n + i];
        for(size_t i=0; i < size; ++i)
            row_inv[i] = inv[n][i];

        for (size_t i = 0; i < size; ++i){
            double multiple = M_copy[i][n];
            for (size_t j = 0; j < size; ++j)
                if (i != n)
                {
                    if (j >= n)
                        M_copy[i][j] += (-1) * multiple * row[j - n];
                    inv[i][j] += (-1) * multiple * row_inv[j];}
        }
    }
    return inv;
}

double matrix_norm(double **M, size_t n, bool type) { //1 = 1 inf = 0
    double max=0; double sum;
        for (size_t i=0;i<n;++i) {
            sum = 0;
            for (size_t j = 0; j < n; ++j)
                sum += type*M[j][i] + (!type)*M[i][j];
            if (sum > max) max = sum;
        }
    return max;
}

double vector_norm(double**x1,double **x2,size_t size) {
    double norm=0;
    for(size_t i=0;i<size;++i)
        norm += (x1[i][0]-x2[i][0])*(x1[i][0]-x2[i][0]);
    norm = sqrt(norm);
    return norm;
}


double** E_min_muA(double mu, double **A, size_t n){
    auto B=dyn_array(n,n);
    for(size_t i=0; i<n;++i)
        for(size_t j=0;j<n;++j)
            B[i][j] = (i==j) - mu*A[i][j];
    return B;
}

void Bx_plus_с(double** B, double **x1, double **x2, double **c, size_t n) {
    for(size_t i=0; i<n; ++i){
            x1[i][0]=x2[i][0];
            x2[i][0]=c[i][0];}
    for(size_t i=0;i<n;++i)
        for(size_t j=0;j<n;++j)
            x2[i][0] += B[i][j]*x1[j][0];
    print_matrix(x2,n,1);
}

bool matrix_pos_def(double **A, size_t n){
    for(size_t i=1; i<n+1; ++i)
        if (det(A,i)<0)
            return false;
    return true;

}

void A_t_A(double**A, size_t n){
    auto *B=dyn_array(n,n);
    for(size_t i=0; i<n; ++i)
        for(size_t j=0;j<n;++j)
            B[i][j]=A[i][j]*A[j][i];
    for(size_t i=0; i<n; ++i)
        for(size_t j=0;j<n;++j)
            A[i][j]=B[i][j];
    dyn_array_destroy(B,n);
}

