#ifndef LUP_H
#define LUP_H

#include <matrix.h>

template <class T>
struct LUPResult {
    Matrix<T> L;
    Matrix<T> U;
    Matrix<T> P;
    int swaps;

    explicit LUPResult(int n) : L(n, n), U(n, n), P(n, n), swaps(0) {}
};

template<class T>
LUPResult<T> LUPDecomposition(const Matrix<T> &A) {
    int n = A.n;
    LUPResult<T> result(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result.P.coefficients[i][j] = (i == j) ? 1 : 0;
        }
    }

    Matrix<T> PA = A;

    for (int k = 0; k < n; k++) {
        T maxVal = 0;
        int maxIndex = k;
        for (int i = k; i < n; i++) {
            T absVal = std::abs(PA.coefficients[i][k]);
            if (absVal > maxVal) {
                maxVal = absVal;
                maxIndex = i;
            }
        }

        if (maxVal == 0)
            throw std::runtime_error("Matrix is singular");

        if (maxIndex != k) {
            for (int j = 0; j < n; j++) {
                std::swap(PA.coefficients[k][j], PA.coefficients[maxIndex][j]);
                std::swap(result.P.coefficients[k][j], result.P.coefficients[maxIndex][j]);
            }
            ++result.swaps;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            T sum = 0;
            for (int k = 0; k < i; k++) {
                sum += result.L.coefficients[i][k] * result.U.coefficients[k][j];
            }
            result.U.coefficients[i][j] = PA.coefficients[i][j] - sum;
        }

        for (int j = i; j < n; j++) {
            if (i == j) {
                result.L.coefficients[i][i] = 1;
            } else {
                T sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += result.L.coefficients[j][k] * result.U.coefficients[k][i];
                }
                result.L.coefficients[j][i] = (PA.coefficients[j][i] - sum) / result.U.coefficients[i][i];
            }
        }
    }

    return result;
}

template<class T>
T *solveLUP(const LUPResult<T> &lup, const T *b) {
    const int n = lup.L.n;
    T *y = new T[n];
    T *x = new T[n];

    T *b_perm = new T[n];
    for (int i = 0; i < n; i++) {
        b_perm[i] = 0;
        for (int j = 0; j < n; j++) {
            b_perm[i] += lup.P.coefficients[i][j] * b[j];
        }
    }

    for (int i = 0; i < n; i++) {
        T sum = 0;
        for (int j = 0; j < i; j++) {
            sum += lup.L.coefficients[i][j] * y[j];
        }
        y[i] = b_perm[i] - sum;
    }

    for (int i = n - 1; i >= 0; --i) {
        T sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += lup.U.coefficients[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / lup.U.coefficients[i][i];
    }

    delete[] y;
    delete[] b_perm;
    return x;
}

template<class T>
T determinantLUP(const LUPResult<T> &lup) {
    T det = 1;
    for (int i = 0; i < lup.U.n; i++) {
        det *= lup.U.coefficients[i][i];
    }
    return (lup.swaps % 2 == 0) ? det : -det;
}


template<class T>
T determinant(const Matrix<T>& A) {
    if (A.n == A.m && A.n == 2) {
        T** coef = A.coefficients;

        return coef[0][0] * coef[1][1] - coef[1][0] * coef[0][1];
    }

    LUPResult<T> lup = LUPDecomposition(A);

    return determinantLUP(lup);
}

template<class T>
Matrix<T> inverseMatrix(const Matrix<T> &A) {
    int n = A.n;

    if (n != A.m)
        throw std::runtime_error("Matrix must be square for inversion");

    LUPResult<T> lup = LUPDecomposition(A);

    Matrix<T> inv(n, n);

    for (int i = n - 1; i >= 0; --i) {
        for (int j = n - 1; j > i - 1; --j) {
            T row_sum = 0;
            for (int k = i + 1; k < n; ++k) {
                row_sum += lup.U.coefficients[i][k] * inv.coefficients[k][j];
            }

            if (i == j)
                inv.coefficients[i][i] = (1.0 - row_sum) / lup.U.coefficients[i][i];
            else
                inv.coefficients[i][j] = -row_sum / lup.U.coefficients[i][i];
        }

        for (int j = i - 1; j >= 0; --j) {
            T row_sum = 0;
            for (int k = j + 1; k < n; ++k) {
                row_sum += inv.coefficients[i][k] * lup.L.coefficients[k][j];
            }
            inv.coefficients[i][j] = -row_sum;
        }
    }

    return inv * lup.P;
}
#endif // LUP_H
