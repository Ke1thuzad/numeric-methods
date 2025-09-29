#include "main.h"


int main(const int argc, char **argv) {
    if (argc == 2) {
        const LinearEquation<double> equation = read_equation_from_file(argv[1]);

        std::cout << "Original matrix A:" << std::endl;
        std::cout << equation.matrix;

        std::cout << "Vector b:" << std::endl;
        for (int i = 0; i < equation.n; i++) {
            std::cout << equation.b[i] << " ";
        }
        std::cout << std::endl << std::endl;

        const LUPResult<double> lup = LUPDecomposition(equation.matrix);

        const double* solution = solveLUP(lup, equation.b);

        std::cout << "Solution x:" << std::endl;
        for (int i = 0; i < equation.n; i++) {
            std::cout << "x[" << i << "] = " << solution[i] << std::endl;
        }
        std::cout << std::endl;

        // Calculate determinant
        const double det = determinantLUP(lup);
        std::cout << "Determinant of A: " << det << std::endl << std::endl;

        const Matrix<double> inv = inverseMatrix(equation.matrix);
        std::cout << "Inverse matrix A^(-1):" << std::endl;
        std::cout << inv;

        std::cout << "Inverse check: " << std::endl;
        std::cout << equation.matrix * inv;
    } else {
        throw std::runtime_error("Incorrect arguments. Add path to your test file");
    }
}

LinearEquation<double> read_equation_from_file(const char *path) {
    std::ifstream stream(path);
    if (!stream.is_open()) {
        throw std::runtime_error("Cannot open file");
    }

    int n;
    stream >> n;

    if (n <= 0) {
        throw std::domain_error("Matrix size cannot be less than 1");
    }

    LinearEquation<double> equation(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (!(stream >> equation.matrix.coefficients[i][j])) {
                throw std::runtime_error("File format is not correct");
            }
        }
    }

    for (int i = 0; i < n; i++) {
        if (!(stream >> equation.b[i])) {
            throw std::runtime_error("File format is not correct");
        }
    }

    return equation;
}

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

        if (maxVal == 0) {
            throw std::runtime_error("Matrix is singular");
        }

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
Matrix<T> inverseMatrix(const Matrix<T> &A) {
    int n = A.n;

    if (n != A.m)
        throw std::runtime_error("Matrix must be square for inversion");

    LUPResult<T> lup = LUPDecomposition(A);

    Matrix<T> Z(n, n);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            T sum = 0;
            for (int k = 0; k < i; k++) {
                sum += lup.L.coefficients[i][k] * Z.coefficients[k][j];
            }
            Z.coefficients[i][j] = lup.P.coefficients[i][j] - sum;
        }
    }

    Matrix<T> U_inv(n, n);

    for (int j = n - 1; j >= 0; --j) {
        T sum_diag = 0;
        for (int k = j + 1; k < n; ++k) {
            sum_diag += lup.U.coefficients[j][k] * U_inv.coefficients[k][j];
        }
        U_inv.coefficients[j][j] = (1.0 - sum_diag) / lup.U.coefficients[j][j];

        for (int i = j - 1; i >= 0; --i) {
            T sum_upper = 0;
            for (int k = i + 1; k < n; ++k) {
                sum_upper += lup.U.coefficients[i][k] * U_inv.coefficients[k][j];
            }
            U_inv.coefficients[i][j] = -sum_upper / lup.U.coefficients[i][i];
        }

        for (int i = n - 1; i > j; --i) {
            T sum_lower = 0;
            for (int k = j + 1; k < n; ++k) {
                sum_lower += U_inv.coefficients[i][k] * lup.U.coefficients[k][j];
            }
            U_inv.coefficients[i][j] = -sum_lower;
        }
    }

    return U_inv * Z;
}
