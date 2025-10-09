#include "main.h"

int main(const int argc, char **argv) {
    if (argc == 2) {
        auto read = read_tridiagonal_matrix_from_file(argv[1]);
        const Matrix<double> matrix = read.first;
        double *d = read.second;

        const double* res = thomas_algorithm(matrix, d);



        for (int i = 0; i < matrix.n; ++i) {
            std::cout << res[i] << ", ";
        }

        std::cout << std::endl;

        delete[] res;
        delete[] d;

    } else if (argc > 2) {
        std::cout << "Invalid arguments. Usage: ./main <matrix filepath>" << std::endl;
    }

    return 0;
}

std::pair<Matrix<double>, double*> read_tridiagonal_matrix_from_file(const char *path) {
    std::ifstream stream(path);
    if (!stream.good())
        throw std::runtime_error("Stream has not been opened");

    int n, read = 0;
    stream >> n;

    Matrix<double> matrix(n, n);

    for (int i = 0; i < n && !stream.eof(); ++i) {
        for (int j = 0; j < 3 && !stream.eof(); ++j) {
            if (i - 1 + j < 0 || i - 1 + j >= n) {
                continue;
            }
            stream >> matrix.coefficients[i][i - 1 + j];
            ++read;
        }
    }

    std::cout << matrix << std::endl;

    if (n > 3 && read != (n - 2) * 3 + 4)
        throw std::runtime_error("File format is not correct");

    read = 0;

    double *d = new double[n];

    for (int i = 0; i < n && !stream.eof(); ++i) {
        stream >> d[i];
        ++read;
    }

    if (read != n)
        throw std::runtime_error("File format is not correct");

    return {matrix, d};
}

double* thomas_algorithm(const Matrix<double>& mat, const double* d) {
    double* P = new double[mat.n];
    double* Q = new double[mat.n];
    double* x = new double[mat.n];

    thomas_propagation(mat, P, Q, d);

    x[mat.n - 1] = Q[mat.n - 1];

    for (int i = mat.n - 2; i >= 0; --i) {
        x[i] = P[i] * x[i + 1] + Q[i];
    }

    delete[] P;
    delete[] Q;

    return x;
}

void thomas_propagation(const Matrix<double>& mat, double* P, double* Q, const double* d) {
    for (int i = 0; i < mat.n; ++i) {
        const double b = mat.coefficients[i][i];
        const double c = i + 1 != mat.n ? mat.coefficients[i][i + 1] : 0;

        if (i == 0) {
            P[i] = -c / b;
            Q[i] = d[i] / b;
        } else {
            const double a = mat.coefficients[i][i - 1];

            P[i] = (-c) / (b + a * P[i - 1]);
            Q[i] = (d[i] - a * Q[i - 1]) / (b + a * P[i - 1]);
        }
    }
}

