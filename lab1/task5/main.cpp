#include "main.h"

#include <format>

int main(int argc, char **argv) {
    if (argc == 2) {
        Matrix<double> A = read_matrix_from_file<double>(argv[1]);

        auto QR = QR_Decomposition(A);

        std::cout << "Q Matrix:\n" << QR.first << "\nR Matrix:\n" << QR.second << std::endl;

        QR_eigenvalues(A, std::numeric_limits<double>::epsilon());
    } else {
        throw std::invalid_argument("Incorrect arguments. Add path to your test file");
    }

    return 0;
}

template<class T>
Matrix<T> read_matrix_from_file(const char *path) {
    std::ifstream stream(path);
    if (!stream.is_open()) {
        throw std::runtime_error("Cannot open file");
    }

    int n;
    stream >> n;

    if (n <= 0) {
        throw std::domain_error("Matrix size cannot be less than 1");
    }

    Matrix<T> matrix(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (!(stream >> matrix.coefficients[i][j])) {
                throw std::runtime_error("File format is not correct");
            }
        }
    }

    return matrix;
}

template<class T>
T Euclid_norm_root(const Matrix<T> &A, int start_index = 0) {
    T result = 0;

    for (int i = start_index; i < A.n; ++i) {
        result += A.coefficients[i][start_index] * A.coefficients[i][start_index];
    }

    return std::sqrt(result);
}

template<class T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<class T>
std::pair<Matrix<T>, Matrix<T> > QR_Decomposition(Matrix<T> &A) {
    const int n = A.n;

    Matrix<T> Q = Matrix<T>::Identity(n), R = A;

    for (int i = 0; i < n - 1; ++i) {
        Matrix<T> v(n, 1);

        v.coefficients[i][0] = R.coefficients[i][i] + sgn(R.coefficients[i][i]) * Euclid_norm_root(R, i);

        for (int j = i + 1; j < n; ++j) {
            v.coefficients[j][0] = R.coefficients[j][i];
        }

        Matrix<T> H = Matrix<T>::Householder(v);

        Q *= H;
        R = H * R;
    }

    return {Q, R};
}

std::pair<std::complex<double>, std::complex<double> > solve_quadratic(double a, double b, double c, double d) {
    const double D = a * a - 2 * a * d + 4 * b * c + d * d;

    std::complex<double> x1, x2;

    if (D > 0) {
        const double D_sqrt = std::sqrt(D);
        x1 = 0.5 * (D_sqrt + a + d);
        x2 = 0.5 * (-D_sqrt + a + d);
    } else if (D == 0) {
        x1 = 0.5 * (a + d);
        x2 = x1;
    } else {
        const double D_sqrt = std::sqrt(std::abs(D));

        x1.real(0.5 * (a + d));
        x1.imag(0.5 * D_sqrt);
        x2.real(0.5 * (a + d));
        x2.imag(-0.5 * D_sqrt);
    }

    return {x1, x2};
}

Matrix<double> QR_eigenvalues(const Matrix<double> &A, const double eps) {
    double norm_value = 1;

    Matrix<double> A_k = A;
    Matrix<double> A_k_prev = A_k;

    int iter = 0;

    while (true) {
        ++iter;
        auto QR = QR_Decomposition(A_k);

        A_k = QR.second * QR.first;

        bool converged = true;
        for (int i = 1; i < A.n; ++i) {
            if (std::abs(A_k.coefficients[i][i - 1]) > eps) {
                converged = false;
                break;
            }
        }
        if (converged) {
            for (int i = 0; i < A.n; ++i) {
                std::cout << std::format("Root #{}: {:.6f}", i + 1, A_k.coefficients[i][i]) << std::endl;
            }

            break;
        }

        for (int i = 0; i < A.n - 1; ++i) {
            double **a = A_k.coefficients;
            auto [x1, x2] = solve_quadratic(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1]);
            if (x1.imag() != 0) {
                a = A_k_prev.coefficients;
                auto [x1_prev, x2_prev] = solve_quadratic(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1]);

                if (std::abs(x1 - x1_prev) <= eps) {
                    int real_pos = (i + 2) % A.n;

                    std::cout << std::format("Complex roots: {:.6f} + {:.6f}i; {:.6f} - {:.6f}i.", x1.real(), x1.imag(),
                                             x2.real(), -x2.imag()) << std::endl;
                    std::cout << std::format("Real root: {:.6f}", A_k.coefficients[real_pos][real_pos]) << std::endl;
                    converged = true;
                }
            }
        }

        if (converged)
            break;

        A_k_prev = A_k;
    }

    std::cout << "Total iterations: " << iter << std::endl;

    return A_k;
}
